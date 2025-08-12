import glob
import re
import csv
import pandas as pd
import collections
import collections.abc
collections.Iterable = collections.abc.Iterable

REF = "reference/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# Pattern: <sampleID>_<panel>_<seqIndex>_R1_001.fastq.gz
R1_files = glob.glob("data/raw_reads/*_R1_001.fastq.gz")

SAMPLES = {}
for f in R1_files:
    basename = f.split("/")[-1]
    match = re.match(r"(F\d+)_(.+)_(S\d+)_R1_001\.fastq\.gz", basename)
    if not match:
        raise ValueError(f"Filename {basename} does not match expected pattern F<id>_<panel>_S##_R1_001.fastq.gz")
    sample_id, panel, seq_index = match.groups()
    sample_name = f"{sample_id}_{panel}"  # logical sample name in outputs
    SAMPLES[sample_name] = {
        "panel": panel,
        "R1": f,
        "R2": f.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz"),
        "bed": f"reference/panel/{panel}.bed"
    }

rule all:
    input:
        expand("data/variant_calls/{sample}.annotated.tsv", sample=SAMPLES.keys()),
        expand("results/{sample}.annotated.vcf.gz", sample=SAMPLES.keys()),
        expand("qc/aligned_post_bqsr/{sample}.flagstat.txt", sample=SAMPLES.keys()),
        expand("qc/aligned_post_bqsr/{sample}.stats.txt", sample=SAMPLES.keys()),
        "qc/fastqc/fastqc_summary.csv",
        "qc/aligned/summary_metrics.csv",
        "qc/aligned_post_bqsr/summary_metrics.csv"

# Indexing reference genome

rule index_reference:
    input:
        fasta=REF
    output:
        fasta_fai=f"{REF}.fai",
        dict=f"{REF}.dict"
    shell:
        """
        samtools faidx {input.fasta}
        gatk CreateSequenceDictionary -R {input.fasta} -O {output.dict}
        """

rule bwa_index:
    input:
        fasta=REF
    output:
        amb=f"{REF}.amb",
        ann=f"{REF}.ann",
        bwt=f"{REF}.bwt",
        pac=f"{REF}.pac",
        sa=f"{REF}.sa"
    shell:
        "bwa index {input.fasta}"

# FastQC - REPLACE WITH FASTP

rule fastqc:
    input:
        fq="data/raw_reads/{sample}_{read}_001.fastq.gz"
    output:
        html="qc/fastqc/{sample}_{read}_001_fastqc.html",
        zip="qc/fastqc/{sample}_{read}_001_fastqc.zip"
    threads: 2
    shell:
        """
        fastqc -t {threads} -o qc/fastqc {input.fq}
        """

rule parse_fastqc_summary:
    input:
        summary="qc/fastqc/{sample}_{read}_001_fastqc/summary.txt"
    output:
        csv="qc/fastqc/{sample}_{read}_metrics.csv"
    run:
        with open(input.summary) as f:
            lines = f.readlines()
        with open(output.csv, "w") as out:
            out.write("status,metric\n")
            for line in lines:
                status, metric, _ = line.strip().split("\t")
                out.write(f"{status},{metric}\n")

rule aggregate_fastqc_metrics:
    input:
        expand("qc/fastqc/{sample}_{read}_metrics.csv", sample=SAMPLES.keys(), read=["R1", "R2"])
    output:
        "qc/fastqc/fastqc_summary.csv"
    run:
        dfs = []
        for csvfile in input:
            df = pd.read_csv(csvfile)
            filename = csvfile.split("/")[-1]
            sample = filename.split("_")[0]
            read = filename.split("_")[1]
            df["sample"] = sample
            df["read"] = read
            dfs.append(df)
        combined = pd.concat(dfs)
        combined.to_csv(output[0], index=False)

# Alignment

rule align_reads:
    input:
        ref=REF,
        ref_amb=f"{REF}.amb",
        ref_ann=f"{REF}.ann",
        ref_bwt=f"{REF}.bwt",
        ref_pac=f"{REF}.pac",
        ref_sa=f"{REF}.sa",
        R1=lambda w: SAMPLES[w.sample]["R1"],
        R2=lambda w: SAMPLES[w.sample]["R2"]
    output:
        bam="data/aligned_reads/{sample}.bam",
        bai="data/aligned_reads/{sample}.bam.bai"
    threads: 8
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.R1} {input.R2} | \
        samtools view -bS - | \
        samtools sort -o {output.bam} && \
        samtools index {output.bam}
        """

# Alignment QC

rule samtools_flagstat:
    input:
        bam="data/aligned_reads/{sample}.bam"
    output:
        flagstat="qc/aligned/{sample}.flagstat.txt"
    shell:
        """
        samtools flagstat {input.bam} > {output.flagstat}
        """

rule samtools_stats:
    input:
        bam="data/aligned_reads/{sample}.bam"
    output:
        stats="qc/aligned/{sample}.stats.txt"
    shell:
        """
        samtools stats {input.bam} > {output.stats}
        """

rule aggregate_flagstat:
    input:
        flagstats=expand("qc/aligned/{sample}.flagstat.txt", sample=SAMPLES.keys()),
        depths=expand("qc/aligned/{sample}.depth_summary.csv", sample=SAMPLES.keys())
    output:
        "qc/aligned/summary_metrics.csv"
    run:
        summary = []
        # Read flagstat metrics
        for f in input.flagstats:
            sample = f.split("/")[-1].replace(".flagstat.txt", "")
            with open(f) as fh:
                lines = fh.readlines()
            total_reads = int(re.match(r"(\d+) \+ \d+ in total", lines[0]).group(1))
            mapped_reads = int(re.match(r"(\d+) \+ \d+ mapped", lines[4]).group(1))
            pct_mapped = mapped_reads / total_reads * 100 if total_reads else 0

            # Read depth summary for the sample
            depth_file = f"qc/aligned/{sample}.depth_summary.csv"
            depth_df = pd.read_csv(depth_file)
            mean_depth = depth_df["mean_depth"].values[0]

            summary.append({
                "sample": sample,
                "total_reads": total_reads,
                "mapped_reads": mapped_reads,
                "pct_mapped": pct_mapped,
                "mean_depth": mean_depth
            })

        keys = summary[0].keys()
        with open(output[0], "w", newline="") as csvfile:
            dict_writer = csv.DictWriter(csvfile, keys)
            dict_writer.writeheader()
            dict_writer.writerows(summary)

# Mark duplicates

rule mark_duplicates:
    input:
        bam="data/aligned_reads/{sample}.bam",
        bai="data/aligned_reads/{sample}.bam.bai"
    output:
        dedup_bam="data/deduplicated/{sample}.dedup.bam",
        dedup_bai="data/deduplicated/{sample}.dedup.bam.bai",
        metrics="dedup/{sample}.metrics.txt"
    shell:
        """
        gatk MarkDuplicates \
          -I {input.bam} \
          -O {output.dedup_bam} \
          -M {output.metrics} \
          --CREATE_INDEX true
        """

# Base quality score recalibration

rule base_recalibrator:
    input:
        bam="data/deduplicated/{sample}.dedup.bam",
        bai="data/deduplicated/{sample}.dedup.bam.bai",
        ref=REF,
        known_sites="reference/variants/All_20180418.vcf.bgz"  # e.g., dbSNP or other trusted variants
    output:
        recal_table="data/recalibrated_reads/{sample}.recal.table"
    shell:
        """
        gatk BaseRecalibrator \
          -I {input.bam} \
          -R {input.ref} \
          --known-sites {input.known_sites} \
          -O {output.recal_table}
        """

rule apply_bqsr:
    input:
        bam="data/deduplicated/{sample}.dedup.bam",
        bai="data/deduplicated/{sample}.dedup.bam.bai",
        ref=REF,
        recal_table="data/recalibrated_reads/{sample}.recal.table"
    output:
        bam_recal="data/recalibrated_reads/{sample}.recal.bam",
        bai_recal="data/recalibrated_reads/{sample}.recal.bam.bai"
    shell:
        """
        gatk ApplyBQSR \
          -R {input.ref} \
          -I {input.bam} \
          --bqsr-recal-file {input.recal_table} \
          -O {output.bam_recal} &&
        samtools index {output.bam_recal}
        """

# Post-BQSR QC

rule samtools_flagstat_post_bqsr:
    input:
        bam="data/recalibrated_reads/{sample}.recal.bam"
    output:
        flagstat="qc/aligned_post_bqsr/{sample}.flagstat.txt"
    shell:
        "samtools flagstat {input.bam} > {output.flagstat}"

rule samtools_stats_post_bqsr:
    input:
        bam="data/recalibrated_reads/{sample}.recal.bam"
    output:
        stats="qc/aligned_post_bqsr/{sample}.stats.txt"
    shell:
        "samtools stats {input.bam} > {output.stats}"

rule depth_summary_post_bqsr:
    input:
        bam="data/recalibrated_reads/{sample}.recal.bam"
    output:
        depth_summary="qc/aligned_post_bqsr/{sample}.depth_summary.csv"
    shell:
        """
        samtools depth -a {input.bam} | \
        awk '{{sum+=$3; count++}} END {{print "mean_depth", sum/count}}' > {output.depth_summary}
        """

rule aggregate_flagstat_post_bqsr:
    input:
        flagstats=expand("qc/aligned_post_bqsr/{sample}.flagstat.txt", sample=SAMPLES.keys()),
        depths=expand("qc/aligned_post_bqsr/{sample}.depth_summary.csv", sample=SAMPLES.keys())
    output:
        "qc/aligned_post_bqsr/summary_metrics.csv"
    run:
        summary = []
        for f in input.flagstats:
            sample = f.split("/")[-1].replace(".flagstat.txt", "")
            with open(f) as fh:
                lines = fh.readlines()
            total_reads = int(re.match(r"(\d+) \+ \d+ in total", lines[0]).group(1))
            mapped_reads = int(re.match(r"(\d+) \+ \d+ mapped", lines[4]).group(1))
            pct_mapped = mapped_reads / total_reads * 100 if total_reads else 0

            # Read depth summary for the sample
            depth_file = f"qc/aligned_post_bqsr/{sample}.depth_summary.csv"
            depth_df = pd.read_csv(depth_file)
            mean_depth = depth_df["mean_depth"].values[0]

            summary.append({
                "sample": sample,
                "total_reads": total_reads,
                "mapped_reads": mapped_reads,
                "pct_mapped": pct_mapped,
                "mean_depth": mean_depth
            })

        keys = summary[0].keys()
        with open(output[0], "w", newline="") as csvfile:
            dict_writer = csv.DictWriter(csvfile, keys)
            dict_writer.writeheader()
            dict_writer.writerows(summary)

# Variant calling

rule call_variants:
    input:
        bam="data/recalibrated_reads/{sample}.recal.bam",
        bai="data/recalibrated_reads/{sample}.recal.bam.bai",
        ref=REF,
        bed=lambda w: SAMPLES[w.sample]["bed"]
    output:
        vcf="data/variant_calls/{sample}.vcf.gz"
    shell:
        """
        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -L {input.bed} \
            -O {output.vcf}
        tabix -p vcf {output.vcf}
        """

# Filter variants

rule filter_variants:
    input:
        vcf="data/variant_calls/{sample}.vcf.gz",
        ref=REF
    output:
        vcf="data/variant_calls/{sample}.filtered.vcf.gz"
    shell:
        """
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.vcf} \
            --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
            --filter-name "FAIL" \
            -O {output.vcf}
        tabix -p vcf {output.vcf}
        """

# Annotate variants

rule annotate_variants:
    input:
        vcf="data/variant_calls/{sample}.filtered.vcf.gz"
    output:
        annotated_vcf="results/{sample}.annotated.vcf.gz",
        annotated_vcf_index="results/{sample}.annotated.vcf.gz.tbi"
    params:
        snpeff_db="GRCh38.86"  # adjust to your SnpEff database name
    shell:
        """
        snpEff {params.snpeff_db} {input.vcf} | bgzip -c > {output.annotated_vcf}
        tabix -p vcf {output.annotated_vcf}
        """

rule extract_variants:
    input:
        vcf="results/{sample}.annotated.vcf.gz"
    output:
        tsv="results/{sample}.variants.tsv"
    shell:
        """
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\n' {input.vcf} > {output.tsv}
        """

# Remove intermediate files

rule clean:
    shell:
        "rm -rf data/aligned_reads/* data/deduplicated/* data/recalibrated_reads/*"