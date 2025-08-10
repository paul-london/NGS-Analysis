import glob
import re
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
        expand("results/vcf/{sample}.vcf.gz", sample=SAMPLES.keys()),
        expand("results/{sample}.variants.tsv", sample=SAMPLES.keys()),
        "qc/aligned/summary_metrics.csv",
        "qc/fastqc/fastqc_summary.csv",
        f"{REF}.fai",
        f"{REF}.dict",
        f"{REF}.amb",
        f"{REF}.ann",
        f"{REF}.bwt",
        f"{REF}.pac",
        f"{REF}.sa"

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

# FastQC

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
        import pandas as pd
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
        bam="output/{sample}.bam",
        bai="output/{sample}.bam.bai"
    threads: 8
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.R1} {input.R2} | \
        samtools view -bS - | \
        samtools sort -o {output.bam}
        samtools index {output.bam}
        """

# Alignment QC

rule samtools_flagstat:
    input:
        bam="output/{sample}.bam"
    output:
        flagstat="qc/aligned/{sample}.flagstat.txt"
    shell:
        """
        samtools flagstat {input.bam} > {output.flagstat}
        """

rule samtools_stats:
    input:
        bam="output/{sample}.bam"
    output:
        stats="qc/aligned/{sample}.stats.txt"
    shell:
        """
        samtools stats {input.bam} > {output.stats}
        """

rule aggregate_flagstat:
    input:
        expand("qc/aligned/{sample}.flagstat.txt", sample=SAMPLES.keys())
    output:
        "qc/aligned/flagstat_summary.csv"
    run:
        import re
        import csv

        summary = []
        for f in input:
            sample = f.split("/")[-1].replace(".flagstat.txt", "")
            with open(f) as fh:
                lines = fh.readlines()
            # Parse total reads and mapped reads from flagstat output
            total_reads = int(re.match(r"(\d+) \+ \d+ in total", lines[0]).group(1))
            mapped_reads = int(re.match(r"(\d+) \+ \d+ mapped", lines[4]).group(1))
            pct_mapped = mapped_reads / total_reads * 100 if total_reads else 0
            summary.append({"sample": sample, "total_reads": total_reads, "mapped_reads": mapped_reads, "pct_mapped": pct_mapped})

        keys = summary[0].keys()
        with open(output[0], "w", newline="") as csvfile:
            dict_writer = csv.DictWriter(csvfile, keys)
            dict_writer.writeheader()
            dict_writer.writerows(summary)

# Variant calling

rule call_variants:
    input:
        bam="output/{sample}.bam",
        bai="output/{sample}.bam.bai",
        ref=REF,
        bed=lambda w: SAMPLES[w.sample]["bed"]
    output:
        vcf="results/vcf/{sample}.vcf.gz"
    shell:
        """
        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -L {input.bed} \
            -O {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule extract_variants:
    input:
        vcf="results/vcf/{sample}.vcf.gz"
    output:
        tsv="results/{sample}.variants.tsv"
    shell:
        """
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\n' {input.vcf} > {output.tsv}
        """