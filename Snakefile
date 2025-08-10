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
        expand("results/{sample}.vcf.gz", sample=SAMPLES.keys())

rule index_reference:
    input:
        fasta=REF
    output:
        fasta_fai="reference/genome/hg38.fai",
        dict="reference/genome/hg38.dict"
    shell:
        """
        samtools faidx {input.fasta}
        gatk CreateSequenceDictionary -R {input.fasta} -O {output.dict}
        """

rule bwa_index:
    input:
        fasta=REF
    output:
        amb=REF + ".amb",
        ann=REF + ".ann",
        bwt=REF + ".bwt",
        pac=REF + ".pac",
        sa=REF + ".sa"
    shell:
        "bwa index {input.fasta}"

rule align_reads:
    input:
        ref=REF,
        ref_amb="reference/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.amb",
        ref_ann="reference/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.ann",
        ref_bwt="reference/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.bwt",
        ref_pac="reference/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.pac",
        ref_sa="reference/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.sa",
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

rule call_variants:
    input:
        bam="output/{sample}.bam",
        bai="output/{sample}.bam.bai",
        ref=REF,
        bed=lambda w: SAMPLES[w.sample]["bed"]
    output:
        vcf="results/{sample}.vcf.gz"
    shell:
        """
        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -L {input.bed} \
            -O {output.vcf}
        tabix -p vcf {output.vcf}
        """
