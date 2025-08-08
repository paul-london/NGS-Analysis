# Config section: adjust paths and sample names here
SAMPLES = ["sample"]  # list your samples here
REF = "references/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

rule all:
    input:
        expand("results/{sample}/variants/{sample}.vcf.gz", sample=SAMPLES)

rule bwa_mem:
    input:
        r1="data/raw_reads/{sample}_R1.fastq.gz",
        r2="data/raw_reads/{sample}_R2.fastq.gz",
        ref=REF + ".fa" if not REF.endswith(".fa") else REF
    output:
        temp("results/{sample}/aligned/{sample}.sam")
    threads: 4
    shell:
        """
        bwa mem -M -t {threads} {input.ref} {input.r1} {input.r2} > {output}
        """

rule sam_to_sorted_bam:
    input:
        "results/{sample}/aligned/{sample}.sam"
    output:
        "results/{sample}/aligned/{sample}_sorted.bam"
    shell:
        """
        samtools view -Sb {input} | samtools sort -o {output}
        """

rule mark_duplicates:
    input:
        bam="results/{sample}/aligned/{sample}_sorted.bam"
    output:
        bam="results/{sample}/aligned/{sample}_dedup.bam",
        metrics="results/{sample}/aligned/{sample}_dedup.metrics.txt"
    shell:
        """
        picard MarkDuplicates I={input.bam} O={output.bam} M={output.metrics} CREATE_INDEX=true
        """

rule base_recalibration:
    input:
        bam="results/{sample}/aligned/{sample}_dedup.bam",
        bai="results/{sample}/aligned/{sample}_dedup.bai",
        ref=REF,
        ref_dict=REF.rsplit(".", 1)[0] + ".dict",
        ref_fai=REF + ".fai"
    output:
        table="results/{sample}/recal/{sample}_recal_data.table"
    params:
        known_sites="references/known_sites.vcf"  # replace with your known sites vcf path
    shell:
        """
        gatk BaseRecalibrator -I {input.bam} -R {input.ref} --known-sites {params.known_sites} -O {output.table}
        """

rule apply_bqsr:
    input:
        bam="results/{sample}/aligned/{sample}_dedup.bam",
        bai="results/{sample}/aligned/{sample}_dedup.bai",
        recal_table="results/{sample}/recal/{sample}_recal_data.table",
        ref=REF
    output:
        bam="results/{sample}/recal/{sample}_recal.bam"
    shell:
        """
        gatk ApplyBQSR -R {input.ref} -I {input.bam} --bqsr-recal-file {input.recal_table} -O {output.bam}
        """

rule haplotype_caller:
    input:
        bam="results/{sample}/recal/{sample}_recal.bam",
        bai="results/{sample}/recal/{sample}_recal.bai",
        ref=REF
    output:
        vcf="results/{sample}/variants/{sample}.vcf.gz"
    shell:
        """
        gatk HaplotypeCaller -R {input.ref} -I {input.bam} -O {output.vcf} -ERC GVCF
        """
rule all:
    input:
        "output.txt"

rule create_output:
    output:
        "output.txt"
    shell:
        "echo 'Hello Snakemake!' > output.txt"

