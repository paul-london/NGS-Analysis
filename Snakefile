import glob
import re

REF = "NGS-Analysis-Raw-Files/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# Pattern: <sampleID>_<panel>_<seqIndex>_R1_001.fastq.gz
R1_files = glob.glob("NGS-Analysis-Raw-Files/reads/*_R1_001.fastq.gz")

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
        expand("results/vcf/{sample}.vcf.gz", sample=SAMPLES.keys())

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

rule align_reads:
    input:
        ref=REF,
        ref_fai=REF + ".fai",
        ref_dict=REF.rsplit(".", 1)[0] + ".dict",
        R1=lambda w: SAMPLES[w.sample]["R1"],
        R2=lambda w: SAMPLES[w.sample]["R2"]
    output:
        temp("output/{sample}.bam")
    shell:
        """
        bwa mem -t 8 {input.ref} {input.R1} {input.R2} | \
        samtools view -bS - | \
        samtools sort -o {output}
        samtools index {output}
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
