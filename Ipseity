import glob

REF = "/home/ipseity-gnt-i9/ACMG/GRCh37.p13.genome.fa"

# Define input files using glob
fastq_files = glob.glob("*.fastq.gz")

# Extract sample names from fastq files
sample_names = [file.split("_R")[0] for file in fastq_files]

# Define output files
FASTP_OUT_1 = "{sample}_R1_fp.fastq"
FASTP_OUT_2 = "{sample}_R2_fp.fastq"
SAM_OUT = "{sample}_aligned.sam"
BAM_OUT = "{sample}_aligned.bam"
BAM_GROUP = "{sample}_grouped.bam"
BAM_INDEX = "{sample}_grouped.bam.bai"
BAM_MARKED = "{sample}_marked.bam"
BAM_MARKED_INDEX = "{sample}_marked.bam.bai"
VCF_OUT = "{sample}_variants.vcf"
VEP_IN = "/home/ipseity-gnt-i9/ACMG/VEP_data/{sample}_variants.vcf"
VEP_OUT = "{sample}_vep.vcf"
BQSR_TABLE = "{sample}_recal.table"
POST_BQSR_TABLE = "{sample}_post_recal.table"
BAM_RECAL = "{sample}_recal.bam"
RECAL_PLOT = "{sample}_recal_plots.pdf"
VCF_NO_STRAND_BIAS = "{sample}_no_strand_bias.vcf"
VCF_FILTERED = "{sample}_filtered.vcf"

rule all:
    input:
        expand(VEP_OUT, sample=sample_names)

rule fastp:
    input:
        fp1=lambda wildcards: f"{wildcards.sample}_R1_001.fastq.gz",
        fp2=lambda wildcards: f"{wildcards.sample}_R2_001.fastq.gz"
    output:
        fo1=FASTP_OUT_1,
        fo2=FASTP_OUT_2
    log:
        "fastp_{sample}.log"
    shell:
        "fastp -i {input.fp1} -I {input.fp2} -o {output.fo1} -O {output.fo2}"

rule bowtie2:
    input:
        bp1=FASTP_OUT_1
        bp2=FASTP_OUT_2
    output:
        SAM_OUT
    threads: 16
    log:
        "bowtie2_{sample}.log"
    shell:
        "bowtie2 -x {config[bowtie_index]} -1 {input.bp1} -2 {input.bp2} -S {output} -p {threads}"
    
rule samtools_view_sort:
    input:
        SAM_OUT
    output:
        BAM_OUT
    threads: 16
    log:
        "samtools_view_sort_{sample}.log"
    shell:
        "samtools view -bS {input} --threads {threads} | "
        "samtools sort -o {output} --threads {threads}"

rule picard_AddOrReplaceReadGroups:
    input:
        BAM_OUT
    output:
        BAM_GROUP
    log:
        "picard_{sample}.log"
    shell:
        "picard AddOrReplaceReadGroups I={input} O={output} RGID={wildcards.sample} RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20"

rule samtools_index:
    input:
        BAM_GROUP
    output:
        BAM_INDEX
    threads: 16
    log:
        "samtools_index_{sample}.log"
    shell:
        "samtools index {input} -@ {threads}"

rule picard_MarkDuplicates:
    input:
        BAM_GROUP
    output:
        BAM_MARKED,
        BAM_MARKED_INDEX
    log:
        "picard_mark_duplicates_{sample}.log"
    shell:
        "picard MarkDuplicates I={input} O={output[0]} M=marked_dup_metrics.txt REMOVE_DUPLICATES=true && "
        "samtools index {output[0]}"
    
rule gatk_BaseRecalibrator:
    input:
        ref=REF,
        bam=BAM_MARKED,
        bai=BAM_MARKED_INDEX
    output:
        BQSR_TABLE
    threads: 16
    log:
        "gatk_BaseRecalibrator_{sample}.log"
    shell:
        "gatk BaseRecalibrator -R {config[ref]} -I {input.bam} --known-sites {config[known_sites]} -O {output}"

rule gatk_ApplyBQSR:
    input:
        ref=REF,
        bam=BAM_MARKED,
        bai=BAM_MARKED_INDEX,
        bqsr_table=BQSR_TABLE
    output:
        recal_bam=BAM_RECAL,
        post_bqsr_table=POST_BQSR_TABLE,
        recal_plot=RECAL_PLOT
    threads: 16
    log:
        "gatk_ApplyBQSR_{sample}.log"
    shell:
        "gatk ApplyBQSR -R {config[ref]} -I {input.bam} -bqsr {input.bqsr_table} -O {output.recal_bam} && gatk BaseRecalibrator -R {config[ref]} -I _______"

rule gatk_HaplotypeCaller:
    input:
        ref=REF,
        bam=BAM_RECAL,
        bai=BAM_MARKED_INDEX
    output:
        VCF_OUT
    threads: 16
    log:
        "gatk_HaplotypeCaller_{sample}.log"
    shell:
        "gatk HaplotypeCaller -R {config[ref]} -I {input.bam} -O {output} --native-pair-hmm-threads {threads} -L {config[bed_file]}"

rule gatk_FilterStrandBias:
    input:
        VCF_OUT
    output:
        VCF_NO_STRAND_BIAS
    log:
        "gatk_FilterStrandBias_{sample}.log"
    shell:
        "bcftools filter -i \'INFO/SOR < 2.0 || INFO/SOR > 0.5' {input} -o {output}"
#        "gatk VariantFiltration -R {config[ref]} -V {input} --filter-expression \"SOR>2\" --filter-name \"StrandBias\" -O {output}"

rule bcftools_filter:
    input:
        VCF_NO_STRAND_BIAS
    output:
        VCF_FILTERED
    log:
        "bcftools_filter_{sample}.log"
    params:
        DP=config["filter_thresholds"]["DP"]
        GQ=config["filter_thresholds"]["GQ"]
    shell:
        "bcftools filter -i '(FORMAT/DP > {params.DP} & FORMAT/GQ >= {params.GQ}) & (GT=\"0/1\" || GT=\"1/0\" || GT=\"2/0\" || GT=\"1/2\")"

rule copy_variants:
    input:
        VCF_FILTERED
    output:
        VEP_IN
    log:
        "copy_variants_{sample}.log"
    shell:
        "cp {input} {output}"

rule ensembl_vep:
    input:
        vcf=VCF_OUT
        vep=VEP_IN
    output:
        touch(VEP_OUT)
    log:
        "ensembl_vep_{sample}.log"
    shell:
        """
        docker run -v {config[docker_data]}:/data:Z -it ensemblorg/ensembl-vep vep --form 50 --cache --sift b --polyphen b --ccds --symbol --canonical
        """