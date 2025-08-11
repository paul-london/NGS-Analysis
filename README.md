# 🧬 NGS Analysis Pipeline
*Work in Progress*

## 🚀 Overview

This repository contains a Next-Generation Sequencing (NGS) data analysis pipeline. The environment leverages Python with bioinformatics libraries running on Linux to enable reproducible and efficient genomic data workflows. 

The current configuration is for **variant calling** following a sequencing run for a variety of genetic panels assessing pathogenicity of any detected variants. 

The core workflow occurs with Snakemake while a Jupyter Notebook houses data analytics and visualizations.

---

## 🛠️ Workflow

| Step Number | Step Name              | Description                                                                                   | Tools                                                                                       |
|:-----------:|------------------------|-----------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------|
| 1           | Data Preparation       | Upload, download, or mount raw sequencing data (FASTQ files) into the project directory.                 | Linux command-line utilities (`rsync`, `scp`), Linux virtual machine, [SRA-tools (NCBI)](https://github.com/ncbi/sra-tools)  |
| 2           | Quality Control        | Perform quality control analyses to assess sequencing data quality.                           | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [MultiQC](https://multiqc.info/)             |
| 3           | Alignment & Processing | Map reads to the reference genome and process alignments, including marking duplicates and base quality score recalibration (BQSR). Collect QC metrics before and after these steps. | [BWA](http://bio-bwa.sourceforge.net/), [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [SAMtools](http://www.htslib.org/), [GATK](https://gatk.broadinstitute.org/hc/en-us)                      |
| 4           | Variant Calling & Annotation | Detect genetic variants and annotate their functional effects. Flag samples with pathogenic variants as Positive. Results will be compiled into a .tsv.               | [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller), [SnpEff](http://snpeff.sourceforge.net/), [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/)                                                       |
| 5           | Visualization & Reporting | Create plots, summaries, and reports to facilitate data interpretation.                        | [Matplotlib](https://matplotlib.org/), [Seaborn](https://seaborn.pydata.org/), [Jupyter notebooks](https://jupyter.org/)                                                     |
| 6           | Iterative Analysis     | Refine analysis parameters, streamline workflow, and document findings for reproducibility.       | [Jupyter notebooks](https://jupyter.org/), [Git](https://git-scm.com/)                                                   |

---

## 📁 Required Files

Before running the pipeline, please ensure the following files and directories are available and correctly placed:

| Data                                       | File Type              |Description                                       | Location and Format                        |
|--------------------------------------------|----------------|-------------------------------------------------|----------------------------------------|
| Raw sequencing data         | FASTQ          | Input reads, paired-end files                    | `data/raw_reads/{sample}_R1_001.fastq.gz` and `{sample}_R2_001.fastq.gz` |
| Reference genome                    | FASTA          | Reference genome for alignment                    | `reference/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa`    |
| Reference genome index files               | Multiple       | Index files for the reference (will be created if missing) | Same directory as reference FASTA      |
| Genes included in panel                | BED            | Target regions (genes) for variant calling                | `reference/panel/{panel}.bed`           |
| Known variants                           | VCF            | Known variant sites for base recalibration       | `reference/variants/All_20180418.vcf.gz`             |

If any file is missing, please obtain or generate it before running the workflow. Make sure all paths match those expected by the pipeline configuration.

---

## ⚙️ Running the Pipeline

### Environment Setup

- Linux Environment: CentOS 7 (or compatible) running Jupyter Notebook server
- Jupyter Notebook: Accessible remotely on port `8888` after starting the server:

     ```bash
   jupyter notebook --ip=0.0.0.0 --port=8888 --no-browser
     
- Python Environment: Includes bioinformatics and data science packages like `Biopython`, `numpy`, `pandas`, and `matplotlib`
- VS Code Integration: Connect VS Code to the remote Jupyter server for data analysis

To ensure all necessary tools and dependencies are installed, create and activate the conda environment (pipeline) using the provided `environment.yaml` file:

```bash
# Create environment from environment.yaml (default name: pipeline)
conda env create -f environment.yaml
# Activate environment
conda activate pipeline
```

Then run the pipeline:
```bash
# Simulate workflow
snakemake -n

# Activate workflow
snakemake -cores #
```


