# 🧬 NGS Analysis Pipeline
*Work in Progress*

## 🚀 Overview

This repository contains a Next-Generation Sequencing (NGS) data analysis pipeline. The environment leverages Python with bioinformatics libraries running on Linux to enable reproducible and efficient genomic data workflows. 

The current configuration is for **variant calling** following a sequencing run for a variety of genetic panels assessing pathogenicity of variants, if present. 

The core workflow occurs with Snakemake while a Jupyter Notebook houses data analytics and visualizations.

---

## 🧬 Workflow

| Step Number | Step Name              | Description                                                                                   | Tools                                                                                       |
|:-----------:|------------------------|-----------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------|
| 1           | Data Preparation       | Upload or mount raw sequencing data (FASTQ files) into the project directory.                 | Linux command-line utilities (`rsync`, `scp`), Linux virtual machine, [SRA-tools (NCBI)](https://github.com/ncbi/sra-tools)  |
| 2           | Quality Control        | Perform quality control analyses to assess sequencing data quality.                           | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [MultiQC](https://multiqc.info/)             |
| 3           | Alignment & Processing | Map reads to the reference genome and process alignments, including marking duplicates and base quality score recalibration (BQSR). Collect QC metrics before and after these steps. | [BWA](http://bio-bwa.sourceforge.net/), [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [SAMtools](http://www.htslib.org/), [GATK](https://gatk.broadinstitute.org/hc/en-us)                      |
| 4           | Variant Calling & Annotation | Detect genetic variants and annotate their functional effects.                                 | [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller), [SnpEff](http://snpeff.sourceforge.net/), [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/)                                                       |
| 5           | Visualization & Reporting | Create plots, summaries, and reports to facilitate data interpretation.                        | [Matplotlib](https://matplotlib.org/), [Seaborn](https://seaborn.pydata.org/), [Jupyter notebooks](https://jupyter.org/)                                                     |
| 6           | Iterative Analysis     | Refine analysis parameters, rerun workflows, and document findings for reproducibility.       | [Jupyter notebooks](https://jupyter.org/), [Git](https://git-scm.com/)                                                   |

---

## 📁 Required Files

Before running the pipeline, please ensure the following files and directories are available and correctly placed:

| File/Directory                              | Description                                       | Example Location                        |
|--------------------------------------------|-------------------------------------------------|----------------------------------------|
| Raw sequencing data (FASTQ files)           | Input reads, paired-end files                    | `data/raw_reads/{sample}_R1_001.fastq.gz` and `_R2_001.fastq.gz` |
| Reference genome FASTA                      | Reference genome for alignment                    | `reference/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa`    |
| Reference genome index files                | Index files for the reference (will be created if missing) | Same directory as reference FASTA      |
| BED files for target panels                  | Target regions for variant calling                | `reference/panel/{panel}.bed`           |
| Known sites VCF                             | Known variant sites for base recalibration       | `reference/known_sites.vcf`             |

If any file is missing, please obtain or generate it before running the workflow. Make sure all paths match those expected by the pipeline configuration.

---

## ⚙️ Environment Setup

- **Linux VM:** CentOS 7 (or compatible) running the Jupyter Notebook server
- **Jupyter Notebook:** Accessible remotely on port `8888` after starting the server:

     ```bash
   jupyter notebook --ip=0.0.0.0 --port=8888 --no-browser
     
- **Python Environment:** Includes bioinformatics and data science packages like `Biopython`, `numpy`, `pandas`, and `matplotlib`
- **VS Code Integration:** Connect VS Code to the remote Jupyter server for seamless editing and execution

To ensure all necessary tools and dependencies are installed, create and activate the conda environment (pipeline) using the provided `environment.yaml` file:

```bash
conda env create -f environment.yaml
conda activate pipeline

