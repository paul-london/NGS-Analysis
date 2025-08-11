# 🧬 NGS Analysis Pipeline
*Work in Progress*

## 🚀 Overview

This repository contains a Next-Generation Sequencing (NGS) data analysis pipeline. The environment leverages Python with bioinformatics libraries running on Linux to enable reproducible and efficient genomic data workflows. 

The current configuration is for variant calling following a sequencing run for a variety of genetic panels. 

The core workflow occurs with Snakemake while a Jupyter Notebook houses data analytics and visualizations.

---

## 🔄 Workflow

| Step Number | Step Name              | Description                                                                                   | Tools                                                                                       |
|:-----------:|------------------------|-----------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------|
| 1           | Data Preparation       | Upload or mount raw sequencing data (FASTQ files) into the project directory.                 | Linux command-line utilities (`rsync`, `scp`), VirtualBox shared folders, SRA-tools (NCBI)  |
| 2           | Quality Control        | Perform quality control analyses to assess sequencing data quality.                           | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [MultiQC](https://multiqc.info/)             |
| 3           | Alignment & Processing | Map reads to the reference genome and process alignments, including marking duplicates and base quality score recalibration (BQSR). Collect QC metrics before and after these steps. | BWA, Bowtie2, SAMtools, [GATK](https://gatk.broadinstitute.org/hc/en-us)                      |
| 4           | Variant Calling & Annotation | Detect genetic variants and annotate their functional effects.                                 | GATK HaplotypeCaller, SnpEff, ANNOVAR                                                       |
| 5           | Visualization & Reporting | Create plots, summaries, and reports to facilitate data interpretation.                        | Matplotlib, Seaborn, Jupyter notebooks                                                     |
| 6           | Iterative Analysis     | Refine analysis parameters, rerun workflows, and document findings for reproducibility.       | Jupyter notebooks, version control (Git)                                                   |


---

## 🖥️ Environment Setup

- **Linux VM:** CentOS 7 (or compatible) running the Jupyter Notebook server
- **Jupyter Notebook:** Accessible remotely on port `8888`
- **Python Environment:** Includes bioinformatics and data science packages like `Biopython`, `numpy`, `pandas`, and `matplotlib`
- **VS Code Integration:** Connect VS Code to the remote Jupyter server for seamless editing and execution
---

## 🔧 How to Run Jupyter Notebook

1. **Start the server on the Linux VM:**

   ```bash
   jupyter notebook --ip=0.0.0.0 --port=8888 --no-browser
