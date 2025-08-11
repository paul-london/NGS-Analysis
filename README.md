# 🧬 NGS Analysis Pipeline
*Work in Progress*

## 🚀 Overview

This repository contains a Next-Generation Sequencing (NGS) data analysis pipeline. The environment leverages Python with bioinformatics libraries running on Linux to enable reproducible and efficient genomic data workflows. The current configuration is for variant calling following a sequencing run for a variety of genetic panels. A Jupyter Notebook houses data analytics and visualizations while the core workflow occurs with Snakemake.

---

## 🔄 Workflow

1. **Data Preparation**  
   Upload or mount raw sequencing data (FASTQ files) into the project directory.  
   _Tools: Linux command-line utilities (`rsync`, `scp`), VirtualBox shared folders, SRA-tools (NCBI)_

2. **Quality Control**  
   Run QC analyses to check sequencing quality.  
   _Tools: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [MultiQC](https://multiqc.info/)_

3. **Alignment & Processing**  
   Map reads to a reference genome and process alignments.  
   _Tools: BWA, Bowtie2, SAMtools, [GATK](https://gatk.broadinstitute.org/hc/en-us)_

4. **Variant Calling & Annotation**  
   Identify genetic variants and annotate their effects.  
   _Tools: GATK HaplotypeCaller, SnpEff, ANNOVAR_

5. **Visualization & Reporting**  
   Generate plots, summaries, and reports for data interpretation.  
   _Tools: Matplotlib, Seaborn, Jupyter notebooks_

6. **Iterative Analysis**  
   Refine parameters, rerun analyses, and document insights.  
   _Tools: Jupyter notebooks, version control (Git)_

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
