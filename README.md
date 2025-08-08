# 🧬 NGS Analysis Project

## 🚀 Overview

This repository contains Jupyter notebooks for **Next-Generation Sequencing (NGS)** data analysis running on a Linux Virtual Machine (VM). The environment leverages Python with bioinformatics libraries to enable reproducible and efficient genomic data workflows.

---

## 🖥️ Environment Setup

- **Linux VM:** CentOS 7 (or compatible) running the Jupyter Notebook server
- **Jupyter Notebook:** Accessible remotely on port `8888`
- **Python Environment:** Includes bioinformatics and data science packages like `Biopython`, `numpy`, `pandas`, and `matplotlib`
- **VS Code Integration:** Connect VS Code to the remote Jupyter server for seamless editing and execution
- **Shared Folders:** Project files available via VirtualBox shared folder mounted at `/media/sf_NGS-Analysis` with a symlink for convenience

---

## 🔧 How to Run Jupyter Notebook

1. **Start the server on the Linux VM:**

   ```bash
   jupyter notebook --ip=0.0.0.0 --port=8888 --no-browser
