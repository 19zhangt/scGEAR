## **SERA** - A Computational Framework for Analyzing Multiple RNA-Based Phenotypes in scRNA-Seq Data

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)  [![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B-blue)](https://www.python.org/)  [![R 4.2.2](https://img.shields.io/badge/R-4.2.2-blue)](https://www.r-project.org/)

<p align="center"><img  src="assets/SERA.png" width="60%" /></p>

Currently, **SERA** quantitatively analyzes three RNA-based phenotypes at single-cell resolution: gene expression, alternative polyadenylation (APA), and RNA editing, using polyA-enriched single-cell RNA-sequencing data. With its efficient design, SERA enables researchers to uncover insights into RNA processing and cellular heterogeneity, providing a comprehensive view of transcriptional diversity at the single-cell level.

**Key Features**:

- **Single-cell resolution**: Works with polyA-enriched single-cell RNA-seq data.
- **Multimodal post-transcriptional modification**: Identifies multiple RNA modification types (e.g., APA and RNA editing).
- **Population-scale analysis**: Scalable to large datasets for population-level studies.

---
## Table of Contents
- [Installation](#installation-)
- [Quick Start](#quick-start-)
- [Data Requirements](#data-requirements-)
- [Output](#output)


## Installation 🔧

1. Clone this repository:
   ```bash
   git clone https://github.com/lilab-bioinfo/SERA.git
   cd SERA
   ```

2. Install dependencies with Mamba
   ```bash
   conda install -n base -c conda-forge mamba
   mamba env create -f environment.yml
   conda activate SERA
   Rscript install_dependencies.R  ## If dependencies are not installed, please run this code.
   ```

3. Download genome
   ```bash
   wget -P reference/ ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
   gunzip reference/GRCh38.primary_assembly.genome.fa.gz
   samtools faidx reference/GRCh38.primary_assembly.genome.fa
   ```

## Quick Start 🚀

**SERA** provides a command-line interface (CLI) for ease of use. Below is a basic example of how to analyze single-cell RNA-seq data.

### For gene expression
```bash
Rscript SERA.R --mode expression \
  --cellinfo example/cell_annotation.csv \
  --baminfo example/bam_list.csv \
  --threads 10 \
  --outdir example/result
```

### For poly(A) sites and alternative polyadenylation events
```bash
Rscript SERA.R --mode APA \
  --cellinfo example/cell_annotation.csv \
  --baminfo example/bam_list.csv \
  --genome_fa reference/GRCh38.primary_assembly.genome.fa \
  --threads 10 \
  --outdir example/result
```

### For RNA editing sites
```bash
Rscript SERA.R --mode editing \
  --cellinfo example/cell_annotation.csv \
  --baminfo example/bam_list.csv \
  --genome_fa reference/GRCh38.primary_assembly.genome.fa \
  --remove_duplicates TRUE \
  --threads 10 \
  --outdir example/result
```

### For all phenotypes
```bash
Rscript SERA.R --mode all \
  --cellinfo example/cell_annotation.csv \
  --baminfo example/bam_list.csv \
  --genome_fa reference/GRCh38.primary_assembly.genome.fa \
  --remove_duplicates TRUE \
  --threads 10 \
  --outdir example/result
```

## Data Requirements 📂

SERA requires outputs from **Cell Ranger** (10x Genomics). Users must specify the following parameters:

- `--mode`: Specify the analysis mode (`expression`, `apa`, `editing`, `all`).
- `--cellinfo`: Path to a CSV file containing cell metadata with four columns: `pool`, `barcode`, `cell_type`, and `donor`. Example:
  ```text
  pool,barcode,cell_type,donor
  S1,AAACCCACACCATTCC-1,CD4,ind1
  S1,AAACCCACATCCGATA-1,CD8,ind2
  S2,AAACCCATCTCTGCTG-1,NK,ind3
  S2,AAACGAACATTGAGGG-1,Mono,ind4
  ```

- `--baminfo`: the output directory of cellranger bam files.
  ```text
  pool,bamfile
  S1,data/cellranger/S1/outs/possorted_genome_bam.bam
  S2,data/cellranger/S2/outs/possorted_genome_bam.bam
  ```

- `--genome_fa`: genome sequence file in fasta format.
- `--outdir`: Path to the output directory where results will be saved.


## Output

The output directory will contain **matrices** of multimodal molecular phenotypes for each cell type, with rows representing phenotypes and columns representing sample IDs.

- ​​Normalized Gene Expression Matrices​​
  - ​Directory​​: `example/result/expression/cell_type_quant/` (​Normalized expression matrices across cell types.)
- ​​Poly(A) Site and APA Event Analysis​​
  - ​Poly(A) Sites​​: `example/result/APA/Annotation/Selected_PAS_sites.txt`
  - ​​Independent APA Events​​: `example/result/APA/Annotation/Independent_APA_events.txt`
  - ​APA Usage Matrices​​: `example/result/APA/cell_type_quant/`
- ​​RNA Editing Analysis Results​​
  - ​Editing Level Matrices​​: `example/result/editing/cell_type_quan/`


## Getting help 📬

To report a bug or request a feature, please [open an Issue](https://github.com/lilab-bioinfo/SERA/issues).

For general discussions about single-cell analysis, please [start a Discussion](https://github.com/lilab-bioinfo/SERA/discussions).

## Release Notes​​ 🎉

v1.0.0 (2025-04-21): Initial release with support for gene expression, APA, and RNA editing analysis.

---
