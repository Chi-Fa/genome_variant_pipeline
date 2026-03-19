# Yeast Variant Calling Pipeline (GATK Best Practices)

This repository contains an automated Bash-based pipeline for processing **Yeast (S288C)** paired-end sequencing data from Aviti-seq. The workflow strictly follows the **GATK Best Practices** for germline short variant discovery, including alignment, preprocessing, BQSR (Base Quality Score Recalibration), variant calling, and functional annotation.

## Pipeline Overview

The script (`[Your_Script_Name].sh`) automatically executes the following stages:

1.  **Reference Indexing**: Generates BWA index for the reference genome.
2.  **Sequence Alignment**: Maps paired-end FASTQ reads using `bwa mem`.
3.  **Mark Duplicates**: Identifies and tags duplicate reads via `gatk MarkDuplicatesSpark`.
4.  **Quality Metrics**: Collects alignment and insert size metrics using `Picard Tools`.
5.  **Coverage Analysis**: Calculates per-base depth using `samtools depth`.
6.  **Initial Variant Calling**: Performs `HaplotypeCaller` (Default ploidy=50 for mixed populations).
7.  **Hard Filtering**: Filters SNPs and INDELs based on GATK standard thresholds.
8.  **BQSR (Base Quality Score Recalibration)**: Uses high-confidence variants from the first round to recalibrate base quality scores.
9.  **Final Variant Calling**: Runs `HaplotypeCaller` on recalibrated BAM files for higher accuracy.
10. **Final Filtering**: Final pass of hard filtering for SNPs/INDELs.
11. **Annotation**: Predicts functional effects of variants using `snpEff`.

## Prerequisites

Ensure the following tools are installed and available in your `PATH`:

* [BWA](https://github.com/lh3/bwa)
* [GATK4 (v4.x+)](https://github.com/broadinstitute/gatk)
* [Picard Tools](https://broadinstitute.github.io/picard/)
* [Samtools](https://github.com/samtools/samtools)
* [snpEff](http://pcingola.github.io/SnpEff/)
    * *Note: Ensure the corresponding database (e.g., `GCF_000146045.2`) is configured.*
* Java JRE/JDK (Required for GATK and Picard)

## Directory Structure

Set up your project directory as follows before running the script:

```text
project_root/
│
├── [Your_Script_Name].sh      # The analysis script
├── README.md                  # Project documentation
├── .gitignore                 # Prevents uploading large data files
│
├── data/
│   ├── reference/
│   │   └── genomic.fna    # Reference genome
│   │
│   └── raw_reads/                             # Raw FASTQ files
│       ├── SampleA_R1.trim.fq.gz
│       ├── SampleA_R2.trim.fq.gz
│       └── SampleB_R1.trim.fq.gz
│       └── SampleB_R2.trim.fq.gz
│
└── results/                                    # Generated automatically


## Post-processing (Matrix Generation)

After variant annotation, you can use the provided Python script to merge individual sample results into a unified matrix:

### Requirements
* Python 3.x
* Pandas (`pip install pandas`)

### Execution
The script `merge_variants.py` will look for `*.out` files in the `results/` directory and generate consolidated CSV matrices for Allele Frequency (AF), Gene names, Codon changes, and Amino Acid variations.

```bash
python3 merge_variants.py


## Evolutionary Analysis (Parallel Evolution)

The script `analyze_parallel_evolution.py` combines mutation data with annotation matrices to identify genes under strong selection pressure.

### Features
* Integrates Codon and Amino Acid change information.
* Calculates a **Selection Score** based on founder groups, replicate hits, and total SNP count.
* Generates ranking tables for both **Gene** and **AA-residue** levels.
* Produces visualization plots for top candidates.

### Usage
```bash
python3 analyze_parallel_evolution.py -i [INPUT_DIR] -m [MATRIX_DIR] -o [OUTPUT_DIR]
