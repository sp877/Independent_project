# NGS Mutation Analysis Pipeline

This repository contains a suite of scripts designed to process and analyse Next-Generation Sequencing (NGS) data. The primary goal is to leverage unique molecular identifiers (UMIs) to reduce sequencing errors and accurately identify true mutations and to use this to quantify mutation frequencies, calculate enrichment of mutations, and visualise these mutations under various experimental conditions.

## Table of Contents

1. [Overview](#overview)
2. [Pipeline Components](#pipeline-components)
   - [1. Consensus Sequence Generation](#1-consensus-sequence-generation)
   - [2. Variant Calling and Processing](#2-variant-calling-and-processing)
   - [3. Enrichment Calculation](#3-enrichment-calculation)
   - [4. Data Cleaning](#4-data-cleaning)
   - [5. Visualisation](#5-visualisation)
3. [Installation](#installation)
4. [Usage](#usage)
   - [Running the Pipeline](#running-the-pipeline)

## Overview

This software package is developed to quantify the frequencies of each mutation in NGS sequences and compare these frequencies across different datasets to calculate mutation enrichment. The pipeline involves the following steps:

1. **Grouping sequences by UMI:** Generates consensus sequences to minimise errors.
2. **Variant Calling:** Identifies mutations from consensus sequences.
3. **Mutation Frequency and Enrichment Analysis:** Compares mutation frequencies across different experimental conditions.
4. **Data Cleaning:** Refines the results by removing artifacts and irrelevant data.
5. **Visualization:** Generates plots to visualise mutation frequencies and enrichments.

## Pipeline Components

### 1. Consensus Sequence Generation

The `consensus_script.py` script processes NGS data by grouping sequences based on UMIs. This step ensures that only real mutations are retained by generating consensus sequences from UMI families.

- **Inputs:** FASTQ files with UMIs.
- **Outputs:** FASTQ files containing consensus sequences.

### 2. Variant Calling and Processing

The `processing_script.sh` script aligns the reads, performs variant calling using the GATK `AnalyzeSaturationMutagenesis` tool, and prepares the data for further analysis.

- **Inputs:** Consensus sequence FASTQ files.
- **Outputs:** BAM files, variant counts.

### 3. Enrichment Calculation

The enrichment scripts (`enrichment__script_nt.py` and `enrichment_script_aa.py`) calculate mutation frequencies and enrichment values by comparing different experimental groups.

- **Inputs:** Variant count files.
- **Outputs:** Enrichment matrices for nucleotides and amino acids.

### 4. Data Cleaning

The cleaning scripts (`clean_frequency_matrix_nt.py`, `clean_enrichment_matrix_nt.py`, `clean_frequency_matrix_aa.py` and `clean_enrichment_matrix_aa.py`) remove irrelevant data points and ensure the matrices contain only relevant data within specified amplicon ranges.

- **Inputs:** Frequency and enrichment matrices.
- **Outputs:** Cleaned matrices ready for visualization.

### 5. Visualisation

Jupyter notebooks (`visualisations_nucleotide.ipynb` and `visualisations_aminoacids.ipynb` ) are used to generate various plots, including frequency graphs, enrichment maps, and violin plots, which help in understanding the mutation landscape.

- **Inputs:** Cleaned matrices.
- **Outputs:** Graphs and charts (e.g., frequency plots, enrichment maps).

## Installation
   
### Install Python dependencies:

Make sure you have Python 3 installed. Then, install the required Python packages using pip:

```bash
pip install -r requirements.txt
``` 

### Set up the Conda environment: 

If you haven't already, install Conda (Miniconda or Anaconda). Then, create and activate the environment needed for bioinformatics tools:

```bash
conda create -n umi_env
conda activate umi_env
conda install -c bioconda bwa samtools gatk4 picard
```
This will install the following tools:

- **BWA**: For aligning reads.
- **Samtools**: For processing BAM files.
- **GATK**: For variant calling.
- **Picard**: For manipulating sequence files.

### Download and prepare reference data:

Ensure you have the necessary reference sequences and index them as needed for BWA and GATK.

## Usage 

### Running the Pipeline

Follow these steps to run the complete NGS mutation analysis pipeline:

1. **Generate Consensus Sequences:**

   Run the script to generate consensus sequences from your FASTQ files:

   ```bash
   python consensus_script.py
   ```
   
2. **Process the Output:**

   Execute the shell script to align reads, perform variant calling, and process the BAM files:

   ```bash
   bash Processing_script.sh
   ```

3. **Calculate Enrichment:**

   Run the enrichment calculation scripts to compare mutation frequencies across different experimental groups:

   For Nulceotides:
   ```bash
   python enrichment_script_nt.py
   ```
   
   For amino acids:
   ```bash
   python enrichment_script_aa.py
   ```
   
4. **Clean the Data:**

   Clean the frequency and enrichment matrices to remove irrelevant data points:

   For Nulceotides:
   ```bash
   python clean_frequency_matrix_nt.py
   ```
   ```bash
   python clean_enrichment_matrix_nt.py
   ```
   
   For amino acids:
   ```bash
   python clean_frequency_matrix_aa.py
   ```
   ```bash
   python clean_enrichment_matrix_aa.py
   ```

5. **Visualisation:**

   Use the provided Jupyter notebook (`visualisations_nucleotide.ipynb` and `visualisations_aminoacids.ipynb` ) to generate various plots. Launch Jupyter Notebook in your terminal:

   ```bash
   Jupyter Notebook
   ```
