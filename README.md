# NGS Mutation Analysis Pipeline

## Overview

This repository contains a series of scripts designed to analyze next-generation sequencing (NGS) data to quantify mutation frequencies, compare mutation frequencies across different NGS datasets, and calculate the enrichment of specific mutations under various experimental conditions. The pipeline is designed to handle data that includes unique molecular identifiers (UMIs) to generate consensus sequences, thereby filtering out PCR and sequencing artifacts, and ensuring accurate mutation detection.

### Key Features:
- **Consensus Sequence Generation**: Group sequences by UMIs and derive consensus sequences to minimize PCR artifacts.
- **Variant Calling**: Use the ASMV software for variant calling on consensus sequences.
- **Mutation Frequency Calculation**: Quantify the frequency of mutations within different NGS datasets.
- **Enrichment Analysis**: Compare mutation frequencies across datasets to calculate the enrichment of mutations.
- **Visualisation**: Generate visualisations, including frequency plots, enrichment graphs, and pie charts, to analyse mutation distributions.

## Scripts

### 1. `consensus_script.py`

**Purpose**: This script processes NGS data to group sequences by UMIs, generate consensus sequences, and prepare data for variant calling.

**Key Functions**:
- `extract_umi_and_sequence(fastq_file)`: Extracts UMIs and corresponding sequences from a FASTQ file.
- `group_by_umi(records, output_file, discarded_output_file, min_family_size=2)`: Groups sequences by UMI and outputs grouped sequences.
- `find_consensus_sequence(sequences, qualities, threshold=1.0, am_threshold=0.4)`: Derives a consensus sequence for each UMI group.
- `process_umi_groups(grouped_records, consensus_output_file, min_family_size=2, consensus_threshold=1.0, am_threshold=0.4)`: Processes UMI groups to find consensus sequences.

**Usage**:
```bash
python consensus_script.py

### 2. `process_output_script.sh`

**Purpose**: This Bash script processes the output from the consensus script by aligning reads, converting file formats, and running the ASMV variant calling tool.

**Key Steps**:

- Activates the required conda environment.
- Aligns reads using BWA and converts SAM to BAM files using Samtools.
- Runs GATK AnalyzeSaturationMutagenesis to identify mutations.

**Usage**:
```bash
./process_output_script.sh

### 3. `enrichment_script_nucleotide.py`

**Purpose**: This script calculates mutation frequency ratios and performs enrichment analysis based on nucleotide sequences.

**Key Functions**:

- `read_and_process_data(file_path)`: Reads and processes variant count files.
- `calculate_ratios(data, max_mutations)`: Calculates mutation frequency ratios.
- `calculate_enrichment()`: Calculates enrichment values by comparing mutation frequencies across datasets.

**Usage**:
```bash
python enrichment_script_nucleotide.py

### 4. `enrichment_script_amino_acid.py`

**Purpose**: This script performs enrichment analysis based on amino acid sequences, similar to the nucleotide script.

**Key Functions**:

- `read_and_process_data(file_path)`: Reads and processes variant count files.
- `calculate_ratios(data, max_nucleotide_changes)`: Calculates mutation frequency ratios based on amino acid changes.

**Usage**:
```bash
python enrichment_script_amino_acid.py

### 5. Cleaning Scripts

**Purpose**: These scripts clean the frequency and enrichment matrices by applying specific filters based on amplicon ranges.

**Scripts**:

- `clean_frequency_matrix_nucleotide.py`
- `clean_enrichment_matrix_nucleotide.py`

**Usage**:
python clean_frequency_matrix_nucleotide.py
python clean_enrichment_matrix_nucleotide.py

### 6. Jupyter Notebook: `visualization_notebook.ipynb`

**Purpose**: A Jupyter Notebook for generating various visualizations such as frequency plots, enrichment graphs, pie charts, and violin plots.

**Usage**:

- Open the notebook using Jupyter and run the cells to generate visualizations.

## Installation

1. Clone this repository:

    ```bash
    git clone https://github.com/yourusername/NGS_Mutation_Analysis_Pipeline.git
    ```

2. Install the required dependencies using conda:

    ```bash
    conda create -n umi_env -f environment.yml
    conda activate umi_env
    ```

## Usage

Follow the steps below to execute the pipeline:

1. **Generate Consensus Sequences**: Run `consensus_script.py` to process your FASTQ files and generate consensus sequences.
2. **Process Outputs**: Use `process_output_script.sh` to align sequences and perform variant calling.
3. **Calculate Enrichment**: Run the enrichment scripts to calculate mutation enrichment based on nucleotides or amino acids.
4. **Clean Data**: Use the cleaning scripts to clean your frequency and enrichment matrices.
5. **Visualize Results**: Use the Jupyter notebook to visualize the results.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For any inquiries or contributions, please contact `yourname`.

