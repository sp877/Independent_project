#!/bin/bash

# Function to check if a command exists
command_exists() {
    command -v "$1" &> /dev/null
}

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate umi_env

# Check if conda activation worked
if [ $? -ne 0 ]; then
    echo "Error: Conda environment 'umi_env' could not be activated"
    exit 1
fi

# Check if BWA is installed
if ! command_exists bwa; then
    echo "Error: BWA is not installed or not found in the PATH."
    exit 1
fi

# Check if Samtools is installed
if ! command_exists samtools; then
    echo "Error: Samtools is not installed or not found in the PATH."
    exit 1
fi

# Define Picard JAR file path
PICARD_JAR_PATH="/home/s/sp877/miniconda3/envs/bioenv/share/picard-3.1.1-0/picard.jar"

# Check if Picard JAR file exists
if [ ! -f "$PICARD_JAR_PATH" ]; then
    echo "Error: Picard JAR file not found at $PICARD_JAR_PATH"
    exit 1
fi

# Set working directory
cd "/scratch/alice/s/sp877/IP/without_umi" || { echo "Error: Directory not found"; exit 1; }

# Define input file
input_file="010014_merged.fq"
base="14test"

# Check if dictionary file exists, if not, create it
if [ ! -f SgeneRBD.dict ]; then
    echo "Creating dictionary for SgeneRBD.fasta"
    java -jar "$PICARD_JAR_PATH" CreateSequenceDictionary R=SgeneRBD.fasta O=SgeneRBD.dict
    if [ $? -ne 0 ]; then
        echo "Error: Failed to create sequence dictionary"
        exit 1
    fi
fi

# Create BWA index and SAMtools FASTA index
if [ ! -f SgeneRBD.fasta.bwt ]; then
    bwa index -a is SgeneRBD.fasta
    if [ $? -ne 0 ]; then
        echo "Error: Failed to create BWA index"
        exit 1
    fi
fi

if [ ! -f SgeneRBD.fasta.fai ]; then
    samtools faidx SgeneRBD.fasta
    if [ $? -ne 0 ]; then
        echo "Error: Failed to create SAMtools FASTA index"
        exit 1
    fi
fi

# Align reads
# Print message before alignment
echo "Aligning reads for ${base}" 
bwa mem SgeneRBD.fasta "${input_file}" > "${base}.sam"
if [ $? -ne 0 ]; then
    echo "Error: BWA alignment failed"
    exit 1
fi

# Make BAM files
# Print message before BAM conversion
echo "Converting SAM to BAM for ${base}"  
samtools view -S -b "${base}.sam" -o "${base}.bam"
if [ $? -ne 0 ]; then
    echo "Error: SAM to BAM conversion failed"
    exit 1
fi

samtools sort "${base}.bam" -o "${base}_sorted.bam"
if [ $? -ne 0 ]; then
    echo "Error: BAM sorting failed"
    exit 1
fi

samtools index "${base}_sorted.bam"
if [ $? -ne 0 ]; then
    echo "Error: BAM indexing failed"
    exit 1
fi

rm "${base}.sam"

# Create output directory for ASMv1.0
mkdir -p ASM_14test

# Convert coordinate sorted BAM to name sorted BAM
echo "Sorting BAM by name for ${base}"  # Print message before name sorting
samtools sort -n -o "${base}_namesorted.bam" "${base}_sorted.bam"
if [ $? -ne 0 ]; then
    echo "Error: Name sorting of BAM failed"
    exit 1
fi

# Run GATK AnalyzeSaturationMutagenesis with name sorted BAM
echo "Running GATK AnalyzeSaturationMutagenesis for ${base}"  
/home/s/sp877/miniconda3/envs/bioenv/share/gatk4-4.5.0.0-0/gatk AnalyzeSaturationMutagenesis -I "${base}_namesorted.bam" -R SgeneRBD.fasta --orf 106-846 -O ASM_14test/"${base}"
if [ $? -ne 0 ]; then
    echo "Error: GATK AnalyzeSaturationMutagenesis failed"
    exit 1
fi

echo "Analysis complete. Results are in ASM_14test/${base}"  
