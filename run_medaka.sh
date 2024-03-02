#!/bin/bash

# Find and source the conda.sh script
CONDA_SH_PATH=$(find ~/ -name conda.sh 2>/dev/null | head -n 1)
if [ -n "$CONDA_SH_PATH" ]; then
    source "$CONDA_SH_PATH"
else
    echo "conda.sh not found, please ensure Conda is installed."
    exit 1
fi

# Activate the medaka conda environment
conda init medaka
conda activate medaka

# Define your variables
reads_to_consensus_fastq="$1"
consensus_ref_fasta="$2"
medaka_output_dir="$3"
cpu_threads="$4"

# Print the command for debugging
echo "Running medaka_consensus with:"
echo "Input FASTQ: $reads_to_consensus_fastq"
echo "Reference FASTA: $consensus_ref_fasta"
echo "Output Directory: $medaka_output_dir"
echo "CPU Threads: $cpu_threads"

# Create the output directory if it doesn't exist
mkdir -p "$medaka_output_dir"

# Run the medaka_consensus command and capture its output
echo "Starting medaka_consensus..."
medaka_consensus -i "$reads_to_consensus_fastq" -d "$consensus_ref_fasta" -o "$medaka_output_dir" -t "$cpu_threads" 2>&1 | tee "${medaka_output_dir}/medaka_output.log"

# Check for errors or success
if [ $? -eq 0 ]; then
    echo "medaka_consensus completed successfully."
else
    echo "Error running medaka_consensus. Check the log file for details."
fi

# Deactivate the conda environment
conda deactivate