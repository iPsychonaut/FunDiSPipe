# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 2023

@author: ian.michael.bollinger@gmail.com with the support of ChatGPT 4.0

This is a module intended to be used as a part of a pipeline.

This can be used individually by calling the command: 
    python /path/to/fundis_haplotype_phaser.py -i /path/to/input_dir -p 80
    python /path/to/fundis_haplotype_phaser.py --input_dir /path/to/input_dir --percent_system_use 80
"""

import os
import glob
import subprocess
import pandas as pd
import multiprocessing
import math
import queue
import pysam
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Function to determine which Operating System the code is being executed in
def check_os():
    global environment_dir
    global environment_cmd_prefix
    # Determine the operating system in use
    # os.name will return 'posix', 'nt', or 'java' 
    os_name = os.name  
    # platform.system() will return 'Linux', 'Windows', 'Java', etc.
    platform_system = platform.system()  
    
    # If the operating system is Windows (identified by 'nt' from os.name or 'Windows' from platform.system())
    if os_name == 'nt' or platform_system == 'Windows':
        # Set the working directory to "E:" (or whatever drive letter is appropriate for your Windows system)
        environment_dir = "E:"
        # For running Linux commands in Windows Subsystem for Linux (WSL), prefix the command with "wsl " 
        environment_cmd_prefix = "wsl "
    
    # If the operating system is Linux (identified by 'posix' from os.name or 'Linux' from platform.system())
    elif os_name == 'posix' or platform_system == 'Linux':
        # Set the working directory to "/mnt/e" (or whatever the corresponding path is in your Linux system)
        environment_dir = "/mnt/e"  
    
    else:
        # If the operating system is neither Windows nor Linux, raise an Exception
        raise Exception("ERROR: OS NOT TESTED WITH THIS CODE.")
    
    # Print out the detected operating system and the determined environment directory
    print(f'Operating System: {platform_system}')
    print(f'Environment Directory: {environment_dir}')
    return environment_dir

# Function to extract sequence id and supporting reads count from a fasta file
def extract_data_from_fasta(fasta_file):
    data = []
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        sequence_id = seq_record.id
        supporting_reads_count = int(sequence_id.split('_supporting_reads_')[1])
        data.append((sequence_id, supporting_reads_count))
    return data

# Main function to process a folder
def analyze_ngspeciesid_folder(ngspeciesid_dir):
    # Find all fasta files in the input directory
    fasta_files = glob.glob(os.path.join(ngspeciesid_dir, "*.fasta"))
    
    # If no fasta files, return
    if not fasta_files:
        return
    
    # # Extract data from all fasta files
    # data = []
    # for fasta_file in fasta_files:
    #     data.extend(extract_data_from_fasta(fasta_file))
    
    # # Create a dataframe from the extracted data
    # df = pd.DataFrame(data, columns=['sequence_id', 'reads_count'])
    
    # # Find the index of the sequence with the maximum reads
    # index_max_reads = df['reads_count'].idxmax()
    
    # # Get the input fasta file
    # input_fasta = fasta_files[index_max_reads]
    
    # Run through each file 
    for input_fastq in fasta_files:
        if os.path.isdir(base_folder):
            print(f"PASS: Skipping haplotype split, {base_folder} already exists\n")
        else:
            # # Define the base folder
            base_folder = input_fasta.rsplit('.', 1)[0]
            os.makedirs(base_folder, exist_ok=True)  # add exist_ok=True
            os.chdir(base_folder)

            input_fastq = input_fasta.replace('consensus_reference','reads_to_consensus').replace('.fasta','.fastq')
            
            fasta_basename = os.path.basename(input_fasta)
            
            sam_file = os.path.join(base_folder, fasta_basename.replace('.fasta','.sam'))
            
            bam_file = os.path.join(base_folder, sam_file.replace('.sam','.bam'))
            
            sorted_bam = os.path.join(base_folder, bam_file.replace('.bam','_sorted.bam'))
                                 
            bcf_file = os.path.join(base_folder, sorted_bam.replace('_sorted.bam','.bcf'))
            
            pileup_output = os.path.join(base_folder, bcf_file.replace(".bcf", "_pileup.bcf"))
            
            vcf_file = os.path.join(base_folder, bcf_file.replace(".bcf",'.vcf'))
            
            phased_vcf = os.path.join(base_folder, vcf_file.replace('.vcf','.phased.vcf'))
            
            haplotype_fasta = os.path.join(base_folder, phased_vcf.replace('.phased.vcf','_haplotypes.fasta'))
            
            # Step 1: Use BWA to align your FASTQ file to your reference FASTA file
            # First, index the reference genome
            subprocess.run(["bwa", "index", input_fasta], check=True)
            
            # Then, align the reads to the reference
            subprocess.run(["bwa", "mem", input_fasta, input_fastq, "-o", sam_file], check=True)
            
            # Step 2: Convert SAM to BAM using SAMtools
            subprocess.run(["samtools", "view", "-S", "-b", sam_file, "-o", bam_file], check=True)
            
            # Step 3: Sort the BAM file
            subprocess.run(["samtools", "sort", bam_file, "-o", sorted_bam], check=True)
            
            # Add this line to create the index file for the sorted BAM file:
            subprocess.run(["samtools", "index", sorted_bam], check=True)
            
            # Step 4: Call variants with bcftools to get a VCF file
            
            subprocess.run(["bcftools", "mpileup", "-Ou", "-f", input_fasta, sorted_bam, "-o", pileup_output], check=True)
            subprocess.run(["bcftools", "call", "-mv", "-Ob", "-o", bcf_file, pileup_output], check=True)
            subprocess.run(["bcftools", "view", bcf_file, "-Ov", "-o", vcf_file], check=True)
            
            # Step 5: Use Whatshap to find haplotypes
            subprocess.run(["whatshap", "phase", "-o", phased_vcf, "--reference", input_fasta, vcf_file, sorted_bam], check=True)
            
            # Step 6: Parse the output VCF file and generate sequences based on the original input FASTA        
            # First, let's load the original reference FASTA file
            original_seq_record = SeqIO.read(input_fasta, "fasta")
            original_seq = str(original_seq_record.seq)
            
            # Then, let's parse the VCF file using pysam
            vcf = pysam.VariantFile(phased_vcf)
            
            new_sequences = []
            
            # Let's iterate over each variant in the VCF file
            for variant in vcf:
                # We'll create a new sequence for each allele
                for allele in variant.alleles:
                    if allele != variant.ref:  # We don't need to create a new sequence for the reference allele
                        new_seq_str = original_seq[:variant.pos-1] + allele + original_seq[variant.pos:]
                        record_id = variant.id if variant.id is not None else f"pos_{variant.pos}_allele_{allele}"
                        new_seq_record = SeqRecord(Seq(new_seq_str), id=record_id, description="")
                        new_sequences.append(new_seq_record)
            
            # Finally, let's write the new sequences to a new FASTA file
            SeqIO.write(new_sequences, haplotype_fasta, "fasta")
            
            # Initialize an empty list to store all haplotype sequences
            haplotype_sequences = []
        
            # Use glob to find all haplotype fasta files
            haplotype_files = glob.glob(os.path.join(base_folder, '**', '*haplotypes.fasta'), recursive=True)
        
            # Take Haplotype File and save individual fasta files for each haplotype that consist of a single degenerate sequence
            for haplotype_file in haplotype_files:
                for i, record in enumerate(SeqIO.parse(haplotype_file, "fasta")):
                    # Create a new fascondata file for each sequence
                    individual_fasta_path = os.path.join(os.path.dirname(haplotype_file), f"{record.id}_individual.fasta")
                    with open(individual_fasta_path, "w") as individual_fasta_file:
                        SeqIO.write(record, individual_fasta_file, "fasta")
                            
        # Remove files that end in '.bcf', '_sorted.vcf', '.amb', '.ann', '.bwt', '.pac', '.sa', and '.bai' from the base folder
        exts = ['.bcf', '_sorted.vcf', '.amb', '.ann', '.bwt', '.pac', '.sa', '.bai']
        
        # Iterate over each extension
        for ext in exts:
            # Use glob to find all files with the current extension, including subfolders
            file_list = glob.glob(os.path.join(ngspeciesid_dir, '**', '*' + ext), recursive=True)
        
            # Iterate over each file and remove it
            for file_path in file_list:
                os.remove(file_path)
        return True
 
# Main entry point of the script
def main(args):
    # Define the main directory based on the argument or use a default value
    input_dir = args.input_dir if args.input_dir else '/mnt/e/Fundis/NGSpeciesID-20230719T211158Z-001/NGSpeciesID'
    percent_system_use = float(args.percent_system_use/100) if args.percent_system_use else 0.8
    
    # Get the number of CPUs and calculate the number of threads
    num_cpus = multiprocessing.cpu_count()
    cpu_threads = int(math.floor(num_cpus * percent_system_use))
    
    # Get a list of all folders in the main directory, excluding folders named '__Summary__' or 'Summary2'
    folder_list = [f.path for f in os.scandir(input_dir) if f.is_dir() and '__Summary__' not in f.path and 'Summary2' not in f.path]
    
    try:
        # Initialize a multiprocessing pool and process all folders
        with multiprocessing.Pool(cpu_threads) as pool:
            pool.map(analyze_ngspeciesid_folder, folder_list)
    except queue.Empty:
        print("Queue is empty. All Pipeline tasks have been processed.")
    finally:
        print('PASS: Haplotypes phased for all RiC >=9 conesensus fasta for each sample in NGSPeciesID input_dir\n')

# If this script is the main entry point, parse the arguments and call the main function
if __name__ == "__main__":
    # Global environment_dir
    environment_dir = ""
    environment_cmd_prefix = ""
    environment_dir = check_os()
    main_working_dir = os.getcwd()
    
    # Parse user arguments
    parser = argparse.ArgumentParser(description="Process NGSpeciesID main folder.")
    parser.add_argument('-i','--input_dir', type=str, help='Path to the NGSpeciesID main directory')
    parser.add_argument('-p','--percent_system_use', type=str, help='Percent system use written as integer.')
    args = parser.parse_args()
    main(args)