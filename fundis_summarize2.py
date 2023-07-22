# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 2023

@author: ian.michael.bollinger@gmail.com with the support of ChatGPT 4.0

This is a module intended to be used as a part of a pipeline.

This can be used individually by calling the command:
    python /path/to/fundis_summarize.py -s /path/to/source_folder -p 80
    python /path/to/fundis_summarize.py --source_folder /path/to/source_folder --percent_system_use 80
"""
import os
import re
import glob
import csv
import shutil
from tqdm import tqdm
import argparse
import math
import multiprocessing
from pathlib import Path

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

# Function to update the name of a given fasta file to include the total supporting reads 
def update_fasta_file(file, parent_folder):
    with open(file, 'r') as f:
        lines = f.readlines()

    first_line = lines.pop(0)
    length = re.search(r'total_supporting_reads_\d+:\d\.\d-(\d+)\.\d', first_line).group(1)
    reads = re.search(r'total_supporting_reads_(\d+):\d\.\d-\d+\.\d', first_line).group(1)

    lines.insert(0, f'>{parent_folder}\n')

    with open(file, 'w') as f:
        f.writelines(lines)

    return length, reads


def get_fasta_reads(fasta_file):
    with open(fasta_file, 'r') as f:
        first_line = f.readline()
        reads = re.search(r'total_supporting_reads_(\d+):\d\.\d-\d+\.\d', first_line).group(1)
    return reads


def process_medaka_folder(folder, medaka_id):
    fasta_file = os.path.join(folder, 'consensus.fasta')

    if os.path.exists(fasta_file):
        p = Path(fasta_file)
        parent_folder = p.parent.parent.name
        parent_folder_name = re.sub(r'^sample_', '', parent_folder)
        reads = get_fasta_reads(fasta_file)
        target_fasta_file = os.path.join(summary_folder, f'{parent_folder_name}-{medaka_id}-RiC{reads}.fasta')
        shutil.copy(fasta_file, target_fasta_file)
        length, reads = update_fasta_file(target_fasta_file, f'{parent_folder_name}-{medaka_id}')

        p = Path(folder)
        medaka_parent_folder = str(p.parent)
        fastq_file_index = 1
        for fastq_file in sorted(glob.glob(f'{medaka_parent_folder}/reads_to_consensus_[0-9]*.fastq')):
            if os.path.exists(fastq_file):
                target_fastq_file = os.path.join(summary_folder, 'FASTQ Files', f'{parent_folder_name}-{fastq_file_index}.fastq')
                shutil.copy(fastq_file, target_fastq_file)
                fastq_file_index += 1
    
        return f'{parent_folder_name}-{medaka_id}', length, reads
    else:
        return None, None, None

def process_folder(folder):
    medaka_id = 1
    stats = []
    items = sorted(glob.iglob(folder + '/*', recursive=False))
    for item in items:
        if not os.path.isfile(item):
            p = Path(item)
            folder_name = p.name
            if folder_name.startswith('medaka'):
                filename, length, reads = process_medaka_folder(item, medaka_id)
                stats.append([filename, length, reads, medaka_id])
                medaka_id += 1
    return stats

def merge_fasta_files(folder, output_file):
    with open(os.path.join(folder, output_file), 'w') as f_out:
        for fasta_file in sorted(glob.glob(f'{folder}/*.fasta')):
            with open(fasta_file, 'r') as f:
                f_out.write(f.read())

def generate_summary(stats, writer):

    summary_totals = {'Total Unique Samples': 0,
                      'Total Consensus Sequences': 0,
                      'Total Reads in Consensus Sequences': 0}
    additional_stats = []
    counted_samples = {}

    for row in stats:
        filename = row[0]
        if filename is not None:
            if row[2] is not None:
                reads = int(row[2])
            else:
                reads = 0
            if filename.count('ONT') > 1:
                additional_stats.append(row)
            else:
                writer.writerow(row)
                unique_part = re.search(r'^([^-]+)', filename).group(1)
                if unique_part not in counted_samples:
                    summary_totals['Total Unique Samples'] += 1
                    counted_samples[unique_part] = 1
            summary_totals['Total Consensus Sequences'] += 1
            summary_totals['Total Reads in Consensus Sequences'] += reads

    writer.writerow([])
    writer.writerow(['Total Unique Samples', summary_totals['Total Unique Samples']])
    writer.writerow(['Total Consensus Sequences', summary_totals['Total Consensus Sequences']])
    writer.writerow(['Total Reads in Consensus Sequences', summary_totals['Total Reads in Consensus Sequences']])
    writer.writerow([])
    for row in additional_stats:
        writer.writerow(row)

def main(args):
    percent_system_use = float(args.percent_system_use/100) if args.percent_system_use else 0.8
    source_folder = args.source_folder if args.source_folder else os.path.dirname(os.path.realpath(__file__))
    
    # Name summary folder according to original source_folder but adding _Summary
    source_folder_split = '/'.join(source_folder.split('/')[:-1])
    source_folder_base = source_folder.split('/')[-1]    
    summary_folder = f"{source_folder_split}/{source_folder_base}_Summary"
    
    if not os.path.exists(summary_folder):
        os.mkdir(summary_folder)
    if not os.path.exists(os.path.join(summary_folder, 'FASTQ Files')):
        os.mkdir(os.path.join(summary_folder, 'FASTQ Files'))
    
    with open(os.path.join(summary_folder, 'summary.txt'), 'w', encoding='utf-8') as f_out:
        writer = csv.writer(f_out, delimiter='\t', lineterminator='\n')
        writer.writerow(['Filename', 'Length', 'Reads in Consensus', 'Multiple'])
        
        try:
            with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
                folders = sorted(glob.iglob(source_folder + '/*', recursive=False))
                results = [pool.apply_async(process_folder, args=(folder,)) for folder in folders]
                stats = [p.get() for p in results]
        except Exception as e:
            print(f"An error occurred during multiprocessing: {e}")
            raise  # re-raise the exception after handling it
        else:
            generate_summary(stats, writer)
    
    merge_fasta_files(summary_folder, 'summary.fasta')
    print('PASS: Successfully summarized NGSpeciesID folder\n')

if __name__ == "__main__":
    environment_dir = ""
    environment_cmd_prefix = ""
    environment_dir = check_os()
    main_working_dir = os.getcwd()
    
    # Parse user arguments
    parser = argparse.ArgumentParser(description="Process NGSpeciesID source folder.")
    parser.add_argument('-s','--source_folder', type=str, help='Path to the source folder')
    parser.add_argument('-p','--percent_system_use', type=str, help='Percent system use written as integer.')
    args = parser.parse_args()    
    main(args)