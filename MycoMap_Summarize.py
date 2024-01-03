# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 12:23:15 2023

@author: ian.michael.bollinger@gmail.com/researchconsultants@critical.consulting

This is the MycoMap Summarizer (summarize.py) developed for a modified protocol by

Stephen Douglas Russell and paid for by the Fungal Diversity Survey (FunDiS).

Protocol Link: https://www.protocols.io/view/primary-data-analysis-basecalling-demultiplexing-a-dm6gpbm88lzp/v3?step=3
"""

# Base Python Imports
import subprocess, os, re, glob, csv, shutil
from pathlib import Path

# Required Python Imports
from termcolor import cprint
from datetime import datetime

# Function to color coded print to console and save to log_file information
def log_print(input_message, log_file=None):
    """
    Logs a message to a file and prints it to the console with appropriate coloring.
    
    This function takes a message and logs it to the specified file. Additionally, the message is printed to the 
    console, potentially with specific coloring depending on the context.
    
    Parameters:
        input_message (str): Message to be logged and printed.
        log_file (str): Path to the log file.

    Notes:
        - The function uses a global default log file if none is specified.
        - Timestamps each log entry for easy tracking.
        - Utilizes color coding in the console to distinguish between different types of messages (e.g., errors, warnings).
        - Supports color coding for specific message types: NOTE, CMD, ERROR, WARN, and PASS.
        - Falls back to default (white) color if the message type is unrecognized.
    """
    # Access the global variable
    global DEFAULT_LOG_FILE
    
    # Use the default log file if none specified
    if log_file is None:
        log_file = DEFAULT_LOG_FILE  
    
    # Establish current date-time
    now = datetime.now()
    message = f'[{now:%Y-%m-%d %H:%M:%S}]\t{input_message}'

    # Determine the print color based on the input_message content
    message_type_dict = {'NOTE': ['blue'],
                         'CMD': ['cyan'],
                         'ERROR': ['red'],
                         'WARN': ['yellow'],
                         'PASS': ['green'],}
    print_color = ['white']  # Default color
    for key, value in message_type_dict.items():
        if key.lower() in input_message.lower():
            print_color = value
            break

    # Writing the message to the log file
    with open(log_file, 'a') as file:
        print(message, file=file)

    # Handling different message types for colored printing
    try:
        cprint(message, print_color[0])
    except (KeyError, IndexError):
        cprint(message, print_color[1] if len(print_color) > 1 else 'white')



def run_subprocess_cmd(cmd_list, shell_check):
    if isinstance(cmd_list, str):
        log_print(f"CMD:\t{cmd_list}")    
        process = subprocess.run(cmd_list, text=True, shell=shell_check, capture_output=True)
        if process.returncode != 0:
            log_print(f"ERROR:\t{process.stderr}")
        else:
            log_print(f"PASS:\tSuccessfully processed command: {cmd_list}")
    else:        
        log_print(f"CMD:\t{' '.join(cmd_list)}")    
        process = subprocess.run(cmd_list, text=True, shell=shell_check, capture_output=True)
        if process.returncode != 0:
            log_print(f"ERROR:\t{process.stderr}")
        else:
            log_print(f"PASS:\tSuccessfully processed command: {' '.join(cmd_list)}")

def update_fasta_file(file, parent_folder):
    lines = []
    with open(file, 'r') as f:
        lines = f.readlines()
    first_line = lines.pop(0)
    length = re.search(r'total_supporting_reads_\d+:\d\.\d-(\d+)\.\d', first_line).group(1)
    reads = re.search(r'total_supporting_reads_(\d+):\d\.\d-\d+\.\d', first_line).group(1)
    lines.insert(0, '>' + parent_folder + '\n')
    with open(file, 'w') as f:
        [f.write(line) for line in lines]

    return length, reads

def get_fasta_reads(fasta_file):
    with open(fasta_file, 'r') as f:
        first_line = f.readline()
        reads = re.search(r'total_supporting_reads_(\d+):\d\.\d-\d+\.\d', first_line).group(1)
    return reads

def process_medaka_folder(folder, medaka_id, summary_folder):
    print(f'Processing Medaka folder: {folder}...')
    fasta_file = os.path.join(folder, 'consensus.fasta')
    if os.path.exists(fasta_file):
        p = Path(fasta_file)
        parent_folder = p.parent.parent.name
        parent_folder_name = re.sub(r'^sample_', '', parent_folder)
        reads = get_fasta_reads(fasta_file)
        target_fasta_file = os.path.join(summary_folder, f'{parent_folder_name}-{medaka_id}-RiC{reads}.fasta')
        # shutil.move(fasta_file, target_fasta_file) # TODO
        shutil.copy(fasta_file, target_fasta_file)
        length, reads = update_fasta_file(target_fasta_file, f'{parent_folder_name}-{medaka_id}')

    p = Path(folder)
    medaka_parent_folder = str(p.parent)
    fastq_file_index = 1
    for fastq_file in sorted(glob.glob(f'{medaka_parent_folder}/reads_to_consensus_[0-9]*.fastq')):
        if os.path.exists(fastq_file):
            target_fastq_file = os.path.join(summary_folder, 'FASTQ Files', f'{parent_folder_name}-{fastq_file_index}.fastq')
            # shutil.move(fastq_file, target_fastq_file) # TODO
            shutil.copy(fastq_file, target_fastq_file)
            fastq_file_index += 1

    return f'{parent_folder_name}-{medaka_id}', length, reads

def process_folder(folder, summary_folder):
    medaka_id = 1
    stats = []
    for item in sorted(glob.iglob(folder + '/*', recursive=False)):
        if not os.path.isfile(item):
            p = Path(item)
            folder_name = p.name
            if folder_name.startswith('medaka'):
                filename, length, reads = process_medaka_folder(item, medaka_id, summary_folder)
                stats.append([filename, length, reads, medaka_id])
                medaka_id += 1

    return stats

def merge_fasta_files(folder, output_file):
    with open(os.path.join(folder, output_file), 'w') as f_out:
        for fasta_file in sorted(glob.glob(f'{folder}/*.fasta')):
            with open(fasta_file, 'r') as f:
                f_out.write(f.read())

def generate_summary(stats, writer):
    # Set placeholders for stats
    summary_totals = {'Total Unique Samples': 0,
                      'Total Consensus Sequences': 0,
                      'Total Reads in Consensus Sequences': 0,}
    additional_stats = []
    counted_samples = {}

    # Compile stats
    for row in stats:
        filename = row[0]
        reads = int(row[2])
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

    # Create final TSV summary
    writer.writerow([])
    writer.writerow(['Total Unique Samples', summary_totals['Total Unique Samples']])
    writer.writerow(['Total Consensus Sequences', summary_totals['Total Consensus Sequences']])
    writer.writerow(['Total Reads in Consensus Sequences', summary_totals['Total Reads in Consensus Sequences']])
    writer.writerow([])
    for row in additional_stats:
        writer.writerow(row)

# Main Summary Function
def summary_main(source_folder, hap_phase_var):
    if hap_phase_var == 0:
        print("Running Default Command Line Version of MycoMap Summarize...")
    
    elif hap_phase_var == 1:
        print("Running Haplotype Phased Version of MycoMap Summarize...")
        # Create Summary folder and subfolder
        summary_folder = os.path.join(source_folder, '__Summary__')
        if not os.path.exists(summary_folder):
            os.mkdir(summary_folder)
        if not os.path.exists(os.path.join(summary_folder, 'FASTQ Files')):
            os.mkdir(os.path.join(summary_folder, 'FASTQ Files'))
        
        # Start Generating a TSV to contain summary data
        with open(os.path.join(summary_folder, 'summary.txt'), 'w', encoding='utf-8') as f_out:
            writer = csv.writer(f_out, delimiter='\t', lineterminator='\n')
            writer.writerow(['Filename', 'Length', 'Reads in Consensus', 'Multiple'])
        
            # Create a variable to contain folder stats
            stats = []
        
            # For each directory in the source_folder generate the stats for that folder and add it to the CSV
            for item in sorted(glob.iglob(source_folder + '/*', recursive=False)):
                if not os.path.isfile(item):
                    folder_stats = process_folder(item, summary_folder)
                    stats = stats + folder_stats
        
            # Generate Final Summary based on the CSV and the compiled Stats
            generate_summary(stats, writer)
        
        merge_fasta_files(summary_folder, 'summary.fasta')
    print('Done!')

# Wrapper function to run ngsid_prep in a separate thread
def run_summary_prep(ngsid_output_dir, hap_phase_bool):
    """
    Wrapper function to run summary preparation in a separate thread.

    This function initiates the summary preparation process, potentially for generating summaries of data processing. 
    It handles the execution in a separate thread and logs the start of the process and any errors encountered.

    Global Variables:
        output_area (str): A global variable used for logging output messages.

    Notes:
        - The actual summary preparation logic (e.g., summary_main) needs to be implemented or called within this function.
        - Catches and logs exceptions that may occur during summary preparation.
    """
    log_print( "MycoMap Summary preparation started...\n")
    
    try:
        summary_main(ngsid_output_dir, hap_phase_bool)
    except Exception as e:
        log_print( f"ERROR:\t{str(e)}")

if __name__ == "__main__":
    # Get the source folder
    source_folder = os.path.dirname(os.path.realpath(__file__))
    summary_main(source_folder)
