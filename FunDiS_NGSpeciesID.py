# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 12:23:15 2023

@author: ian.michael.bollinger@gmail.com/researchconsultants@critical.consulting

This is the Next-Generation Barcoding Species Identifier (NGSpeciesID) main wrapper developed for a modified

protocol by Stephen Douglas Russell and paid for by the Fungal Diversity Survey (FunDiS).

Protocol Link: https://www.protocols.io/view/primary-data-analysis-basecalling-demultiplexing-a-dm6gpbm88lzp/v3?step=3
"""

# Base Python Imports
import subprocess, os

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


from FunDiS_hap_phase import phase_haplotypes

# Wrapper function to run ngsid_prep in a separate thread
def run_ngsid_prep(ngsid_primers_path, input_fastq_file, sample_size, min_length_bp, max_std_dev_bp, hap_phase_bool, ngsid_output_dir):
    """
    Wrapper function to run NGSpeciesID preparation in a separate thread.

    This function handles the NGSpeciesID preparation process. It executes the preparation steps in a separate thread 
    and logs any errors encountered during execution.

    Parameters:
        ngsid_primers_path (str): Path to the NGSpeciesID primers file.
        input_fastq_file (str): Path to the input FASTQ file.
        sample_size (int): The sample size parameter for NGSpeciesID.
        min_length_bp (int): Minimum length in base pairs for NGSpeciesID.
        max_std_dev_bp (int): Maximum standard deviation in base pairs for NGSpeciesID.
        hap_phase_bool (bool): Boolean to indicate if haplotype phasing is to be performed.
        ngsid_output_dir (str): Directory path for storing NGSpeciesID output files.

    Global Variables:
        output_area (str): A global variable used for logging output messages.

    Notes:
        - Captures and logs exceptions during NGSpeciesID preparation.
    """
    global DEFAULT_LOG_FILE
    try:
        ngsid_prep(ngsid_primers_path, input_fastq_file, sample_size, min_length_bp, max_std_dev_bp, hap_phase_bool, ngsid_output_dir)
    except Exception as e:
        log_print( f"ERROR:\t{str(e)}")

# Function to perform NGSpeciesID preparation
def ngsid_prep(ngsid_primers_path, input_fastq_file, sample_size, min_length_bp, max_std_dev_bp, hap_phase_bool, ngsid_output_dir):
    """
    Performs NGSpeciesID preparation.

    This function carries out the necessary steps for NGSpeciesID preparation. It processes FASTQ files, 
    executes NGSpeciesID, and optionally performs haplotype phasing.

    Parameters:
        ngsid_primers_path (str): Path to the NGSpeciesID primers file.
        input_fastq_file (str): Path to the input FASTQ file.
        sample_size (int): The sample size parameter for NGSpeciesID.
        min_length_bp (int): Minimum length in base pairs for NGSpeciesID.
        max_std_dev_bp (int): Maximum standard deviation in base pairs for NGSpeciesID.
        hap_phase_bool (bool): Boolean to indicate if haplotype phasing is to be performed.
        ngsid_output_dir (str): Directory path for storing NGSpeciesID output files.

    Global Variables:
        output_area (str): A global variable used for logging output messages.
        cpu_threads (int): Number of CPU threads to be used in processing.

    Notes:
        - Processes each .fastq file in the specified output directory.
        - Logs progress and errors during the preparation process.
    """    
    global cpu_threads, DEFAULT_LOG_FILE
    # For any .fastq file in the ngsid_output_dir run the NGSID function on it
    fastq_files_list = []
    
    # Iterate over all files in the given folder
    for file in os.listdir(ngsid_output_dir):
        # Check if the file ends with .fastq
        if file.endswith('.fastq'):
            # Add the file to the list, with full path
            fastq_files_list.append(os.path.join(ngsid_output_dir, file))
 
    for input_fastq_file in fastq_files_list:
        ngsid_output_dir = ngsid_fastq(ngsid_primers_path, input_fastq_file, sample_size, min_length_bp, max_std_dev_bp)
        if hap_phase_bool == True:
            try:
                phased_fasta_file = phase_haplotypes(ngsid_output_dir, cpu_threads)
            except TypeError:
                log_print( f"Skipping Haplotype Phasing, unable to process: {input_fastq_file} with NGSpeciesID")
                continue
        log_print(f"PASS:\tSuccessfully processed {ngsid_output_dir}")

def ngsid_fastq(ngsid_primers_path, input_fastq_file, sample_size, min_length_bp, max_std_dev_bp):
    """
    Processes a FASTQ file with NGSpeciesID.

    This function handles the processing of a single FASTQ file using the NGSpeciesID tool. It constructs and executes 
    the NGSpeciesID command with the provided parameters.

    Parameters:
        ngsid_primers_path (str): Path to the NGSpeciesID primers file.
        input_fastq_file (str): Path to the input FASTQ file.
        sample_size (int): The sample size parameter for NGSpeciesID.
        min_length_bp (int): Minimum length in base pairs for NGSpeciesID.
        max_std_dev_bp (int): Maximum standard deviation in base pairs for NGSpeciesID.

    Global Variables:
        output_area (str): A global variable used for logging output messages.
        cpu_threads (int): Number of CPU threads to be used in processing.

    Returns:
        str: The directory path where the NGSpeciesID output is stored.

    Notes:
        - Logs the constructed command and the outcome of the process.
        - Handles and logs medaka consensus calls.
    """
    global cpu_threads, DEFAULT_LOG_FILE
    log_print( f"Processing {input_fastq_file} with NGSpeciesID...")
    ngsid_output_dir = input_fastq_file.replace(".fastq","")
    
    NGSID_cmd = ["NGSpeciesID", "--ont", "--consensus", 
                   "--sample_size", str(sample_size),
                   "--m", str(min_length_bp),
                   "--s", str(max_std_dev_bp),
                   "--racon", "--racon_iter", str(3),
                   "--primer_file", ngsid_primers_path,
                   "--fastq", input_fastq_file,
                   "--outfolder", ngsid_output_dir]
    
    log_print( f"CMD:\t{' '.join(NGSID_cmd)}")
    process = subprocess.run(NGSID_cmd, text=True)
    if process.returncode != 0:
        log_print( f"ERROR:\t{process.stderr}")
        return None
    else:
        log_print( f"PASS:\tSuccessfully processed {input_fastq_file} with NGSpeciesID")
        
    racon_folders = [folder for folder in os.listdir(ngsid_output_dir)
                     if os.path.isdir(os.path.join(ngsid_output_dir, folder)) and 'racon' in folder]
    
    for folder in racon_folders:
        base_number = folder.split("_")[-1]
        medaka_output_dir = f"{ngsid_output_dir}/medaka_cl_id_{base_number}"

        if not os.path.exists(medaka_output_dir):
            os.makedirs(medaka_output_dir)
            
        reads_to_consensus_fastq = f"{ngsid_output_dir}/reads_to_consensus_{base_number}.fastq"
        consensus_ref_fasta = f"{ngsid_output_dir}/consensus_reference_{base_number}.fasta"        
        call_medaka_script(reads_to_consensus_fastq, consensus_ref_fasta, medaka_output_dir)
    
    return ngsid_output_dir

def call_medaka_script(reads_to_consensus_fastq, consensus_ref_fasta, medaka_output_dir):
    """
    Calls a script to run Medaka on generated consensus reads.

    This function handles the execution of a shell script to run Medaka, a tool for generating consensus sequences from 
    sequencing data. It logs the process and captures the output and errors of the script execution.

    Parameters:
        reads_to_consensus_fastq (str): Path to the FASTQ file containing reads for consensus.
        consensus_ref_fasta (str): Path to the FASTA file used as a reference for consensus calling.
        medaka_output_dir (str): Directory path for storing Medaka output files.

    Global Variables:
        output_area (str): A global variable used for logging output messages.
        cpu_threads (int): Number of CPU threads to be used in processing.

    Notes:
        - Executes a shell script for running Medaka and handles the subprocess communication.
        - Logs the output and errors from the Medaka script execution.
    """
    global cpu_threads, DEFAULT_LOG_FILE
    log_print(f"Running medaka on generated consensus reads: {reads_to_consensus_fastq}...")
    medaka_shell_path = "/mnt/d/FunDiS/run_medaka.sh"
    script_command = [medaka_shell_path, reads_to_consensus_fastq, consensus_ref_fasta, medaka_output_dir, str(cpu_threads)]

    process = subprocess.Popen(script_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate()

    if stdout:
        log_print("NOTE:\n" + stdout)
    if stderr:
        log_print("ERROR:\n" + stderr)

    if process.returncode != 0:
        log_print("ERROR:\tSubprocess returned with error")
    else:
        log_print("PASS:\tSuccessfully ran medaka on generated consensus reads")
        
if  __name__ == "__main__":
    print("UNLOGGED DEBUG:\tUNDEVELOPED FunDiS_NGSpeciesID.py DEBUG AREA")
