# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 12:23:15 2023

@author: ian.michael.bollinger@gmail.com/researchconsultants@critical.consulting

This is the Next-Generation Barcoding Species Identifier (NGSpeciesID) main wrapper developed for a modified

protocol by Stephen Douglas Russell and paid for by the Fungal Diversity Survey (FunDiS).

Protocol Link: https://www.protocols.io/view/primary-data-analysis-basecalling-demultiplexing-a-dm6gpbm88lzp/v3?step=3
"""

# Base Python Imports
import os

# Custom Pythyon Imports
from FunDiS_hap_phase import phase_haplotypes
from FunDiS_Tools import log_print, run_subprocess_cmd

# Global output_area variable
PERCENT_RESOURCES = 0.75

# Wrapper function to run ngsid_prep in a separate thread
def run_ngsid_prep(ngsid_primers_path, input_fastq_file, sample_size, min_length_bp, max_std_dev_bp, hap_phase_bool, ngsid_output_dir, summary_ric_cutoff, CPU_THREADS):
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
        summary_ric_cutoff (int): Minimum Reads in Consensus (RiC) allowed for processing

    Global Variables:
        DEFAULT_LOG_FILE (str): A global variable used for logging output messages.

    Notes:
        - Captures and logs exceptions during NGSpeciesID preparation.
    """
    global DEFAULT_LOG_FILE
    try:
        ngsid_prep(ngsid_primers_path, input_fastq_file, sample_size, min_length_bp, max_std_dev_bp, hap_phase_bool, ngsid_output_dir, summary_ric_cutoff, CPU_THREADS)
    except Exception as e:
        log_print( f"ERROR:\t{str(e)}")

# Function to perform NGSpeciesID preparation
def ngsid_prep(ngsid_primers_path, input_fastq_file, sample_size, min_length_bp, max_std_dev_bp, hap_phase_bool, ngsid_output_dir, summary_ric_cutoff, CPU_THREADS):
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
        summary_ric_cutoff (int): Minimum Reads in Consensus (RiC) allowed for processing

    Global Variables:
        DEFAULT_LOG_FILE (str): A global variable used for logging output messages.
        CPU_THREADS (int): Number of CPU threads to be used in processing.

    Notes:
        - Processes each .fastq file in the specified output directory.
        - Logs progress and errors during the preparation process.
    """    
    global DEFAULT_LOG_FILE
    # For any .fastq file in the ngsid_output_dir run the NGSID function on it
    fastq_files_list = []
    
    # Iterate over all files in the given folder
    for file in os.listdir(ngsid_output_dir):
        # Check if the file ends with .fastq
        if file.endswith('.fastq'):
            # Add the file to the list, with full path
            fastq_files_list.append(os.path.join(ngsid_output_dir, file))
 
    for input_fastq_file in fastq_files_list:
        ngsid_output_dir = ngsid_fastq(ngsid_primers_path, input_fastq_file, sample_size, min_length_bp, max_std_dev_bp, CPU_THREADS)
        if hap_phase_bool == True:
            try:
                phased_fasta_file = phase_haplotypes(ngsid_output_dir, summary_ric_cutoff, CPU_THREADS)
            except TypeError:
                log_print( f"Skipping Haplotype Phasing, unable to process: {input_fastq_file} with NGSpeciesID")
                continue
        log_print(f"PASS:\tSuccessfully processed {ngsid_output_dir}")

def check_racon_consensus_exists(ngsid_output_dir):
    """
    Check if each racon output directory contains a 'consensus.fasta' file.
    
    Parameters:
        ngsid_output_dir (str): Directory containing racon output folders.
    
    Returns:
        bool: True if every racon folder has a 'consensus.fasta' file, False otherwise.
    """
    racon_folders = [folder for folder in os.listdir(ngsid_output_dir)
                     if os.path.isdir(os.path.join(ngsid_output_dir, folder)) and 'racon' in folder]
    for folder in racon_folders:
        consensus_path = os.path.join(ngsid_output_dir, folder, "consensus.fasta")
        if not os.path.isfile(consensus_path):
            return False
    return True

def ngsid_fastq(ngsid_primers_path, input_fastq_file, sample_size, min_length_bp, max_std_dev_bp, CPU_THREADS):
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
        DEFAULT_LOG_FILE (str): A global variable used for logging output messages.

    Returns:
        str: The directory path where the NGSpeciesID output is stored.

    Notes:
        - Logs the constructed command and the outcome of the process.
        - Handles and logs medaka consensus calls.
    """
    global DEFAULT_LOG_FILE
    log_print(f"Processing {input_fastq_file} with NGSpeciesID...")
    ngsid_output_dir = input_fastq_file.replace(".fastq", "")
    
    # Start with the default racon iterations
    racon_iter = 3
    
    while racon_iter > 0:
        NGSID_cmd = ["NGSpeciesID", "--ont", "--consensus", 
                     "--sample_size", str(sample_size),
                     "--m", str(min_length_bp),
                     "--s", str(max_std_dev_bp),
                     "--racon", "--racon_iter", str(racon_iter),
                     "--primer_file", ngsid_primers_path,
                     "--fastq", input_fastq_file,
                     "--outfolder", ngsid_output_dir]
        
        run_subprocess_cmd(NGSID_cmd, False)
        
        if check_racon_consensus_exists(ngsid_output_dir):
            log_print("All racon folders contain a consensus.fasta file.")
            break  # Exit loop if all consensus files are found
        else:
            log_print(f"Missing consensus.fasta files in racon output. Reducing racon iterations to {racon_iter-1}.")
            racon_iter -= 1  # Reduce racon iterations and retry
    
    if racon_iter == 0:
        log_print("Failed to generate consensus.fasta files with racon polishing.")
        # Handle the failure case as needed
    
    racon_folders = [folder for folder in os.listdir(ngsid_output_dir)
                     if os.path.isdir(os.path.join(ngsid_output_dir, folder)) and 'racon' in folder]
    
    for folder in racon_folders:
        base_number = folder.split("_")[-1]
        medaka_output_dir = f"{ngsid_output_dir}/medaka_cl_id_{base_number}"

        if not os.path.exists(medaka_output_dir):
            os.makedirs(medaka_output_dir)
            
        reads_to_consensus_fastq = f"{ngsid_output_dir}/reads_to_consensus_{base_number}.fastq"
        consensus_ref_fasta = f"{ngsid_output_dir}/consensus_reference_{base_number}.fasta"        
        call_medaka_script(reads_to_consensus_fastq, consensus_ref_fasta, medaka_output_dir, CPU_THREADS)
    
    return ngsid_output_dir

def call_medaka_script(reads_to_consensus_fastq, consensus_ref_fasta, medaka_output_dir, CPU_THREADS):
    """
    Calls a script to run Medaka on generated consensus reads.

    This function handles the execution of a shell script to run Medaka, a tool for generating consensus sequences from 
    sequencing data. It logs the process and captures the output and errors of the script execution.

    Parameters:
        reads_to_consensus_fastq (str): Path to the FASTQ file containing reads for consensus.
        consensus_ref_fasta (str): Path to the FASTA file used as a reference for consensus calling.
        medaka_output_dir (str): Directory path for storing Medaka output files.

    Global Variables:
        DEFAULT_LOG_FILE (str): A global variable used for logging output messages.
        CPU_THREADS (int): Number of CPU threads to be used in processing.

    Notes:
        - Executes a shell script for running Medaka and handles the subprocess communication.
        - Logs the output and errors from the Medaka script execution.
    """
    global DEFAULT_LOG_FILE
    log_print(f"Running medaka on generated consensus reads: {reads_to_consensus_fastq}...")
    medaka_shell_path = "/mnt/d/FunDiS/run_medaka.sh"
    script_command = ["bash", medaka_shell_path, reads_to_consensus_fastq, consensus_ref_fasta, medaka_output_dir, str(CPU_THREADS)]
    run_subprocess_cmd(script_command, False)

if  __name__ == "__main__":
    print("UNLOGGED DEBUG:\tUNDEVELOPED FunDiS_NGSpeciesID.py DEBUG AREA")
