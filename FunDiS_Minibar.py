# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 12:23:15 2023

@author: ian.michael.bollinger@gmail.com/researchconsultants@critical.consulting

This is the Mini-Barcoder (minibar.py) main wrapper developed for a modified

protocol by Stephen Douglas Russell and paid for by the Fungal Diversity Survey (FunDiS).

Protocol Link: https://www.protocols.io/view/primary-data-analysis-basecalling-demultiplexing-a-dm6gpbm88lzp/v3?step=3
"""

# Base Python Imports
import os, gzip

# Custom Python Imports
from FunDiS_Tools import log_print, generate_log_file, initialize_logging_environment, run_subprocess_cmd

# Global output_area variable
CPU_THREADS = 1
PERCENT_RESOURCES = 0.75

# Function to unzip .fastq.gz file with gzip.
def unzip_fastq(ngsid_fastq_gz_path, ngsid_fastq_path):
    """
    Unzips a .fastq.gz file into a .fastq file.

    This function unzips a gzip compressed FASTQ file (.fastq.gz) into a regular FASTQ file (.fastq). 
    It reads from the compressed file and writes the uncompressed data to a new file.

    Parameters:
        ngsid_fastq_gz_path (str): The file path for the .fastq.gz file.
        ngsid_fastq_path (str): The target file path for the uncompressed .fastq file.
    """
    with gzip.open(ngsid_fastq_gz_path, 'rb') as f_in:
        with open(ngsid_fastq_path, 'wb') as f_out:
            f_out.write(f_in.read())

# Wrapper function to run ngsid_minibar_prep in a separate thread.
def run_minibar_prep(minibar_path, minibar_index_path, ngsid_fastq_gz_path, ngsid_output_dir, chopper_command_dict, minibar_command_dict):
    """
    Wrapper function to run minibar preparation in a separate thread.
    
    This function is designed to execute the minibar preparation process in an isolated thread. 
    It handles exceptions and logs any errors encountered during the execution.
    
    Parameters:
        minibar_path (str): Path to the minibar script.
        minibar_index_path (str): Path to the minibar index file.
        ngsid_fastq_gz_path (str): Path to the compressed FASTQ file.
        ngsid_output_dir (str): Directory path for storing output files.
        chopper_command_dict (dict): Dictionary containing chopper settings.
        minibar_command_dict (dict): Dictionary containing MiniBar settings.
    
    Notes:
        - Logs an error message if an exception is encountered during the minibar preparation.
    """
    global DEFAULT_LOG_FILE
    try:
        minibar_prep(minibar_path, minibar_index_path, ngsid_fastq_gz_path, ngsid_output_dir, chopper_command_dict, minibar_command_dict)
    except Exception as e:
        log_print( f"ERROR:\t{str(e)}")

# Function to perform minibar preparation  
def minibar_prep(minibar_path, minibar_index_path, input_file_path, ngsid_output_dir, chopper_command_dict, minibar_command_dict):
    """
    Performs minibar preparation for NGSpeciesID processing.
    
    This function handles the preparation process for running minibar, including file extraction, 
    command execution, and file cleanup. It logs various steps and outcomes of the process.
    
    Parameters:
        minibar_path (str): Path to the minibar script.
        minibar_index_path (str): Path to the minibar index file.
        input_file_path (str): Path to the compressed FASTQ file.
        ngsid_output_dir (str): Directory path for storing output files.
        chopper_command_dict (dict): Dictionary containing chopper settings.
        minibar_command_dict (dict): Dictionary containing MiniBar settings.
    
    Global Variables:
        output_area (str): A global variable used for logging output messages.
    
    Notes:
        - Extracts .fastq.gz files if they don't already exist.
        - Constructs and executes the minibar command.
        - Handles and logs file removal errors.
    """
    global DEFAULT_LOG_FILE
    main_dir = os.getcwd()
    log_print(f"Prepping {input_file_path} for NGSpeciesID with minibar...")

    # Determine if the file is compressed (.gz) and set the correct chopper command
    if input_file_path.endswith(".gz"):
        log_print(f"File is compressed: {input_file_path}")
        chopper_input_cmd = f"gunzip -c {input_file_path}"
        uncompressed_fastq_path = input_file_path.replace(".gz","")  # Remove .gz extension
    else:
        log_print(f"File is uncompressed: {input_file_path}")
        chopper_input_cmd = f"cat {input_file_path}"
        uncompressed_fastq_path = input_file_path

    # Ensure the output directory exists
    if not os.path.exists(ngsid_output_dir):
        os.makedirs(ngsid_output_dir)
    os.chdir(ngsid_output_dir)

    # Construct the chopper command
    chopper_cmd_str = f"{chopper_input_cmd} | chopper -q {chopper_command_dict['-q']} -l {chopper_command_dict['--minlength']} --maxlength {chopper_command_dict['--maxlength']} --threads {CPU_THREADS}"
    
    # Construct minibar command
    minibar_cmd_str = f"{minibar_path} -F {minibar_index_path} -i - --outfolder {ngsid_output_dir}"
    for key, value in minibar_command_dict.items():
        if value != "":
            minibar_cmd_str += f" {key} {value}"
        elif value == True:  # Assuming some keys might be flags without explicit values
            minibar_cmd_str += f" {key}"

    # Combine the chopper and minibar commands
    combined_cmd_str = f"{chopper_cmd_str} | {minibar_cmd_str}"

    run_subprocess_cmd(combined_cmd_str, shell_check=True)
    
    # TODO: Do not remove these files until future notice
    # for file in [f"{ngsid_output_dir}/sample_Multiple_Matches.fastq", f"{ngsid_output_dir}/sample_unk.fastq", uncompressed_fastq_path]:
    #     try:
    #         os.remove(file)
    #         log_print( f"PASS:\tRemoved file: {file}")
    #     except OSError as e:
    #         log_print( f"ERROR:\tDuring file removal {file}: {e}")

    os.chdir(main_dir)

if  __name__ == "__main__":
    log_print("UNLOGGED DEBUG:\tUNDEVELOPED FunDiS_Minibar.py DEBUG AREA")
