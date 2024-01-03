# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 12:23:15 2023

@author: ian.michael.bollinger@gmail.com/researchconsultants@critical.consulting

This is the Mini-Barcoder (minibar.py) main wrapper developed for a modified

protocol by Stephen Douglas Russell and paid for by the Fungal Diversity Survey (FunDiS).

Protocol Link: https://www.protocols.io/view/primary-data-analysis-basecalling-demultiplexing-a-dm6gpbm88lzp/v3?step=3
"""

# Base Python Imports
import subprocess, os, gzip

# Required Python Imports
from termcolor import cprint
from datetime import datetime

# Set Global File
DEFAULT_LOG_FILE = None

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
def run_minibar_prep(minibar_path, minibar_index_path, ngsid_fastq_gz_path, ngsid_output_dir):
    """
    Wrapper function to run minibar preparation in a separate thread.
    
    This function is designed to execute the minibar preparation process in an isolated thread. 
    It handles exceptions and logs any errors encountered during the execution.
    
    Parameters:
        minibar_path (str): Path to the minibar script.
        minibar_index_path (str): Path to the minibar index file.
        ngsid_fastq_gz_path (str): Path to the compressed FASTQ file.
        ngsid_output_dir (str): Directory path for storing output files.
    
    Notes:
        - Logs an error message if an exception is encountered during the minibar preparation.
    """
    global DEFAULT_LOG_FILE
    try:
        minibar_prep(minibar_path, minibar_index_path, ngsid_fastq_gz_path, ngsid_output_dir)
    except Exception as e:
        log_print( f"ERROR:\t{str(e)}")

# Function to perform minibar preparation  
def minibar_prep(minibar_path, minibar_index_path, ngsid_fastq_gz_path, ngsid_output_dir):
    """
    Performs minibar preparation for NGSpeciesID processing.
    
    This function handles the preparation process for running minibar, including file extraction, 
    command execution, and file cleanup. It logs various steps and outcomes of the process.
    
    Parameters:
        minibar_path (str): Path to the minibar script.
        minibar_index_path (str): Path to the minibar index file.
        ngsid_fastq_gz_path (str): Path to the compressed FASTQ file.
        ngsid_output_dir (str): Directory path for storing output files.
    
    Global Variables:
        output_area (str): A global variable used for logging output messages.
    
    Notes:
        - Extracts .fastq.gz files if they don't already exist.
        - Constructs and executes the minibar command.
        - Handles and logs file removal errors.
    """
    global DEFAULT_LOG_FILE
    log_print( f"Prepping {ngsid_fastq_gz_path} for NGSpeciesID with minibar...")
    ngsid_fastq_path = ngsid_fastq_gz_path.replace(".gz", "")
    
    if os.path.exists(ngsid_fastq_path):
        log_print( "PASS:\tSkipping extraction, files already exist")
    else:
        unzip_fastq(ngsid_fastq_gz_path, ngsid_fastq_path)

    main_dir = os.getcwd()
    if not os.path.exists(ngsid_output_dir):
        os.makedirs(ngsid_output_dir)
    os.chdir(ngsid_output_dir)
        
    minibar_cmd = [minibar_path, "-F", minibar_index_path, ngsid_fastq_path] # TODO: Be able to adjust the number of differences allowed for in primer sequence
    log_print( f"CMD:\t{' '.join(minibar_cmd)}")
    process = subprocess.run(minibar_cmd, capture_output=True, text=True)
    log_print( process.stdout)
    log_print( process.stderr)
    
    if process.returncode != 0:
        log_print( f"ERROR:\t{process.stderr}")
    else:
        log_print( f"PASS:\tSuccessfully processed {ngsid_fastq_path} with minibar and cleaned up files")
    
    for file in [f"{ngsid_output_dir}/sample_Multiple_Matches.fastq", "sample_unk.fastq"]:
        try:
            os.remove(file)
            log_print( f"PASS:\tRemoved file: {file}")
        except OSError as e:
            log_print( f"ERROR:\tDuring file removal {file}: {e}")

    os.chdir(main_dir)

if  __name__ == "__main__":
    global DEFAULT_LOG_FILE
    log_print("UNLOGGED DEBUG:\tUNDEVELOPED FunDiS_Minibar.py DEBUG AREA")
