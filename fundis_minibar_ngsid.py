# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 2023

@author: ian.michael.bollinger@gmail.com with the support of ChatGPT 4.0

This is a module intended to be used as a part of a pipeline.

This can be used individually by calling the command:
    python /path/to/fundis_minibar_ngsid.py -i /path/to/input.fastq -x /path/to/index.txt -t /path/to/primers.txt -p 80
    python /path/to/fundis_minibar_ngsid.py --input_fastq /path/to/input.fastq ---minbar_index_path /path/to/index.txt --primers_text_path /path/to/primers.txt --percent_system_use 80
"""

import psutil
import os
import platform
import sys
import subprocess
import multiprocessing
import math
import glob
import concurrent.futures
from tqdm import tqdm
import argparse

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

# Function to index a given fastq file with minibar.py provided by FunDiS:
def index_fastq_with_minibar(input_fastq, path_to_minibar, path_to_minibar_index):
    print(f'\nIndexing {input_fastq} with minibar...')
    minbar_dir = input_fastq.replace('.fastq','_minibar')
    if os.path.isdir(minbar_dir):
        print(f"PASS: Skipping minibar, {minbar_dir} already exists")
        return minbar_dir
    else:
        os.makedirs(minbar_dir, exist_ok=True)  # add exist_ok=True
    os.chdir(minbar_dir)
    # Run command
    cmd = f'{path_to_minibar} -F {path_to_minibar_index} {input_fastq}'
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, error = process.communicate()
    print(output.decode())
    if process.returncode != 0:
        print(f"ERROR: There was a problem running minibar. Return code: {process.returncode}")
        return None
    else:
        print("PASS: Successfully indexed fastq with minibar")
        return minbar_dir

# Function to generate consensus 
def consensus_align_with_ngspeciesid(input_fastq, primers_path):
    print(f'\nGenerating raw alignment and consensus reads for {input_fastq} with NGSpeciesID...')     
    output_folder = input_fastq.replace('.fastq','_NGSequenceID')
    # Check if the output file already exists
    if os.path.isdir(output_folder):
        print(f"PASS: Skipping NGSpeciesID, {output_folder} already exists")
        return output_folder
    else:
        os.makedirs(output_folder, exist_ok=True)  # add exist_ok=True
    # If the file doesn't exist, run Pilon
    cmd = f"NGSpeciesID --ont --consensus --sample_size 500 --m 730 --s 400 --medaka --primer_file {primers_path} --fastq {input_fastq} --outfolder {output_folder}"
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, error = process.communicate()
    print(output.decode())
    if process.returncode != 0:
        print(f"ERROR: There was a problem running NGSpeciesID. Return code: {process.returncode}")
    else:
        print("PASS: Successfully genreated NGSpeciesID output folder...")

def main(args):   
    # Parse user arguments
    input_fastq = args.input_fastq if args.input_fastq else f"{environment_dir}/Fundis/TEST/combined2.fastq"
    path_to_minibar_index = args.minbar_index_path if args.minbar_index_path else "/mnt/e/Fundis/Programs-20230719T043926Z-001/Programs/Index.txt"
    primers_path = args.primers_path if args.primers_path else f"{environment_dir}/Fundis/Programs-20230719T043926Z-001/Programs/primers.txt"
    percent_system_use = float(args.percent_system_use/100) if args.percent_system_use else 0.8
    
    # Get the number of CPUs available on the system
    num_cpus = multiprocessing.cpu_count()
    
    # Get the amount of RAM (GB) currently available
    mem_info = psutil.virtual_memory()
    
    # Calculate the number of threads as 80% of available CPUs & RAM
    cpu_threads = int(math.floor(num_cpus * ))
    ram_gb = int(mem_info.total / (1024.0 ** 3) * percent_system_use)    
        
    # TODO: Update minibar.py to be its own module
    path_to_minibar = f"{environment_dir}/Fundis/Programs-20230719T043926Z-001/Programs/minibar.py"
    minbar_dir = index_fastq_with_minibar(input_fastq, path_to_minibar, path_to_minibar_index)
    os.chdir(main_working_dir)
   
    # Get a list of all .fastq files in the directory
    fastq_files = glob.glob(os.path.join(minbar_dir, "*.fastq"))
    
    # Use a ThreadPoolExecutor to run the NGSpeciesID command in parallel for each .fastq file
    with concurrent.futures.ThreadPoolExecutor(max_workers=cpu_threads) as executor:
        futures = {executor.submit(consensus_align_with_ngspeciesid, fastq_file, primers_path) for fastq_file in fastq_files}
        for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc="Processing", unit="task"):
            output_folder = future.result()
            if output_folder is not None:
                print(f"Finished processing: {output_folder}")

    print("PASS: FunDiS Minibar.py & NGSequenceID complete\n")
                
# If this script is the main entry point, parse the arguments and call the main function
if __name__ == "__main__":
    # Global environment_dir
    environment_dir = ""
    environment_cmd_prefix = ""
    environment_dir = check_os()
    main_working_dir = os.getcwd()
    
    # Parse user arguments
    parser = argparse.ArgumentParser(description="Process ONT FASTQ file with minibar.py.")
    parser.add_argument('-i', '--input_fastq', type=str, help='Path to the FASTQ file containing ONT nrITS data')
    parser.add_argument('-x', '--minbar_index_path', type=str, help='Path to Text file containing the minibar index to parse input_fastq.')
    parser.add_argument('-t', '--primers_text_path', type=str, help='Path to Text file containing the Primers used to generate input_fastq.')
    parser.add_argument('-p', '--percent_system_use', type=str, help='Percent system use written as integer.')
    args = parser.parse_args()
    main(args)