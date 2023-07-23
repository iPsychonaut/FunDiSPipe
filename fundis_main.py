# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 2023

@author: ian.michael.bollinger@gmail.com with the support of ChatGPT 4.0

This is a module intended to be used as a part of a pipeline.

This can be used individually by calling the command:
    python /path/to/fundis_main.py -i /path/to/input.fastq -x /path/to/index.txt -t /path/to/primers.txt -p 50
    python /path/to/fundis_main.py --input /path/to/input.fastq ---minbar_index_path /path/to/index.txt --primers_text_path /path/to/primers.txt --percent_system_use 50
"""

import os
import platform
import sys
import subprocess
import pkg_resources
import multiprocessing
import math
import glob
import concurrent.futures
import argparse
import queue
import shutil
from fundis_minibar_ngsid import minibar_ngsid
from fundis_haplotype_phaser import haplotype_phaser
from fundis_summarize import summarize

# Function to check that all Python Libraries are installed, and if not installs them
def libraries_check(libraries):
    print("Checking Python Library Prerequisites...")
    installed = {pkg.key for pkg in pkg_resources.working_set}
    missing = [lib for lib in libraries if lib not in installed]
    if missing:
        print(f"NOTE: Missing libraries: {missing}. Installing...")
        python = sys.executable
        for lib in missing:
            process = subprocess.Popen([python, '-m', 'pip', 'install', lib], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            output, error = process.communicate()
            print(output.decode())
            if process.returncode != 0:
                print(f"ERROR: There was a problem installing {lib}. Return code: {process.returncode}")
    else:
        print("PASS: All libraries are installed\n")

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
    print(f'Environment Directory: {environment_dir}\n')
    return environment_dir

# Function to check if File is installed in the current working environment
def check_file_installed(input_file, environment_dir):
    # check the os currently working in
    # Run check if file exists
    cmd = f'[ -f $(find . -name "{input_file}" -print -quit) ] && echo "{input_file} is installed" || echo "{input_file} is not installed"'
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, error = process.communicate()
    print(output.decode())
    if process.returncode != 0:
        print(f"ERROR: There was a problem checking the installation of {input_file}. Return code: {process.returncode}")
    return(return_code)

def check_tools_installed(tools, environment_dir):
    print("Checking Pipeline Prerequisites...")   
    
    missing_tools = []

    if platform.system() == 'Windows':
        cmd_tools = ' && '.join([f'{tool} --version 2>NUL || echo {tool} is not installed' for tool in tools])
    else:
        cmd_tools = ' && '.join([f'command -v {tool} >/dev/null 2>&1 || echo {tool} is not installed' for tool in tools])

    process = subprocess.Popen(cmd_tools, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, error = process.communicate()
    output_lines = output.decode().splitlines()
    if process.returncode != 0:
        print(f"ERROR: There was a problem checking the installation of tools. Return code: {process.returncode}")
        for line in output_lines:
            if "is not installed" in line:
                missing_tool = line.split(" ")[0]
                missing_tools.append(missing_tool)
        print("Missing tools: ", missing_tools)
    else:
        print("PASS: All tools checked successfully\n")
    
    for tool in missing_tools:
        if tool == 'medaka':
            tool = "medaka==0.11.5"
        elif tool == 'openblas':
            tool = 'openblas==0.3.3'
        print(f"Trying to install {tool} using conda...")
        cmd = f"{environment_cmd_prefix}conda install -y -c bioconda {tool}"
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output, error = process.communicate()
        if process.returncode != 0:
            print(f"ERROR: There was a problem installing {tool}. Return code: {process.returncode}")
            print(output.decode())
        else:
            print(f"PASS: {tool} installed successfully\n")

# Global environment_dir
environment_dir = ""
environment_cmd_prefix = ""
environment_dir = check_os()
main_working_dir = os.getcwd()

# Check Pipeline Prerequisites
libraries = ["tqdm", "pandas", "pysam", "biopython"]
libraries_check(libraries)
tools = ['NGSpeciesID', 'bwa', 'samtools', 'bcftools', 'whatshap', 'medaka', 'openblas', 'spoa']
check_tools_installed(tools, environment_dir)

from tqdm import tqdm
import pandas as pd
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main(args):
    # Parse user arguments
    input_fastq = args.input if args.input else f"{environment_dir}/Fundis/TEST/combined2.fastq"
    path_to_minibar_index = args.minbar_index_path if args.minbar_index_path else "/mnt/e/Fundis/Programs-20230719T043926Z-001/Programs/Index.txt"
    primers_text_path = args.primers_text_path if args.primers_text_path else f"{environment_dir}/Fundis/Programs-20230719T043926Z-001/Programs/primers.txt"
    percent_system_use = float(args.percent_system_use)/100 if args.percent_system_use else 0.5
    minbar_dir = input_fastq.replace('.fastq','_minibar')

    # Get the number of CPUs available on the system
    num_cpus = multiprocessing.cpu_count()
   
    # Calculate the number of threads as 80% of available CPUs & RAM
    cpu_threads = int(math.floor(num_cpus * percent_system_use))
    
    # Try to run the MiniBar & NGSpeciesID functions
    try:
        minibar_ngsid(args)
    except Exception as e:
        print(f"ERROR: There was a problem running minibar_ngsid: {str(e)}")

    # Try to run the Haplotype Phaser function
    try:
        haplotype_phaser(args)
    except Exception as e:
        print(f"ERROR: There was a problem: {str(e)}")
    
    # Try to run the Summarize function
    try:
        summarize(args)
    except Exception as e:
        print(f"ERROR: There was a problem: {str(e)}")
    
# If this script is the main entry point, parse the arguments and call the main function
if __name__ == "__main__":
    # Parse user arguments
    parser = argparse.ArgumentParser(description="Process ONT FASTQ file with minibar.py.")
    parser.add_argument('-i', '--input', type=str, help='Path to the FASTQ file containing ONT nrITS data')
    parser.add_argument('-x', '--minbar_index_path', type=str, help='Path to Text file containing the minibar index to parse input_fastq.')
    parser.add_argument('-t', '--primers_text_path', type=str, help='Path to Text file containing the Primers used to generate input_fastq.')
    parser.add_argument('-p', '--percent_system_use', type=str, help='Percent system use written as integer.')
    args = parser.parse_args()
    main(args)
