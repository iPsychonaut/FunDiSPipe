# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 2023

@author: ian.michael.bollinger@gmail.com with the support of ChatGPT 4.0

This is a module intended to be used as a part of a pipeline.

This can be used individually by calling the command:
    python /path/to/fundis_main.py -i /path/to/input.fastq -t /path/to/primers.txt -p 80
    python /path/to/fundis_main.py --input_fastq /path/to/input.fastq --primers_text_path /path/to/primers.txt --percent_system_use 80
"""

import psutil
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
# from tqdm import tqdm

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
            
def main(args):
    # Parse user arguments
    input_fastq = args.input_fastq if args.input_fastq else f"{environment_dir}/Fundis/TEST/combined2.fastq"
    primers_text_path = args.primers_text_path if args.primers_text_path else ''
    minbar_dir = input_fastq.replace('.fastq','_minibar')
    percent_system_use = float(args.percent_system_use/100) if args.percent_system_use else 0.8

    # Check Pipeline Prerequisites
    libraries = ["psutil", "tqdm", "pandas", "pysam", "biopython"]
    libraries_check(libraries)
    tools = ['NGSpeciesID', 'bwa', 'samtools', 'bcftools', 'whatshap', 'medaka', 'openblas', 'spoa']
    check_tools_installed(tools, environment_dir)

    # Get the number of CPUs available on the system
    num_cpus = multiprocessing.cpu_count()
    
    # Get the amount of RAM (GB) currently available
    mem_info = psutil.virtual_memory()
    
    # Calculate the number of threads as 80% of available CPUs & RAM
    cpu_threads = int(math.floor(num_cpus * percent_system_use))
    ram_gb = int(mem_info.total / (1024.0 ** 3) * percent_system_use)    
    
    # TODO: python fundis_minibar_ngsid.py args.input_fastq, args.primers_text_path and args.percent_system_use
    print("TODO: python fundis_minibar_ngsid.py args.input_fastq, args.primers_text_path and args.percent_system_use")
    minibar_ngsid_cmd = f'python fundis_minibar_ngsid.py -i {input_fastq} -t {args.primers_text_path} -p {args.percent_system_use}'
    minibar_ngsid_process = subprocess.Popen(minibar_ngsid_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    minibar_ngsid_output, minibar_ngsid_error = minibar_ngsid_process.communicate()
    print(minibar_ngsid_output.decode())
    if minibar_ngsid_process.returncode != 0:
        print(f"ERROR: There was a problem: {minibar_ngsid_process.returncode}")
    else:
        print("PASS: Successfully parsed input FASTQ with minibar & NGSequenceID")
    
    # TODO: ptyhon fundis_haplotype_phaser.py PATH FOR HAPLOTYPE PARSER and args.percent_system_use
    print("TODO: python fundis_haplotype_phaser.py minbar_dir and args.percent_system_use")    
    hyplotype_parser_cmd = f'python fundis_haplotype_phaser.py -i {minbar_dir} -p {args.percent_system_use}'
    hyplotype_parser_process = subprocess.Popen(hyplotype_parser_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    hyplotype_parser_output, hyplotype_parser_error = hyplotype_parser_process.communicate()
    print(hyplotype_parser_output.decode())
    if hyplotype_parser_process.returncode != 0:
        print(f"ERROR: There was a problem: {hyplotype_parser_process.returncode}")
    else:
        print("PASS: Successfully phased hyplotypes for NGSequenceID folder")
    
    # TODO: python fundis_summarize.py PATH FOR HAPLOTYPE PARSER and args.percent_system_use
    print("TODO: python fundis_summarize.py minbar_dir and args.percent_system_use")
    summarize_cmd = f'python fundis_summarize2.py -s {minbar_dir} -p {args.percent_system_use}'
    summarize_process = subprocess.Popen(summarize_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    summarize_output, summarize_error = summarize_process.communicate()
    print(summarize_output.decode())
    if summarize_process.returncode != 0:
        print(f"ERROR: There was a problem: {summarize_process.returncode}")
    else:
        print("PASS: Successfully summarized NGSequenceID folder for MycoMap")
    
    print("\n*****PASS: FunDiS Pipeline complete*****")
    
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
    parser.add_argument('-t', '--primers_text_path', type=str, help='Path to Text file containing the Primers used to generate input_fastq.')
    parser.add_argument('-p','--percent_system_use', type=str, help='Percent system use written as integer.')
    args = parser.parse_args()
    main(args)