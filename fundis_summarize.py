# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 2023

@author: ian.michael.bollinger@gmail.com with the support of ChatGPT 4.0

This is a module intended to be used as a part of a pipeline.

This can be used individually by calling the command:
    python /path/to/fundis_summarize.py -i /path/to/input_dir -p 50
    python /path/to/fundis_summarize.py --input /path/to/input_dir --percent_system_use 50
"""
from Bio import SeqIO
import pandas as pd
import os
import glob
import shutil
from tqdm import tqdm
import argparse
import platform

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

# Function to take in a folder containined processed NGSequenceID folders and generate a summary folder for MycoMap
def mycomap_summarize_ngsid_dir(ngsid_dir):
    if '.fastq' in ngsid_dir:
        ngsid_dir = ngsid_dir.replace('.fastq','_minibar_NGSID')
    print(f'\nGenerating MycoMap summary folder for {ngsid_dir}...')
    # Create the summary and FASTQ directories if they don't exist
    summary_dir = ngsid_dir.replace("_minibar_NGSID","_Summary")
    
    # Create pandas DataFrame for sequence information
    stats_df = pd.DataFrame(columns=['Filename', 'Length', 'Reads in Consensus', 'Multiple'])
    
    # Get all directories in the ngsid_dir
    sample_dirs = [d for d in os.listdir(ngsid_dir) if os.path.isdir(os.path.join(ngsid_dir, d))]    
    
    # Generate storage folders
    fastq_dir = f"{summary_dir}/FASTQ Files"
    os.makedirs(summary_dir, exist_ok=True)
    os.makedirs(fastq_dir, exist_ok=True)
    
    # Look through current_sample_dir in sample_dirs
    for current_sample_dir in tqdm(sample_dirs):
        current_sample_dir = f'{ngsid_dir}/{current_sample_dir}'
        base_name = current_sample_dir.split("sample_")[-1]
        
        # Establish main variables and directories for processing
        consensus_dirs = [entry.path for entry in os.scandir(current_sample_dir) if entry.is_dir() and entry.name.startswith('consensus_reference_') and any(os.scandir(entry.path))]
        consensus_dirs = [entry.replace("\\","/") for entry in consensus_dirs]
        consensus_fastq_list = []
        medaka_count = 1
        reads_in_consensus = 0
        
        for consensus_dir in consensus_dirs:
            # Use glob to get a list of all files that match the pattern "_individual.fasta"
            individual_fasta_files = glob.glob(f"{consensus_dir}/**/*_individual.fasta", recursive=True)
            haplotype_count = len(individual_fasta_files)
            
            # Generate filenames for fasta and fastq associated with the current consensus
            fasta_file = f'{consensus_dir}.fasta'
            fastq_file = f'{consensus_dir}.fastq'
            fastq_file = fastq_file.replace('consensus_reference','reads_to_consensus')
            
            # Get number of reads in consensus
            for record in SeqIO.parse(fasta_file, "fasta"):
                sequence_length = len(record.seq)
                record_id = record.id
                reads_in_consensus = record_id.split('total_supporting_reads_')[-1]
            new_fastq_filename = base_name.replace('_NGSequenceID', f'-{medaka_count}-RiC{reads_in_consensus}.fastq')
            new_fastq_filepath = os.path.join(fastq_dir, new_fastq_filename)        
            shutil.move(fastq_file, new_fastq_filepath)
            if haplotype_count == 1:
                # Rename the consensus_reference_#.fasta in current_sample_dir to use the base_name but replace '_NGSequenceID' with '-{medaka_count}-RiC{reads_in_consensus}.fasta'
                haplotype_letter = chr(65)
                new_fasta_filename = base_name.replace('_NGSequenceID', f'-{medaka_count}-RiC{reads_in_consensus}-{haplotype_letter}.fasta')
                new_fasta_filepath = os.path.join(summary_dir, new_fasta_filename)
                
                # Relocate the '-{medaka_count}-RiC{reads_in_consensus}.fasta' file to the summary_dir
                shutil.move(fasta_file, new_fasta_filepath)
                
                # Rename the new_fastq_filepath file to include the haplotype_letter
                new_fastq_filepath_1 = new_fasta_filepath.replace('.fasta','.fastq')
                new_fastq_filepath_1 = new_fastq_filepath_1.replace(f'{summary_dir}',f'{fastq_dir}')
                os.rename(new_fastq_filepath, new_fastq_filepath_1)
                new_base_name = f'{base_name}-{medaka_count}-{haplotype_letter}'
                stats = [new_base_name, sequence_length, reads_in_consensus, medaka_count]
                stats_df.loc[len(stats_df)] = stats
                for i, haplotype_fasta_file in enumerate(individual_fasta_files):
                    # Rename the individual fasta files to use the base_name but replace '_NGSequenceID' with '-{medaka_count}-RiC{reads_in_consensus}-A.fasta' or '-{medaka_count}-RiC{reads_in_consensus}-B.fasta'
                    haplotype_letter = chr(66+i)
                    new_fasta_filename = base_name.replace('_NGSequenceID', f'-{medaka_count}-RiC{reads_in_consensus}-{haplotype_letter}.fasta')
                    new_fasta_filepath = os.path.join(summary_dir, new_fasta_filename)
                    
                    # Relocate the '-{medaka_count}-RiC{reads_in_consensus}-A.fasta' & '-{medaka_count}-RiC{reads_in_consensus}-B.fasta' files to the summary_dir
                    shutil.move(haplotype_fasta_file, new_fasta_filepath)
                    
                    # Duplicate the new_fastq_filepath file to include the the updated haplotype_letter
                    new_fastq_filepath_i = new_fasta_filepath.replace('.fasta','.fastq')
                    new_fastq_filepath_i = new_fastq_filepath_i.replace(f'{summary_dir}',f'{fastq_dir}')
                    
                    # Duplicate and rename the file
                    shutil.copy(new_fastq_filepath_1, new_fastq_filepath_i)
                    new_base_name = f'{base_name}-{medaka_count}-{haplotype_letter}'
                    stats = [new_base_name, sequence_length, reads_in_consensus, medaka_count]
                    stats_df.loc[len(stats_df)] = stats
            elif haplotype_count == 2:
                for i, haplotype_fasta_file in enumerate(individual_fasta_files):
                    # Rename the individual fasta files to use the base_name but replace '_NGSequenceID' with '-{medaka_count}-RiC{reads_in_consensus}-A.fasta' or '-{medaka_count}-RiC{reads_in_consensus}-B.fasta'
                    haplotype_letter = chr(65+i)
                    new_fasta_filename = base_name.replace('_NGSequenceID', f'-{medaka_count}-RiC{reads_in_consensus}-{haplotype_letter}.fasta')
                    new_fasta_filepath = os.path.join(summary_dir, new_fasta_filename)
                    
                    # Relocate the '-{medaka_count}-RiC{reads_in_consensus}-A.fasta' & '-{medaka_count}-RiC{reads_in_consensus}-B.fasta' files to the summary_dir
                    shutil.move(haplotype_fasta_file, new_fasta_filepath)
                    
                    # Duplicate the new_fastq_filepath file to include the the updated haplotype_letter
                    new_fastq_filepath_i = new_fasta_filepath.replace('.fasta','.fastq')
                    new_fastq_filepath_i = new_fastq_filepath_i.replace(f'{summary_dir}',f'{fastq_dir}')
                    shutil.copy(new_fastq_filepath, new_fastq_filepath_i)
                    new_base_name = f'{base_name}-{medaka_count}-{haplotype_letter}'
                    stats = [new_base_name, sequence_length, reads_in_consensus, medaka_count]
                    stats_df.loc[len(stats_df)] = stats
                os.remove(new_fastq_filepath)
            else:
                # Rename the consensus_reference_#.fasta & fastq in current_sample_dir to use the base_name but replace '_NGSequenceID' with '-{medaka_count}-RiC{reads_in_consensus}.fasta'
                new_fasta_filename = base_name.replace('_NGSequenceID', f'-{medaka_count}-RiC{reads_in_consensus}.fasta')
                new_fasta_filepath = os.path.join(summary_dir, new_fasta_filename)
                
                # Relocate the -{medaka_count}-RiC{reads_in_consensus}.fasta & .fastq file to the summary_dir
                shutil.move(fasta_file, new_fasta_filepath)
                new_base_name = f'{base_name}-{medaka_count}'
                stats = [new_base_name, sequence_length, reads_in_consensus, medaka_count]
                stats_df.loc[len(stats_df)] = stats
            medaka_count += 1
        
        # Write the tab separated Summary.txt file based on the columns in the stats dataframe
        summary_name = summary_dir.split("/")[-1]
        stats_df.to_csv(os.path.join(summary_dir, f'{summary_name}.txt'), sep='\t', index=False)
        
        # Get a list of all fasta files in the directory
        fasta_files = glob.glob(f"{summary_dir}/*.fasta")
        
        # Create the new combined fasta file
        combined_fasta_filename = os.path.join(summary_dir, f"{summary_name}.fasta")
        with open(combined_fasta_filename, 'w') as combined_fasta_file:
            for fasta_file in fasta_files:
                # Get the base name of the file (strip off the .fasta extension)
                new_id = os.path.basename(fasta_file).replace('.fasta', '')
        
                # Parse the fasta file and update the id of each record
                records = list(SeqIO.parse(fasta_file, 'fasta'))
                for record in records:
                    record.id = new_id
                    record.description = ''
        
                # Write the updated records to the combined fasta file
                SeqIO.write(records, combined_fasta_file, 'fasta')
                
def summarize(args):
    try:
        # Global environment_dir
        environment_dir = ""
        environment_cmd_prefix = ""
        environment_dir = check_os()
        main_working_dir = os.getcwd()
        
        # Set the path to the folder containing the fastq files
        percent_system_use = float(args.percent_system_use)/100 if args.percent_system_use else 0.5
        input_dir = args.input if args.input else os.path.dirname(os.path.realpath(__file__))
        mycomap_summarize_ngsid_dir(input_dir)
        print('PASS: Successfully summarized NGSequenceID Folder for MycoMap upload')
        return True
    
    except Exception as e:
        print(f"ERROR: There was a problem in summarize: {str(e)}")
        return False

if __name__ == "__main__":
    # Parse user arguments
    parser = argparse.ArgumentParser(description="Process NGSpeciesID source folder.")
    parser.add_argument('-i','--input', type=str, help='Path to the NGSpeciesID source folder')
    parser.add_argument('-p','--percent_system_use', type=str, help='Percent system use written as integer.')
    args = parser.parse_args()    
    summarize(args)
