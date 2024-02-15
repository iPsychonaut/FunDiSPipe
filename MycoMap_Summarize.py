# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 12:23:15 2023

@author: ian.michael.bollinger@gmail.com/researchconsultants@critical.consulting

This is the MycoMap Summarizer (summarize.py) developed for a modified protocol by

Stephen Douglas Russell and paid for by the Fungal Diversity Survey (FunDiS).

Protocol Link: https://www.protocols.io/view/primary-data-analysis-basecalling-demultiplexing-a-dm6gpbm88lzp/v3?step=3
"""

# Base Python Imports
import csv, os, re, glob, shutil, time
from pathlib import Path
from Bio import SeqIO

# Custom Python Imports
from FunDiS_Tools import log_print, generate_log_file, initialize_logging_environment, run_subprocess_cmd

# Global output_area variable
PERCENT_RESOURCES = 0.75

def process_sample_folder(sample_folder, summary_folder, summary_ric_cutoff, hap_phase_bool):
    """
    Processes a sample folder by iterating through its 'medaka' subfolders and, if applicable, 
    processes the phased FASTA file. Collects and appends statistics from each processed folder.

    Args:
        sample_folder (str): The path of the sample folder to be processed.
        summary_folder (str): The path to the summary folder.
        summary_ric_cutoff (str): Minimum number of reads needed to be processed.
        hap_phase_bool (bool): True/False on whether Haplotype Phasing was performed.

    Returns:
        folder_stats (list): NEEDS DEESCRIPTION.
        
    Notes:
        - Assumes the presence of 'medaka' subfolders for processing.
        - If a phased FASTA file is present and valid, its data will also be processed.
        - Skips processing of 'medaka' folders that match the phased ID, if applicable.

    Example:
        process_sample_folder("/path/to/sample_folder")
    """
    # set phased_id = None
    phased_id = None
    
    # Generate folder_stats list
    folder_stats = []  
    
    # Collect all 'racon' subfolders in the given directory
    racon_folders = glob.glob(os.path.join(sample_folder, "racon*"))

    # Collect all 'medaka' subfolders in the given directory
    medaka_folders = glob.glob(os.path.join(sample_folder, "medaka*"))

    # Check if hap_phase_bool == True
    if hap_phase_bool == True:
        # if True: Process Phased Haplotype FASTA
        for index, racon_folder in enumerate(racon_folders):
            try:
                # Phased Haplotype FASTA Function
                phased_id, filename, length, reads = process_phased_fasta(racon_folder, summary_folder, summary_ric_cutoff)
                
                if phased_id != None:
                    # Get and Append Stats Function
                    folder_stats.append([filename, length, reads, phased_id])
                    
                    log_print(f"PASS:\tProcessed Phased Haplotype FASTA: {filename}, Length: {length}, Reads: {reads}, ID: {phased_id}")
            except:
                medaka_id, filename, length, reads = process_medaka_folder(medaka_folders[index], summary_folder, summary_ric_cutoff)
                
                if filename != None:
                    # Get and Append Stats Function
                    folder_stats.append([filename, length, reads, medaka_id])
                    
                    log_print(f"PASS:\tProcessed Medaka Folder: {filename}, Length: {length}, Reads: {reads}, ID: {medaka_id}")
                
    else:        
        # Iterate over each medaka_folder in the list
        for medaka_folder in medaka_folders:
            medaka_id, filename, length, reads = process_medaka_folder(medaka_folder, summary_folder, summary_ric_cutoff)
            
            if filename != None:
                # Get and Append Stats Function
                folder_stats.append([filename, length, reads, medaka_id])
                
                log_print(f"PASS:\tProcessed Medaka Folder: {filename}, Length: {length}, Reads: {reads}, ID: {medaka_id}")
    
    return folder_stats

def process_phased_fasta(sample_folder, summary_folder, summary_ric_cutoff):
    """
    Finds and processes a phased FASTA file within the given directory. Extracts relevant 
    information like the phased ID, total supporting reads, and sequence length.

    Args:
        sample_folder (str): The directory where the phased FASTA file is located.
        summary_folder (str): The path to the summary folder.
            
    Returns:
        tuple: A tuple containing phased_id, filename, length, and reads information.

    Notes:
        - The function looks for a file ending with '_phased.fasta'.
        - Extracts the phased ID from the filename and supporting reads and length from the sequence name.

    Example:
        filename, length, reads = process_phased_fasta("/path/to/sample_folder")
    """
    log_print(f"Attempting to process Phased FASTA file in {sample_folder}...")
    
    # Find phased fasta file in given directory with "_phased.fasta" in its name
    temp_path = sample_folder.replace("racon_cl_id","reads_to_consensus")
    phased_fasta_path = f'{temp_path}_phased.fasta'
    # phased_fasta_paths = glob.glob(os.path.join(sample_folder, "*_phased.fasta"))
    # phased_fasta_path = phased_fasta_paths[0]

    if not phased_fasta_path:
        log_print(f"ERROR:\tNo Phased FASTA file exists: {phased_fasta_path}")
        return None, None, None, None
    
    reads = None
    length = None
    filename = None
    parent_folder_name = None
    
    # Split the phased_fasta_path string by "_" and save the last item in the list as "phased_id" variable
    phased_id = phased_fasta_path.replace('_phased.fasta', '').split("_")[-1]

    # Use Biopython to open the fasta file
    with open(phased_fasta_path, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            # Check to see if the sequence name is as follows: "consensus_cl_id_{medaka_id}_total_supporting_reads_[0-9*] Haplotype Phased total supporting reads [0-9*]"
            if "Haplotype Phased total supporting reads" in record.description:
                # Save the last item in the seuqence name split by " " list as the "reads" variable
                reads = record.description.split()[-1]
                print(reads)
                if int(reads) >= int(summary_ric_cutoff):
                    # Get the length of the sequence and save it as "length" variable        
                    length = len(record.seq)
    
                    # Generate the variable parent_folder_name = Path(phased_fasta_path).parent.name.replace("sample_","")
                    parent_folder_name = Path(phased_fasta_path).parent.name.replace("sample_", "")
    
                    # Generate the variable filename = f"{parent_folder_name}-{phased_id}"
                    filename = f"{parent_folder_name}-phased-{phased_id}"
    
                    log_print(f"PASS:\tCorrect tags found in {phased_fasta_path}: reads: {reads}; length: {length}")
                    
                    updated_fasta_path = os.path.join(summary_folder, f'{filename}-RiC{reads}.fasta')
                    
                    move_updated_fasta(updated_fasta_path, summary_folder, phased_fasta_path, parent_folder_name, filename, reads)
                
                    return phased_id, filename, length, reads
                    
                else:
                    log_print(f"NOTE:\tSkipping additition to MycoMmap; Reads in Consensus (RiC) <{summary_ric_cutoff}: {reads}")
                    return None, None, None, None 
            else:
                log_print("ERROR:\tPhased FASTA file does not have the expected format.")
                return None, None, None, None
    

def process_medaka_folder(medaka_folder, summary_folder, summary_ric_cutoff):
    """
    Processes a given 'medaka' folder. Extracts the medaka ID and processes the consensus FASTA file 
    to get relevant data like supporting reads count and sequence length.

    Args:
        medaka_folder (str): The path of the 'medaka' folder to be processed.
        summary_folder (str): The path to the summary folder.

    Returns:
        tuple: A tuple containing filename, length, and reads information.

    Notes:
        - The function looks for a 'consensus.fasta' file within the given folder.
        - Extracts the medaka ID from the folder name and supporting reads and length from the consensus sequence name.

    Example:
        filename, length, reads = process_medaka_folder("/path/to/medaka_folder")
    """
    # Split the medaka_folder string by "_" and save the last item in the list as "medaka_id" variable
    medaka_id = medaka_folder.split("_")[-1]

    # Generate a path to the consensus file based on os.path.join(medaka_folder + "consensus.fasta")
    consensus_path = os.path.join(medaka_folder, "consensus.fasta")
    print(consensus_path)
    reads = None
    length = None
    filename = None
    parent_folder_name = None
    
    # Use Biopython to open the fasta file
    if not os.path.exists(consensus_path):
        log_print(f"ERROR:\tConsensus FASTA file not found in {medaka_folder}")
        return None, None, None, None
    with open(consensus_path, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            # Check to see if the sequence name is as follows: "consensus_cl_id_{medaka_id}_total_supporting_reads_[0-9*]"
            if f"consensus_cl_id_{medaka_id}_total_supporting_reads" in record.description:
                # Save the last item in the seuqence name split by "_" list as the "reads" variable
                reads = record.description.split("_")[-1]
                print(reads)
                if int(reads) >= int(summary_ric_cutoff):
                    # Get the length of the sequence and save it as "length" variable
                    length = len(record.seq)
                    
                    # Generate the variable parent_folder = Path(consensus_path).parent.parent.name.replace("sample_","")
                    parent_folder_name = Path(consensus_path).parent.parent.name.replace("sample_", "")
    
                    # Generate the variable filename = f'{parent_folder_name}-{medaka_id}'
                    filename = f'{parent_folder_name}-{medaka_id}'
                    
                    updated_fasta_path = os.path.join(summary_folder, f'{filename}-RiC{reads}.fasta')
                    
                    if os.path.exists(updated_fasta_path):
                        log_print("NOTE:\tSkipping move")
                    else:
                        move_updated_fasta(updated_fasta_path, summary_folder, consensus_path, parent_folder_name, filename, reads)
                        
                    return medaka_id, filename, length, reads
                else:
                    log_print(f"NOTE:\tSkipping additition to MycoMmap; Reads in Consensus (RiC) <{summary_ric_cutoff}: {reads}")
                    return None, None, None, None 
            else:
                log_print(f"ERROR:\tConsensus FASTA file in {medaka_folder} does not have the expected format")
                return None, None, None, None


def move_updated_fasta(updated_fasta_path, summary_folder, consensus_path, parent_folder_name, filename, reads):
    """
    DESCRIPTION NEEDED.

    Args:
        summary_folder (TYPE NEEDED): DESCRIPTION NEEDED.
        consensus_path (TYPE NEEDED): DESCRIPTION NEEDED.
        parent_folder_name (TYPE NEEDED): DESCRIPTION NEEDED.
        filename (TYPE NEEDED): DESCRIPTION NEEDED.
        reads (TYPE NEEDED): DESCRIPTION NEEDED.

    Notes:
        -
        -
        -
            
    Example:
        generate_summary(stats, writer)
    """
    # Wait for a moment before moving the file
    time.sleep(0.5)
    
    # Move consensus_path to updated_fasta_path
    if os.access(consensus_path, os.W_OK):
        shutil.move(consensus_path, updated_fasta_path)
    else:
        log_print(f"ERROR:\tCannot move file {consensus_path} - it may be in use by another process")

    consensus_parent_folder = str(Path(consensus_path).parent)
    fastq_file_index = 1
    for fastq_file in sorted(glob.glob(f'{consensus_parent_folder}/reads_to_consensus_[0-9]*.fastq')):
        if os.path.exists(fastq_file):
            target_fastq_file = os.path.join(summary_folder, 'FASTQ Files', f'{parent_folder_name}-{fastq_file_index}.fastq')
            # shutil.move(fastq_file, target_fastq_file) # TODO
            shutil.copy(fastq_file, target_fastq_file)
            fastq_file_index += 1

def merge_fasta_files(summary_folder, output_file):
    """
    DESCRIPTION NEEDED.

    Args:
        summary_folder (TYPE NEEDED): DESCRIPTION NEEDED.
        output_file (TYPE NEEDED): DESCRIPTION NEEDED.
        
    Notes:
        - 
        -
        -
            
    Example:
        merge_fasta_files(summary_folder, 'summary.fasta')
    """
    with open(os.path.join(summary_folder, output_file), 'w') as f_out:
        for fasta_file in sorted(glob.glob(f'{summary_folder}/*.fasta')):
            with open(fasta_file, 'r') as f:
                f_out.write(f.read())

def generate_summary(stats, writer):
    """
    DESCRIPTION NEEDED.

    Args:
        stats (TYPE NEEDED): DESCRIPTION NEEDED.
        writer (TYPE NEEDED): DESCRIPTION NEEDED.
        
    Notes:
        -
        -
        -
            
    Example:
        generate_summary(stats, writer)
    """
    # Set placeholders for stats
    summary_totals = {'Total Unique Samples': 0,
                      'Total Consensus Sequences': 0,
                      'Total Reads in Consensus Sequences': 0,}
    additional_stats = []
    counted_samples = {}

    # Compile stats
    for row in stats:
        filename = row[0]
        reads = row[2]
        if reads is not None:
            reads = int(reads)
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
        else:
            #print(f"Warning: Reads value is None for filename {filename}. Skipping this row.")
            continue  # Skip processing this row        

    # Create final TSV summary
    writer.writerow([])
    writer.writerow(['Total Unique Samples', summary_totals['Total Unique Samples']])
    writer.writerow(['Total Consensus Sequences', summary_totals['Total Consensus Sequences']])
    writer.writerow(['Total Reads in Consensus Sequences', summary_totals['Total Reads in Consensus Sequences']])
    writer.writerow([])
    for row in additional_stats:
        writer.writerow(row)

# Main Summary Function
def summary_main(source_folder, summary_ric_cutoff, hap_phase_bool):
    """
    DESCRIPTION NEEDED.

    Args:
        source_folder (str): Path to the output of NGSpeciesID.
        hap_phase_bool (bool): True/False on whether Haplotype Phasing was performed.
        
    Notes:
        -
        -
        -
            
    Example:
        summary_main(source_folder, False)
    """
    print(source_folder)
    print(summary_ric_cutoff)
    print(hap_phase_bool)
    
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
        for sample_folder in sorted(glob.iglob(source_folder + '/*', recursive=False)):
            log_print(sample_folder)
            if sample_folder == "/mnt/d/FunDiS/ONT04/combined/__Summary__":
                pass
            elif not os.path.isfile(sample_folder):
                folder_stats = process_sample_folder(sample_folder, summary_folder, summary_ric_cutoff, hap_phase_bool)
                stats = stats + folder_stats
    
        # Generate Final Summary based on the CSV and the compiled Stats
        generate_summary(stats, writer)
        
    # Run final FASTA merging function
    merge_fasta_files(summary_folder, 'summary.fasta')
    log_print('MycoMap Summarize Complete!')

# Wrapper function to run ngsid_prep in a separate thread
def run_summary_prep(source_folder, summary_ric_cutoff, hap_phase_bool):
    """
    Wrapper function to run summary preparation in a separate thread.

    This function initiates the summary preparation process, potentially for generating summaries of data processing. 
    It handles the execution in a separate thread and logs the start of the process and any errors encountered.

    Args:
        source_folder (str): Path to the output of NGSpeciesID.
        hap_phase_bool (bool): True/False on whether Haplotype Phasing was performed.

    Notes:
        - The actual summary preparation logic (e.g., summary_main) needs to be implemented or called within this function.
        - Catches and logs exceptions that may occur during summary preparation.
            
    Example:
        run_summary_prep(INPUT_FOLDER, "40", False)
    """
    log_print("MycoMap Summary preparation started...\n")
    
    try:
        summary_main(source_folder, summary_ric_cutoff, hap_phase_bool)
    except Exception as e:
        print( f"ERROR:\t{str(e)}")

# Debug/Testing Area
if __name__ == "__main__":
    # INPUT_FOLDER = "D:/FunDiS/ONT04/combined"
    INPUT_FOLDER = "/mnt/d/FunDiS/ONT02/combined"

    # Generate log file with the desired behavior
    initialize_logging_environment(INPUT_FOLDER)
       
    run_summary_prep(INPUT_FOLDER, "40", True)
