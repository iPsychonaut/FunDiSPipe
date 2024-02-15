# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 12:23:15 2023

@author: ian.michael.bollinger@gmail.com/researchconsultants@critical.consulting

Can be run on it's own using the following command:
    python MycoMap_Summarize.py -i /path/to/folder/to/summarize -c 40 -p True
    
    -i, --input_folder (str): The path to the folder containing NGSpeciesID output folders needing processing for MycoMap.
    -c, --min_cutoff (int): Minimum number of Reads in Consensus (RiC) allowed for processing for MycoMap.
    -p, --haplotype_phase (bool): True/False of whether or not Haplotype Phasing data should be processed.

This is the MycoMap Summarizer (summarize.py) developed for a modified protocol by

Stephen Douglas Russell and paid for by the Fungal Diversity Survey (FunDiS).

Protocol Link: https://www.protocols.io/view/primary-data-analysis-basecalling-demultiplexing-a-dm6gpbm88lzp/v3?step=3
"""

# Base Python Imports
import csv, os, re, glob, shutil, time, argparse
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
    Moves the updated FASTA file to a specified summary folder and copies associated FASTQ files.

    Args:
        updated_fasta_path (str): The target path for the updated FASTA file after moving.
        summary_folder (str): The path to the summary folder where the FASTA file is stored.
        consensus_path (str): The original path of the consensus FASTA file.
        parent_folder_name (str): The name of the parent folder for organization in the summary folder.
        filename (str): The name of the file, used to generate the final path.
        reads (int): The number of reads in the consensus sequence, used in naming the FASTQ files.

    Notes:
        - This function moves the consensus FASTA file to the updated path and copies related FASTQ files into a
          designated 'FASTQ Files' directory within the summary folder.
        - The function ensures the correct organization and accessibility of sequence data for downstream analysis.
        - Uses a brief delay before moving files to ensure file accessibility, especially in networked or slow file systems.

    Example:
        move_updated_fasta("/path/to/summary_folder/updated.fasta", "/path/to/summary_folder", "/path/to/consensus.fasta", "sample_01", "sample_01-100", 200)
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
    Merges all FASTA files in the summary folder into a single output file.

    Args:
        summary_folder (str): The directory containing individual FASTA files to be merged.
        output_file (str): The name of the file to which the merged FASTA contents will be written.

    Notes:
        - This function is useful for creating a comprehensive FASTA file that aggregates all sequences from a project or experiment.
        - The resulting file is suitable for further analyses or submissions to sequence databases.

    Example:
        merge_fasta_files("/path/to/summary_folder", "combined_sequences.fasta")
    """
    with open(os.path.join(summary_folder, output_file), 'w') as f_out:
        for fasta_file in sorted(glob.glob(f'{summary_folder}/*.fasta')):
            with open(fasta_file, 'r') as f:
                f_out.write(f.read())

def generate_summary(stats, writer):
    """
    Generates a summary table from provided statistics and writes it using a CSV writer.

    Args:
        stats (list): A list of lists, where each sublist contains information about a sequence (e.g., filename, length, reads, ID).
        writer (csv.writer): An initialized CSV writer object to write the summary table.

    Notes:
        - This function calculates the total unique samples, total consensus sequences, and total reads in consensus sequences.
        - It handles both standard and additional stats, segregating entries that contain 'ONT' more than once for special consideration.
        - The summary includes a breakdown by unique sample as well as cumulative totals.

    Example:
        generate_summary(stats, csv_writer)
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
    Main function to generate a summary of sequence data processing results.

    Args:
        source_folder (str): The root directory containing processed sequence data from different samples.
        summary_ric_cutoff (int): The read-in-consensus (RiC) cutoff for including sequences in the summary.
        hap_phase_bool (bool): Indicates whether haplotype phasing was considered in processing.

    Notes:
        - This function orchestrates the creation of a summary folder, collection of statistics from each sample,
          and the generation of a comprehensive summary and merged FASTA file.
        - It is designed to provide an overview and facilitate the review of sequence processing outcomes.

    Example:
        summary_main("/path/to/processed_data", 40, True)
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
    log_print('PASS:\tMycoMap Summarize Complete!')

# Wrapper function to run ngsid_prep in a separate thread
def run_summary_prep(source_folder, summary_ric_cutoff, hap_phase_bool):
    """
    Wrapper function to initiate the summary preparation process in a separate execution context.

    Args:
        source_folder (str): The directory where sequence data processing results are stored.
        summary_ric_cutoff (int): The minimum number of reads in consensus required to include a sequence in the summary.
        hap_phase_bool (bool): Flag to indicate whether haplotype phasing was performed.

    Notes:
        - This function is designed to be run as a separate thread or process, allowing for non-blocking execution of the summary generation.
        - It logs the start of the process and captures any errors that may occur, facilitating debugging and ensuring robust operation.

    Example:
        run_summary_prep("/path/to/sequence_data", 40, False)
    """
    log_print("MycoMap Summary preparation started...\n")
    
    try:
        summary_main(source_folder, summary_ric_cutoff, hap_phase_bool)
    except Exception as e:
        print( f"ERROR:\t{str(e)}")

def parse_arguments():
    """
    Parses command line arguments.

    Returns:
        argparse.Namespace: An object containing the parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description="MycoMap Summarizer")
    parser.add_argument('-i', '--input_folder', type=str, required=True, default="/mnt/d/FunDiS/ONT02/combined",
                        help='The path to the folder containing NGSpeciesID output folders needing processing for MycoMap.')
    parser.add_argument('-c', '--min_cutoff', type=int, required=True, default="40",
                        help='Minimum number of Reads in Consensus (RiC) allowed for processing for MycoMap.')
    parser.add_argument('-p', '--haplotype_phase', type=lambda x: (str(x).lower() in ['true', '1', 'yes']), default=False,
                        required=True, help='True/False of whether or not Haplotype Phasing data should be processed.')

    return parser.parse_args()

if __name__ == "__main__":
    """
    Main argument if run directly from command line, parses inputs as arguments as well.
    """
    args = parse_arguments()

    INPUT_FOLDER = args.input_folder
    SUMMARY_RIC_CUTOFF = args.min_cutoff
    HAPLOTYPE_PHASE = args.haplotype_phase

    # Generate log file with the desired behavior
    initialize_logging_environment(INPUT_FOLDER)
       
    run_summary_prep(INPUT_FOLDER, "40", True)