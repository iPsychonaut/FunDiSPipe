# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 12:23:15 2023

@author: ian.michael.bollinger@gmail.com/researchconsultants@critical.consulting

This is the Haplotype Phaser developed for a modified protocol by

Stephen Douglas Russell and paid for by the Fungal Diversity Survey (FunDiS).

Protocol Link: https://www.protocols.io/view/primary-data-analysis-basecalling-demultiplexing-a-dm6gpbm88lzp/v3?step=3
"""

# Base Python Imports
import os, re, subprocess, multiprocessing, argparse

# Required Python Imports
from termcolor import cprint
from datetime import datetime
import pandas as pd
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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

def run_subprocess_cmd(cmd_list, shell_check):
    if isinstance(cmd_list, str):
        log_print(f"CMD:\t{cmd_list}")    
        process = subprocess.run(cmd_list, text=True, shell=shell_check, capture_output=True)
        if process.returncode != 0:
            log_print(f"ERROR:\t{process.stderr}")
        else:
            log_print(f"PASS:\tSuccessfully processed command: {cmd_list}")
    else:        
        log_print(f"CMD:\t{' '.join(cmd_list)}")    
        process = subprocess.run(cmd_list, text=True, shell=shell_check, capture_output=True)
        if process.returncode != 0:
            log_print(f"ERROR:\t{process.stderr}")
        else:
            log_print(f"PASS:\tSuccessfully processed command: {' '.join(cmd_list)}")


def run_blast(query_file, db_file, add_headers=True):
    log_print("BLASTING query against database...")
    blast_output = query_file.replace('.fasta', '_blast_results.tsv')
    blastn_cline = NcbiblastnCommandline(query=query_file, db=db_file, evalue=0.001, outfmt=6, out=blast_output)
    os.system(str(blastn_cline))

    if add_headers:
        headers = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                   "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
        with open(blast_output, 'r') as original: data = original.read()
        with open(blast_output, 'w') as modified: modified.write("\t".join(headers) + "\n" + data)
    
    log_print("PASS:\tSuccessfully BLAST'd query against database")
    return blast_output

def append_sequence_to_fasta(source_fasta, destination_fasta):
    with open(source_fasta, 'r') as source_file:
        source_content = source_file.read()

    with open(destination_fasta, 'a') as dest_file:
        dest_file.write(source_content)
        dest_file.write("\n")  # Add a newline for separation

def find_latest_paf_file(directory):
    log_print("Searching for the latest PAF file...")
    paf_files = [f for f in os.listdir(directory) if re.match(r'read_alignments_it_\d+\.paf', f)]
    if not paf_files:
        raise ValueError("ERROR:\tNo PAF files found in the directory")

    try:
        highest_iteration = max(int(re.search(r'(\d+)', f).group()) for f in paf_files)
    except Exception as e:
        log_print(f"ERROR:\tError finding the highest iteration PAF file: {e}")
        raise

    latest_paf_file = f'read_alignments_it_{highest_iteration}.paf'
    log_print(f"PASS:\tLatest PAF file found: {latest_paf_file}")
    return os.path.join(directory, latest_paf_file)

def check_and_add_rg_to_bam(bam_file, reference_seq_ids):
    modified_bam_file = bam_file.replace(".bam", "_rg_added.bam")
    existing_rgs = set()
    cmd = f"samtools view -H {bam_file}"
    process = subprocess.run(cmd, shell=True, text=True, stdout=subprocess.PIPE)
    if process.returncode != 0:
        log_print(f"ERROR:\t{process.stderr}")
        return None
    for line in process.stdout.split('\n'):
        if line.startswith('@RG'):
            rg_id = line.split('\t')[1].replace("ID:", "")
            existing_rgs.add(rg_id)

    for seq_id in reference_seq_ids:
        if seq_id not in existing_rgs:
            log_print(f"Adding Read Group for {seq_id} to BAM file")
            cmd = f"samtools addreplacerg -r \"@RG\\tID:{seq_id}\\tSM:{seq_id}\" -o {modified_bam_file} {bam_file}"
            process = subprocess.run(cmd, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if process.returncode != 0:
                log_print(f"ERROR:\t{process.stderr}")
                return None
            bam_file = modified_bam_file  # Update BAM file for next iteration

    if bam_file != modified_bam_file:
        log_print("NOTE:\tNo Read Groups added. Using original BAM file")
        return bam_file

    log_print(f"PASS:\tRead Groups added to BAM file: {modified_bam_file}")
    return modified_bam_file

def extract_read_names_from_paf(paf_file, min_matching_bases=10, max_divergence=0.01):
    log_print(f"Extracting read names from PAF file: {paf_file}...")
    read_names = set()
    with open(paf_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) < 16:
                log_print(f"Skipping line with unexpected format: {line}")
                continue

            # Find the 'cm:i:' and 'dv:f:' fields
            cm_field = next((field for field in parts if field.startswith('cm:i:')), None)
            dv_field = next((field for field in parts if field.startswith('dv:f:')), None)

            if not cm_field or not dv_field:
                log_print(f"ERROR:\tRequired fields not found in line: {line}")
                continue

            try:
                matching_bases = int(cm_field.split(':')[2])
                divergence = float(dv_field.split(':')[2])
            except ValueError as e:
                log_print(f"ERROR:\tError parsing line: {line}\nError: {e}")
                continue

            if matching_bases >= min_matching_bases and divergence <= max_divergence:
                read_name = parts[0]
                read_names.add(read_name)
    log_print(f"PASS:\tNumber of reads extracted: {len(read_names)}")
    return read_names

def extract_sequences_from_fastq(fastq_file, read_names):
    log_print(f"Number of read IDs provided: {len(read_names)}...")
    sequences = []
    not_found_ids = []

    for record in SeqIO.parse(fastq_file, "fastq"):
        if record.id in read_names:
            sequences.append(record)
        else:
            not_found_ids.append(record.id)

    if not_found_ids:
        log_print(f"ERROR:\tRead IDs not found in the set: {not_found_ids[:10]} (first 10)")
    else:
        log_print("PASS:\tAll read IDs were successfully found.")

    log_print(f"PASS:\tNumber of sequences extracted: {len(sequences)}")

    return sequences

def determine_medaka_consensus_seqs(fastq_file, seq):
    log_print("Determinine sequences used to generate medaka consensus...")
    directory = f"{fastq_file.split('reads')[0]}/medaka_cl_id_{seq}"
    output_file = fastq_file.replace(".fastq","_medaka_filtered.fastq") 
    read_names_bam = os.path.join(directory,"calls_to_draft.bam")
    read_names_file = read_names_bam.replace(".bam",".txt")

    # Running samtools view and capturing output
    cmd = f"samtools view {read_names_bam}"
    process = subprocess.run(cmd, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if process.returncode != 0:
        log_print(f"ERROR:\t{process.stderr}")
        return None

    # Writing output to file
    with open(read_names_file, 'w') as file:
        file.write('\n'.join(line.split('\t')[0] for line in process.stdout.splitlines()))

    # Debug: Print the first few read names from the BAM file
    first_few_reads = []
    with open(read_names_file, 'r') as file:
        try:
            for _ in range(10):
                line = next(file).strip()
                first_few_reads.append(line)
        except StopIteration:
            pass  # Handle the case where the file has fewer than 10 lines

    log_print("NOTE:\tFirst few read names extracted from BAM file:", first_few_reads)

    sequences = extract_sequences_from_fastq(fastq_file, read_names_file)
    if not sequences:
        log_print("ERROR:\tNo sequences were extracted. There might be a mismatch in read names.")

    # Debug: Print the first few read names from the FASTQ file
    with open(fastq_file, "r") as file:
        for _ in range(10):
            log_print(file.readline().strip())  # Print the title line of each read

    SeqIO.write(sequences, output_file, "fastq")
    log_print(f"PASS:\tExtracted sequences written to {output_file}")
    return output_file

def determine_racon_consensus_seqs(fastq_file, seq):
    log_print("Determinine sequences used to generate racon consensus...")
    directory = f"{fastq_file.split('reads')[0]}/racon_cl_id_{seq}"
    output_file = fastq_file.replace(".fastq", "_racon_filtered.fastq")
    try:
        latest_paf_file = find_latest_paf_file(directory)
        if not os.path.exists(latest_paf_file):
            log_print(f"ERROR:\tLatest PAF file not found: {latest_paf_file}")
            return None

        read_names = extract_read_names_from_paf(latest_paf_file)
        if not read_names:
            log_print("ERROR:\tNo read names extracted from PAF file. Aborting racon consensus sequence extraction.")
            return None

        if not os.path.exists(fastq_file):
            log_print(f"ERROR:\tFASTQ file not found: {fastq_file}")
            return None

        sequences = extract_sequences_from_fastq(fastq_file, read_names)
        if not sequences:
            log_print("ERROR:\tNo sequences extracted from FASTQ file. Aborting racon consensus sequence extraction.")
            return None

        SeqIO.write(sequences, output_file, "fastq")
        log_print(f"PASS:\tExtracted sequences written to {output_file}")
        return output_file
    except Exception as e:
        log_print(f"ERROR:\tAn error occurred in determine_racon_consensus_seqs: {e}")
        return None

def combine_fastq_files(racon_fastq, medaka_fastq, output_fastq):
    log_print("Combining racon and medaka consensus sequences into a single FASTQ file without duplicates...")
    combined_sequences = {}

    # Read sequences from Racon FASTQ file
    with open(racon_fastq, 'r') as file:
        for title, seq, qual in FastqGeneralIterator(file):
            combined_sequences[title.split()[0]] = (seq, qual)

    # Read sequences from Medaka FASTQ file
    with open(medaka_fastq, 'r') as file:
        for title, seq, qual in FastqGeneralIterator(file):
            combined_sequences[title.split()[0]] = (seq, qual)

    # Write combined sequences to a new FASTQ file
    with open(output_fastq, 'w') as file:
        for title, (seq, qual) in combined_sequences.items():
            file.write(f"@{title}\n{seq}\n+\n{qual}\n")

    log_print(f"PASS:\tCombined FASTQ file written to {output_fastq}")

def create_phased_fasta(reference_seq_file, phased_vcf_file, seq, consensus_seq_count):
    log_print("Generating phased FASTA file from phased VCF data...")
    phased_fasta_file = phased_vcf_file.replace(".vcf", ".fasta")

    # Define the IUPAC ambiguity codes
    iupac_ambiguity = {('A', 'C'): 'M', ('A', 'G'): 'R', ('A', 'T'): 'W',
                       ('C', 'G'): 'S', ('C', 'T'): 'Y', ('G', 'T'): 'K',
                       ('C', 'A'): 'M', ('G', 'A'): 'R', ('T', 'A'): 'W',
                       ('G', 'C'): 'S', ('T', 'C'): 'Y', ('T', 'G'): 'K'}

    # Load the reference sequence
    reference_seq_record = next(SeqIO.parse(reference_seq_file, "fasta"))
    reference_seq_list = list(reference_seq_record.seq)

    # Load VCF data into a DataFrame
    vcf_df = pd.read_csv(phased_vcf_file, sep='\t', comment='#', 
                         names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE_DATA"])
    
    # Iterate through VCF DataFrame
    for index, row in vcf_df.iterrows():
        pos = int(row['POS']) - 1
        ref = row['REF']
        alt = row['ALT'].split(',')[0]  # Consider the first alternate allele

        if len(ref) == 1 and len(alt) == 1:  # SNP
            ambiguity_code = iupac_ambiguity.get((ref, alt), alt)
            reference_seq_list[pos] = ambiguity_code
            log_print(f"PASS:\tSNP at Position: {pos + 1}, Seq.Nucleotide: {reference_seq_list[pos]} = VCF.Nucleotides: ({ref}, {alt}) => {ambiguity_code}")
        elif len(ref) != len(alt):  # Indel
            if len(ref) < len(alt):  # Insertion
                insert_seq = alt[len(ref):].lower()
                reference_seq_list.insert(pos + 1, insert_seq)
                log_print(f"PASS:\tInsertion at Position: {pos + 1}, Inserted Nucleotides: {insert_seq}")
            else:  # Deletion
                for i in range(len(ref) - len(alt)):
                    reference_seq_list[pos + len(alt) + i] = reference_seq_list[pos + len(alt) + i].lower()
                log_print(f"PASS:\tDeletion at Position: {pos + 1}, Lowercased Nucleotides: {ref[len(alt):].lower()}")
        else:
            log_print(f"ERROR:\tPosition: {row['POS']}, VCF.Nucleotides: ({row['REF']}, {row['ALT']})")
            
    # Convert the list back to a string to get the updated sequence
    updated_reference_sequence = ''.join(reference_seq_list)

    # Save the new sequence into a new fasta file with the same ID but with "_haplotype_phased"
    new_record = SeqRecord(Seq(updated_reference_sequence), id=reference_seq_record.id + "_haplotype_phased",
                           description=f"Haplotype Phased total supporting reads {consensus_seq_count}")

    # Save the new sequence into a new fasta file
    with open(phased_fasta_file, 'w') as output_handle:
        SeqIO.write(new_record, output_handle, 'fasta')

    log_print(f"PASS:\tSuccessfully generated haplotype phased FASTA file: {phased_fasta_file}")

    return phased_fasta_file

def get_single_sequence_length(fasta_file):
    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            return len(record.seq)
    return None  # In case the FASTA file is empty or not properly formatted

def get_sequence_ids_from_fasta(fasta_file):
    sequence_ids = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence_ids.append(record.id)
    return sequence_ids

def phase_consensus_seqs(seq, racon_consensus_file, medaka_consensus_file, ngsid_output_dir, cpu_threads, sanger_sequence_file=None):
    log_print("Phasing BLAST output...")
    sequence_threshold = 30
    consensus_fastq_file = f"{ngsid_output_dir}/reads_to_consensus_{seq}.fastq"
    racon_consensus_seqs_fastq = determine_racon_consensus_seqs(consensus_fastq_file, seq)
    medaka_consensus_seqs_fastq = determine_medaka_consensus_seqs(consensus_fastq_file, seq)
    
    # Combine unique racon and medaka fastq sequences into one
    if racon_consensus_seqs_fastq and medaka_consensus_seqs_fastq:                    
        combine_fastq_files(racon_consensus_seqs_fastq, medaka_consensus_seqs_fastq, consensus_fastq_file)
        fastqc_dir = consensus_fastq_file.replace(".fastq","_fastqc")
        if not os.path.exists(fastqc_dir):
            os.makedirs(fastqc_dir)
        consensus_seq_count = 0
        for record in SeqIO.parse(consensus_fastq_file, "fastq"):
            consensus_seq_count += 1
        
        if consensus_seq_count < sequence_threshold:
            log_print(f"NOTE:\tConsensus sequence count did not meet required threshold of 40 unique sequences: {consensus_seq_count}")
            return None
        
        # Run FastQC to check for desireable qualities in the consensus sequences
        run_subprocess_cmd(["fastqc",
                          consensus_fastq_file,
                          "--threads",
                          str(cpu_threads),
                          "-o",
                          fastqc_dir],
                           shell_check=False)
        """
        TODO: Perform FastQC output check for desirable qualities:
            Adapter content: Check if the adapter content is within acceptable limits.
            Overrepresented sequences: Accept the presence of numerous overrepresented sequences as normal.
            Duplication levels: Check if there's a high duplication level of a specific sequence.
            Sequence length distribution: Verify that the sequence length is within ±15% of the mean.
            Per-base N content: Ensure that it's low at all data points.
            Per-base sequence quality: Check that the quality "passes" the FastQC threshold.
        """

        # Convert FASTQ reads into VCF format for Haplotype Phasing

         # Run Minimap2 using sanger sequence if available, if not use medaka consensus as reference fasta
        aligned_sam_file = consensus_fastq_file.replace(".fastq","_aligned.sam")
        reference_seq_file = sanger_sequence_file if os.path.exists(sanger_sequence_file) else medaka_consensus_file

        run_subprocess_cmd(f"samtools faidx {reference_seq_file}", shell_check=True)
        run_subprocess_cmd(f"minimap2 -ax map-ont {reference_seq_file} {consensus_fastq_file} > {aligned_sam_file}", shell_check=True)
        
         # Convert SAM to BAM
        aligned_bam_file = aligned_sam_file.replace("_aligned.sam", "_aligned.bam")
        run_subprocess_cmd(f"samtools view -bS {aligned_sam_file} > {aligned_bam_file}", shell_check=True)
        
         # Extract sequence IDs from reference FASTA file
        reference_seq_ids = get_sequence_ids_from_fasta(reference_seq_file)
    
        # Check and add read groups to BAM file
        aligned_bam_file_with_rg = check_and_add_rg_to_bam(aligned_bam_file, reference_seq_ids)
        if aligned_bam_file_with_rg is None:
            log_print("ERROR:\tFailed to add read groups to BAM file.")
            return None
        
         # Sort and index the modified BAM file
        sorted_bam_file = aligned_sam_file.replace("_aligned.sam", "_sorted.bam")
        run_subprocess_cmd(f"samtools sort {aligned_bam_file_with_rg} -o {sorted_bam_file}", shell_check=True)
        run_subprocess_cmd(f"samtools index {sorted_bam_file}", shell_check=True)
        
        # Run bcftools and WhatsHap with the modified BAM file
        variants_vcf_file = consensus_fastq_file.replace(".fastq","_variants.vcf")
        run_subprocess_cmd(f"bcftools mpileup -Ou -f {reference_seq_file} {sorted_bam_file} | bcftools call -mv -Ov -o {variants_vcf_file}", shell_check=True)

        phased_vcf_file = consensus_fastq_file.replace(".fastq","_phased.vcf")
        run_subprocess_cmd(f"whatshap phase -o {phased_vcf_file} -r {reference_seq_file} {variants_vcf_file} {sorted_bam_file}", shell_check=True)
        
        # Create phased FASTA with ambiguities from phased_vcf
        phased_fasta_file = create_phased_fasta(reference_seq_file, phased_vcf_file, seq, consensus_seq_count)
        
        return phased_fasta_file
    else:
        log_print("ERROR:\tError occurred in processing consensus sequences. Could not combine FASTQ files.")
        return None

def phase_haplotypes(ngsid_output_dir, cpu_threads):
    # Initialize an empty set to store all the extracted numbers
    all_numbers = set()

    # Regular expression pattern to match "cl_id_" followed by digits
    pattern = r'cl_id_(\d+)'

    # Iterate through the directory
    for root, dirs, files in os.walk(ngsid_output_dir):
        for dir_name in dirs:
            # Check if "racon" is in the folder name or "medaka"
            if 'racon' in dir_name or 'medaka' in dir_name:
                match = re.search(pattern, dir_name)
                if match:
                    all_numbers.add(int(match.group(1)))

    # Convert the set to a sorted list if needed
    seq_list = sorted(list(all_numbers))

    for seq in seq_list:
        concatenated_fasta = f"{ngsid_output_dir}/concatentated_reads_{seq}.fasta"

        # Blast racon consensus against medaka consensus
         # Use Medaka consensus as db
        medaka_consensus_file = f"{ngsid_output_dir}/medaka_cl_id_{seq}/consensus.fasta"
        medaka_db_path = medaka_consensus_file.replace(".fasta", "")
        run_subprocess_cmd(["makeblastdb",
                            "-in", medaka_consensus_file,
                            "-dbtype", "nucl",
                            "-out", medaka_db_path],
                           shell_check=False)
        
         # Add Medaka consensus sequence to the concatenated sequence fasta file
        append_sequence_to_fasta(medaka_consensus_file, concatenated_fasta)

         # Blast Racon consensus against Medaka consensus
        racon_consensus_file = f"{ngsid_output_dir}/racon_cl_id_{seq}/consensus.fasta"
        blast_output_tsv = run_blast(racon_consensus_file, medaka_db_path, add_headers=True)
        blast_output_df = pd.read_csv(blast_output_tsv, sep="\t")

         # Add Racon consensus sequence to the concatenated sequence fasta file
        append_sequence_to_fasta(racon_consensus_file, concatenated_fasta)

         # Include Sanger sequence if available
        sanger_sequence_file = f"{ngsid_output_dir}/{ngsid_output_dir.split('/')[-1].split('_iNat')[0]}_sanger.fasta"
        log_print(sanger_sequence_file)
        if os.path.exists(sanger_sequence_file):
            log_print("Sanger sequence exists")
            # Add Sanger sequence to concatenated sequence fasta file
            append_sequence_to_fasta(sanger_sequence_file, concatenated_fasta)
        else:
            log_print("No Sanger sequence found")
        
        # Determine if phasing is needed based on BLAST output
        identity_threshold=100.0
        evalue_threshold= 0 #1e-5
        
        medaka_seq_len = get_single_sequence_length(medaka_consensus_file)
        racon_seq_len = get_single_sequence_length(racon_consensus_file)
        log_print(f"PASS:\tConsensus Sequence Lenghts: medaka: {medaka_seq_len}, racon: {racon_seq_len}")
        
        for index, row in blast_output_df.iterrows():
            log_print(row)
        
        for index, row in blast_output_df.iterrows():
            if row['pident'] < identity_threshold or row['evalue'] > evalue_threshold:
                phase_fasta_file = phase_consensus_seqs(seq, racon_consensus_file, medaka_consensus_file, ngsid_output_dir, cpu_threads, sanger_sequence_file=sanger_sequence_file)
                return phase_fasta_file
            else:
                log_print(f"NOTE:\tNo discernable differences exist between the Medaka and Racon files at the following thresholds: %ID: {identity_threshold}, e-value: {evalue_threshold}")
                phase_fasta_file = phase_consensus_seqs(seq, racon_consensus_file, medaka_consensus_file, ngsid_output_dir,cpu_threads, sanger_sequence_file=sanger_sequence_file)
                return phase_fasta_file
                #return None
            
if __name__ == "__main__":
    # Set Default values
    default_ngsid_output_dir = "/mnt/d/FunDiS/combined_NGSID/sample_HS_ONT02_01_41_HAY-F-000397_iNat148667504_Minores"   
    # default_ngsid_output_dir = "/mnt/d/FunDiS/combined_NGSID/sample_HS_ONT02_01_16_HAY-F-000306_iNat147376929_Xerocomellus"
    # default_ngsid_output_dir = "/mnt/d/FunDiS/combined_NGSID/sample_HS_ONT02_01_15_HAY-F-000312_iNat147376930_Coprinellus"
    default_percent_resources = 0.8
    
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Phase an NGSpeciesID Folder.")
    parser.add_argument("--input_folder", "-i",
                        type=str, default=default_ngsid_output_dir,
                        help="Path to the NGSPeciesID folder to parse.")
    parser.add_argument("--percent_resources", "-pr",
                        type=float, default=default_percent_resources,
                        help="Percentage of CPU resources to use.")

    # Parse arguments
    args = parser.parse_args()
    NGSID_OUTPUT_DIR = args.input_folder
    PERCENT_RESOURCES = args.percent_resources

    # Get number of total cpu threads available on system and calculate the number to use
    total_cpu_threads = multiprocessing.cpu_count()
    # cpu_threads = int(round(total_cpu_threads * PERCENT_RESOURCES, 0))
    cpu_threads = 10

    # Phase an NGSID Folder
    phased_fasta_file = phase_haplotypes(default_ngsid_output_dir, cpu_threads)
