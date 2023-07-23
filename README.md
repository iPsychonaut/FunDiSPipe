
# FunDiS Pipeline

FunDiS Pipeline is a suite of scripts intended to streamline the processing of Next-Generation Sequencing (NGS) data. The scripts can be run individually or as a whole to form a complete pipeline.

## Prerequisites

This application is designed to be run on a Linux/WSL environment and requires the following Python libraries:

- psutil
- tqdm
- pandas
- pysam
- biopython
- multiprocessing
- math
- queue
- glob
- shutil

The application also relies on the following tools:

- NGSpeciesID
- bwa
- samtools
- bcftools
- whatshap
- medaka
- openblas
- spoa

Note: The application checks for the required Python libraries and tools during execution and attempts to install any missing dependencies.

## Running the Pipeline

Each module of the pipeline can be run individually or as a whole.

### Running the Whole Pipeline

To run the whole pipeline, use the `fundis_main.py` script. For example:

```
python /path/to/fundis_main.py -i /path/to/input.fastq -x /path/to/index.txt -t /path/to/primers.txt -p 80
```

### Running Individual Modules

Each module can also be run individually. Here's what each module does:

- **fundis_minibar_ngsid.py**: This script processes the input FASTQ file with MiniBar and NGSpeciesID. MiniBar is a tool for demultiplexing barcoded read data and NGSpeciesID is a tool used for the identification of specimens in NGS datasets. The script starts by checking the operating system, installing missing libraries, and setting up the working environment. It then moves on to demultiplexing and identifying species from the input FASTQ data. The results are output in a directory specified by the user.

```
python /path/to/fundis_minibar_ngsid.py -i /path/to/input.fastq -x /path/to/index.txt -t /path/to/primers.txt -p 80
```

- **fundis_haplotype_phaser.py**: This script takes the output from the `fundis_minibar_ngsid.py` script and phases the haplotypes for each sample. Phasing is the process of determining the specific set of variants found on each physical copy of a particular gene or genomic region. The phased haplotypes are output in the NGSpeciesID output directory.

```
python /path/to/fundis_haplotype_phaser.py -i /path/to/input_dir -p 80
```

- **fundis_summarize.py**: This script summarizes the output from the `fundis_haplotype_phaser.py` script. It provides a summary of the results, including counts of unique samples, total consensus sequences, and total reads in consensus sequences. It also copies and updates the names of all FASTQ and consensus FASTA files. The results are output in a summary directory named after the source directory.

```
python /path/to/fundis_summarize.py -i /path/to/input_dir -p 80
```

## Arguments

- `-i`, `--input_fastq` or `--input_dir`: Path to the FASTQ file containing ONT nrITS data or path to the directory containing the data.
- `-t`, `--primers_text_path`: Path to Text file containing the Primers used to generate input_fastq.
- `-x`, `--minbar_index_path`: Path to Text file containing the minibar index to parse input_fastq.
- `-p`, `--percent_system_use`: Percent system use written as an integer.

## Author

Ian Michael Bollinger (ian.michael.bollinger@gmail.com/researchconsultants@critical.consulting)
