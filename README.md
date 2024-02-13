# FunDiSPipe: Fungal Diversity Survey Pipeline

FunDiSPipe is a comprehensive bioinformatics pipeline designed for the Fungal Diversity Survey (FunDiS), specifically tailored for analyzing fungal ITS data from Oxford Nanopore Technologies sequencing. This pipeline streamlines the process from sequencing data to species identification and summarization. This is the main Graphical User Interface for a modified protocol devloped by Stephen Douglas Russell (https://www.protocols.io/view/primary-data-analysis-basecalling-demultiplexing-a-dm6gpbm88lzp/v3?step=3); this pipeline was paid for by the Fungal Diversity Survey (FunDiS).

## Prerequisites

This application is designed to be run on a Linux/WSL environment and requires the following Python libraries:
- numpy==1.26.3
- pandas==2.2.0
- matplotlib==3.8.2
- seaborn==0.13.1
- pyyaml==6.0.1
- statsmodels==0.14.1
- pandas==2.0.3
- openblas==0.3.3
- biopython==1.81
- samtools==1.18
- minimap2==2.26
- bcftools==1.17
- bwa==0.7.17
- whatshap==2.1
- spoa==4.1.3
- racon==1.5.0
- psutil==5.9.8
- blast==2.15.0
- pyvcf==0.6.8
- fastqc==0.12.1

The application also relies on the following tools:
- NGSpeciesID (https://github.com/ksahlin/NGSpeciesID)
- medaka (https://github.com/nanoporetech/medaka)

Note: Running fundis_setup.sh attempts to provide the setup to install all python libraries and dependencies.

## Installation

To install FunDiSPipe, follow these steps after cloning the GitHub repository:
```bash
bash ./fundis_setup.sh
```

Close this terminal; open a new one and run the following commands:
```bash
# base installs
mamba install -y -c bioconda -c conda-forge -c agbiome -c prkrekel numpy==1.26.3 pandas==2.2.0 matplotlib==3.8.2 seaborn==0.13.1 pyyaml==6.0.1 statsmodels==0.14.1
```

```bash
# fundis_env installs
conda activate fundis_env && mamba install -y -c bioconda -c conda-forge -c agbiome -c prkrekel pandas==2.0.3 openblas==0.3.3 biopython==1.81 samtools==1.18 minimap2==2.26 bcftools==1.17 bwa==0.7.17 whatshap==2.1 spoa==4.1.3 racon==1.5.0 psutil==5.9.8 blast==2.15.0 pyvcf==0.6.8 fastqc==0.12.1 && pip install NGSpeciesID && conda deactivate
```

```bash
# medaka env installs
conda activate medaka && mamba install -y -c bioconda -c conda-forge -c agbiome -c prkrekel medaka==1.11.1 && conda deactivate
```

## Modules and Their Functionalities

1. **GUI (FunDiS_GUI.py)**:
   - Acts as the central interface for the pipeline.
   - Facilitates file selection, process initiation, and result visualization.
   - Integrates other modules for a seamless workflow.

2. **Mini-Barcoder (FunDiS_Minibar.py)** (https://github.com/calacademy-research/minibar): 
   - Prepares `.fastq.gz` files for species identification.
   - Extracts and processes sequences from raw data.
   - Essential for initial data preparation and quality control.

3. **NGSpeciesID (FunDiS_NGSpeciesID.py)** (https://github.com/ksahlin/NGSpeciesID):
   - Identifies species from processed sequencing data.
   - Utilizes advanced algorithms for accurate species matching.
   - Outputs detailed reports on identified species and their characteristics.

4. **Haplotype Phaser (FunDiS_hap_phase.py)**:
   - Resolves haplotype variations in sequencing data.
   - Enhances species identification accuracy.
   - Critical for detailed genetic analysis and research.

5. **MycoMap Summarizer (MycoMap_Summarize.py)** (https://www.protocols.io/view/primary-data-analysis-basecalling-demultiplexing-a-dm6gpbm88lzp/v3?step=3):
   - Aggregates results from the entire pipeline.
   - Produces comprehensive summary reports for analysis and interpretation.
   - Simplifies data review and sharing.

## Inputs and Outputs

- **Input**: `.fastq.gz` file containing Oxford Nanopore Guppy Basecalled sequences.
- **Outputs**:
  - Processed and quality-checked sequence data.
  - Species identification reports and detailed analysis.
  - Summarized outputs and aggregated data for further study.

## Usage

Start by navigating the GUI to select your input files. Here’s a brief guide on using each module:

- **GUI**: Launch the GUI script to access the pipeline's functionalities. The interface is intuitive and guides you through the process.
- **Mini-Barcoder**: After selecting your `.fastq.gz` file, this module will prepare it for the NGSpeciesID analysis.
- **NGSpeciesID**: Once the data is prepped, use this module for species identification. The output will include detailed species information.
   - **Haplotype Phaser**: An advanced setting that will perform a dual analysis and preserve haplotypes with IUPAC ambiguities and lower-case letters for insertions/deletions.
- **MycoMap Summarizer**: Finally, to aggregate and summarize your results, use this module. It consolidates the data into an easy-to-interpret format.
   - FEATURE PENDING **Haplotype Phaser**: An advanced setting that will use the Phased Haplotype file instead of the traditionally used medaka consensus file.

For detailed instructions and options for each module, refer to the comments and documentation within each script file. These instructions provide guidance on executing the scripts and customizing the analysis to your requirements.

## Contributing

Contributions are welcome. Please follow standard coding practices and clearly document any changes or enhancements.

## License

Please see the LICENSE file in the GitHub repository for detailed licensing information.

## Author

ian.michael.bollinger@gmail.com/researchconsultants@critical.consulting
