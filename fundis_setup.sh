#!/bin/bash
# FunDiS Nanopore Barcoding Pipeline Installation Bash Script

# Check if Conda is installed
if ! command -v conda &> /dev/null; then
    echo -e "\e[31mConda is not installed. Please install conda before running this script.\e[0m"
    # TODO: Add command for user to run and install conda
    exit 1
else
    # Update Conda
    echo -e "\e[36mUpdating Conda...\e[0m"
    conda update -n base -c defaults conda -y
fi

# Install and update Mamba in the base environment
echo -e "\e[36mInstalling and updating mamba...\e[0m"
conda install mamba -n base -c conda-forge -y
conda update mamba -n base -c conda-forge -y

# Update pip in the base environment
echo -e "\e[36mUpdating pip...\e[0m"
conda run -n base pip install --upgrade pip

# Create a custom Conda environment with Python 3.8 for main FunDiS Nanopore Barcoding Pipeline GUI
echo -e "\e[36mCreating custom Conda environment 'fundis_env'...\e[0m"
conda create -n fundis_env python=3.8 -y

# Install necessary python libraries with mamba
echo -e "\e[36mInstall necessary python libraries into fundis_env with mamba...\e[0m"
conda run -n fundis_env mamba install -y -c bioconda -c conda-forge -c agbiome -c prkrekel openblas==0.3.3 biopython==1.81 samtools==1.18 minimap2==2.26 bcftools==1.17 bwa==0.7.17 whatshap==2.1 spoa==4.1.3 racon==1.5.0 pyvcf==0.6.8 termcolor=2.3.0 gdown==4.7.1

# Load the Conda environment and install NGSpeciesID
echo -e "\e[36mLoad the Conda environment and install NGSpeciesID into fundis_env...\e[0m"
conda run -n fundis_env pip install NGSpeciesID

# Create a custom Conda environment with Python 3.10.8 for Medaka branch
echo -e "\e[36mCreating custom Conda environment 'medaka' with Python 3.10.8...\e[0m"
conda create -n medaka python=3.10.8 -y
conda run -n medaka mamba install -y -c bioconda -c conda-forge -c agbiome -c prkrekel medaka==1.11.1 termcolor=2.3.0

echo -e "\e[33mRun 'conda activate fundis_env' to start working\e[0m"
echo -e "\e[32mFunDiS Nanopore Barcoding Pipeline set-up complete!\e[0m"
echo "@author: ian.michael.bollinger@gmail.com/researchconsultants@critical.consulting"
echo "This is a modified protocol by Stephen Douglas Russell and paid for by the Fungal Diversity Survey - FunDiS"
echo "Protocol Link: https://www.protocols.io/view/primary-data-analysis-basecalling-demultiplexing-a-dm6gpbm88lzp/v3?step=3"
