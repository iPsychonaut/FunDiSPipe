#!/bin/bash
# FunDiS Nanopore Barcoding Pipeline Installation Bash Script

# Function to find and source the Conda configuration file
find_conda_path() {
    local conda_paths=("$HOME/miniforge3" "$HOME/anaconda3" "$HOME/miniconda3" "/opt/miniforge3" "/opt/anaconda3" "/opt/miniconda3" "/usr/local/miniforge3" "/usr/local/anaconda3" "/usr/local/miniconda3")
    for path in "${conda_paths[@]}"; do
        if [ -f "$path/etc/profile.d/conda.sh" ]; then
            echo "$path/etc/profile.d/conda.sh"
            return
        fi
    done
    echo ""
}

# Stop script if any command fails
set -e

# Find and source Conda configuration file
CONDA_PATH=$(find_conda_path)
if [ -n "$CONDA_PATH" ]; then
    source "$CONDA_PATH"
    conda init bash
else
    echo -e "\e[31mConda configuration file not found. Please check your Conda installation.\e[0m"
    exit 1
fi

# Determine OS Type
OS_TYPE=$(uname -s)
if [ "$OS_TYPE" == "Linux" ] || [ "$OS_TYPE" == "Darwin" ]; then
    echo "Detected $OS_TYPE OS."
elif [ "$OS_TYPE" == "MSYS_NT-10.0" ] || [ "$OS_TYPE" == "MSYS" ]; then
    echo -e "\e[31mThis pipeline was not built for Windows, please use WSL\e[0m"
    exit 1
fi

# Check for Conda
echo -e "\e[36mChecking for Conda...\e[0m"
if ! command -v conda &> /dev/null; then
    echo -e "\e[31mConda is not installed. Please install Conda before running this script.\e[0m"
    exit 1
fi

# Update Conda
echo -e "\e[36mUpdating Conda...\e[0m"
conda update -n base -c defaults conda -y

# Update pip in the base environment
echo -e "\e[36mUpdating pip...\e[0m"
conda run -n base pip install --upgrade pip

# Create environments
echo -e "\e[36mCreating custom Conda environment 'fundis_env' with Python 3.8...\e[0m"
conda create -n fundis_env python=3.8 -y

echo -e "\e[36mCreating custom Conda environment 'medaka' with Python 3.10.8...\e[0m"
conda create -n medaka python=3.10.8 -y

# Check and install Mamba if needed
echo -e "\e[36mChecking for Mamba...\e[0m"
if ! command -v mamba &> /dev/null; then
    echo -e "\e[36mMamba is not installed. Attempting to install Mamba...\e[0m"
    wget --no-check-certificate https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
    bash Mambaforge-Linux-x86_64.sh -b -p $HOME/mambaforge
    source $HOME/mambaforge/etc/profile.d/conda.sh
    conda init
    rm Mambaforge-Linux-x86_64.sh
fi

# Install packages using Mamba in the base and specific environments
echo -e "\e[36mInstalling base packages with Mamba...\e[0m"
mamba install -y -c bioconda -c conda-forge -c agbiome -c prkrekel numpy==1.26.3 pandas==2.2.0 matplotlib==3.8.2 seaborn==0.13.1 pyyaml==6.0.1 statsmodels==0.14.1

echo -e "\e[36mInstalling packages in fundis_env...\e[0m"
conda activate fundis_env
mamba install -y -c bioconda -c conda-forge -c agbiome -c prkrekel pandas==2.0.3 openblas==0.3.3 biopython==1.81 samtools==1.18 minimap2==2.26 bcftools==1.17 bwa==0.7.17 whatshap==2.1 spoa==4.1.3 racon==1.5.0 psutil==5.9.8 blast==2.15.0 pyvcf==0.6.8 fastqc==0.12.1 chopper==0.7.0
pip install NGSpeciesID
conda deactivate

echo -e "\e[36mInstalling packages in Medaka environment...\e[0m"
conda activate medaka
mamba install -y -c bioconda -c conda-forge -c agbiome -c prkrekel medaka==1.11.1
conda deactivate

# About Section
echo ""
echo "FunDiS ONT Demultiplexing Pipeline with Haplotype Phasing GUI set-up complete"
echo "@author: ian.michael.bollinger@gmail.com/researchconsultants@critical.consulting"
echo "This is a modified protocol by Stephen Douglas Russell and paid for by the Fungal Diversity Survey - FunDiS"
echo "Protocol Link: https://www.protocols.io/view/primary-data-analysis-basecalling-demultiplexing-a-dm6gpbm88lzp/v3?step=3"
