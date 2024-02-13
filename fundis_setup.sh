#!/bin/bash
# FunDiS Nanopore Barcoding Pipeline Installation Bash Script

# Function to check if a command was successful
check_success() {
    if [ $? -ne 0 ]; then
        echo -e "\e[31mERROR: $1 failed.\e[0m"
        exit 1
    fi
}

# Determine OS Type
OS_TYPE=$(uname -s 2>/dev/null || echo "This pipeline was not built for Windows, please use WSL")

# Checking and Updating Conda
echo -e "\e[36mChecking for Conda...\e[0m"
if ! command -v conda &> /dev/null; then
    echo -e "\e[31mConda is not installed. Please install Conda before running this script.\e[0m"
    exit 1
else
    echo -e "\e[36mUpdating Conda...\e[0m"
    conda update -n base -c defaults conda -y
    check_success "Conda update"
fi

# Checking for Mamba
echo -e "\e[36mChecking for Mamba...\e[0m"
if ! command -v mamba &> /dev/null; then
    echo -e "\e[36mMamba is not installed. Attempting to install Mamba...\e[0m"
    if [ "$OS_TYPE" != "Windows" ]; then
        wget --no-check-certificate https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
        bash Mambaforge-Linux-x86_64.sh
        rm Mambaforge-Linux-x86_64.sh
    else
        echo "Please manually install Mambaforge for Windows."
    fi
fi

# Update pip in the base environment
echo -e "\e[36mUpdating pip...\e[0m"
conda run -n base pip install --upgrade pip

# Create a custom Conda environment with Python 3.8 for main FunDiS Nanopore Barcoding Pipeline GUI
echo -e "\e[36mCreating custom Conda environment 'fundis_env'...\e[0m"
conda create -n fundis_env python=3.8 -y

# Create a custom Conda environment with Python 3.10.8 for Medaka branch
echo -e "\e[36mCreating custom Conda environment 'medaka' with Python 3.10.8...\e[0m"
conda create -n medaka python=3.10.8 -y


# Give Final Instructions
echo -e "\e[32mFunDiS Nanopore Barcoding Pipeline set-up done!\e[0m"
echo -e "\e[33mFinal Instructions - Close this terminal; open a new one and run the following commands:\e[0m"

# base installs
echo -e "\e[36m'mamba install -y -c bioconda -c conda-forge -c agbiome -c prkrekel numpy==1.26.3 pandas==2.2.0 matplotlib==3.8.2 seaborn==0.13.1 pyyaml==6.0.1 statsmodels==0.14.1\e[0m"

# fundis_env installs
echo -e "\e[36m'conda activate fundis_env && mamba install -y -c bioconda -c conda-forge -c agbiome -c prkrekel pandas==2.0.3 openblas==0.3.3 biopython==1.81 samtools==1.18 minimap2==2.26 bcftools==1.17 bwa==0.7.17 whatshap==2.1 spoa==4.1.3 racon==1.5.0 psutil==5.9.8 blast==2.15.0 pyvcf==0.6.8 fastqc==0.12.1 && pip install NGSpeciesID && conda deactivate'\e[0m"

# medaka env installs
echo -e "\e[36m'conda activate medaka && mamba install -y -c bioconda -c conda-forge -c agbiome -c prkrekel medaka==1.11.1 && conda deactivate'\e[0m"

echo -e "\e[32mOnce the above are done Set-up is Complete!\e[0m"
echo -e "\e[36mRun 'conda activate fundis_env && python FunDiS_GUI.py' to execute the program!\e[0m"

# About Section
echo "@author: ian.michael.bollinger@gmail.com/researchconsultants@critical.consulting"
echo "This is a modified protocol by Stephen Douglas Russell and paid for by the Fungal Diversity Survey - FunDiS"
echo "Protocol Link: https://www.protocols.io/view/primary-data-analysis-basecalling-demultiplexing-a-dm6gpbm88lzp/v3?step=3"
