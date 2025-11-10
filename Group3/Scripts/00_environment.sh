#!/bin/bash
#SBATCH --job-name=install_bisulfite_env
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=install_bisulfite_%j.out
#SBATCH --error=install_bisulfite_%j.err


echo "[$(date)] Creating conda environment 'bisulfite'..."
conda create -n bisulfite -c conda-forge -c bioconda python=3.9 -y

# Activate environment
eval "$(conda shell.bash hook)"
conda activate bisulfite

# UMI-tools
echo "Installing UMI-tools..."
python -m pip install --no-cache-dir "umi_tools==1.1.3"

# Trim Galore
echo "Installing Trim Galore..."
conda config --set channel_priority flexible
conda install -c bioconda trim-galore -y
conda config --set channel_priority strict


# Demultiplexing
echo "Installing idemp..."
conda install -c conda-forge zlib -y
conda install -c conda-forge make gcc gxx -y
git clone https://github.com/yhwu/idemp.git
cd idemp
make
cp idemp $CONDA_PREFIX/bin/


# Alignment and BAM tools
echo "Installing Bismark, Bowtie2, Samtools..."
conda install -c bioconda -c conda-forge bismark bowtie2 samtools -y
conda install -c conda-forge ncurses=6.4 -y
ln -s $CONDA_PREFIX/lib/libncurses.so.6 $CONDA_PREFIX/lib/libncurses.so.5


echo "Installing BAM/bed utilities..."
conda install -c bioconda bamtools bedtools -y

# Duplicate marking
echo "Installing Picard..."
conda install -c bioconda picard -y

