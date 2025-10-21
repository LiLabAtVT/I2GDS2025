#!/bin/bash
#SBATCH --job-name=Fastqc_G5
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --output=fastqc.%j.out
#SBATCH --error=fastqc.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aegreen@vt.edu

# ============================================================

# Load Module

module reset
module load FastQC/0.12.1-Java-11

# Change directory

cd /projects/intro2gds/I2GDS2025/G5_MG_AMR/03_BWA

# Run fastqc on all files

fastqc -o /projects/intro2gds/I2GDS2025/G5_MG_AMR/04_fastqc *.fq.gz

