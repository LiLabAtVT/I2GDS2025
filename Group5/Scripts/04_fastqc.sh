#!/bin/bash
#SBATCH --job-name=Fastqc_G5
#SBATCH --account=introtogds #Make sure this is an account you have access to
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --output=fastqc.%j.out
#SBATCH --error=fastqc.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=userid@vt.edu #insert your email here

# ============================================================

# Load Module

module reset
module load FastQC/0.12.1-Java-11

# Change directory

cd /insert/your/input/directory/here #make sure you change this to your output directory from BWA

# Run fastqc on all files

#change filepath below to your preferred output directory
fastqc -o /insert/your/output/directory/here *.fq.gz

