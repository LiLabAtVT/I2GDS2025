#!/bin/bash
#SBATCH --job-name=index_hg_G5
#SBATCH --account=introtogds #Make sure this is an account you have access to
#SBATCH --partition=normal_q 
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=index_hg.%j.out
#SBATCH --error=index_hg.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=userid@vt.edu #insert your email here

# ============================================================

# 1. Load Module
module load BWA/0.7.18-GCCcore-13.3.0

# 2. Change directory
cd /insert/human/ref/genome/filepath/here #change to whatever directory you downloaded your human reference genome in

# 3. Index the human reference genome for BWA

bwa index hg19.fa

