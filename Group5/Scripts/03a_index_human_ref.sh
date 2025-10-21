#!/bin/bash
#SBATCH --job-name=index_hg_G5
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=BWA.%j.out
#SBATCH --error=BWA.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aegreen@vt.edu

# ============================================================

# 1. Load Module
module load BWA/0.7.18-GCCcore-13.3.0

# 2. Change directory
cd /projects/intro2gds/I2GDS2025/G5_MG_AMR/03a_Human_ref

# 3. Index the human reference genome for BWA

bwa index hg19.fa

