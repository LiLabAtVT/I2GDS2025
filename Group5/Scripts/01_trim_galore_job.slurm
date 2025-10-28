#!/bin/bash
#SBATCH --job-name=trim_galore_G5 
#SBATCH --account=introtogds #Make sure this is an account you have access to
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=trim_galore.%j.out
#SBATCH --error=trim_galore.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=userid@vt.edu #insert your email here

# ============================================================

# 1. Load Module

module load Trim_Galore/0.6.10-GCCcore-12.3.0

# ------------------------------

# 2. Define directories

# ------------------------------

RAW_DIR= "/projects/intro2gds/I2GDS2025/TestData_LinuxPeerEval/G5_testdata"

OUT_DIR= "insert/your/trimgalore/output/filepath/here"  #make sure you change this to your preferred directory

# ------------------------------

# 3. Run Trim Galore for each sample

# ------------------------------

cd "$RAW_DIR"


for R1 in *_1.fastq

do

    SAMPLE=${R1%_1.fastq}

    R2=${SAMPLE}_2.fastq



    echo "Processing sample: $SAMPLE"



    trim_galore --paired --cores 8 --quality 30 --length 60 --fastqc --fastqc_args "--outdir $OUT_DIR" --output_dir "$OUT_DIR" "$RAW_DIR/$R1" "$RAW_DIR/$R2"



done




# ------------------------------

# 4. Done

# ------------------------------

echo "All samples processed successfully."