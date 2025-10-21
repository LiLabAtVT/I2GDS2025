#!/bin/bash
#SBATCH --job-name=trim_galore_G5
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=trim_galore.%j.out
#SBATCH --error=trim_galore.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hibyeon@vt.edu

# ============================================================

# 1. Activate environment

source /apps/common/software/Miniconda3/24.7.1-0/etc/profile.d/conda.sh

conda activate /home/hibyeon/.conda/envs/trimgalore_env

# ------------------------------

# 2. Define directories

# ------------------------------

RAW_DIR="/projects/intro2gds/I2GDS2025/G5_MG_AMR/rawdata"

OUT_DIR="/projects/intro2gds/I2GDS2025/G5_MG_AMR/Trim_Galore"



mkdir -p "$OUT_DIR"



# ------------------------------

# 3. Run Trim Galore for each sample

# ------------------------------

cd "$RAW_DIR"


for R1 in *_1.fastq

do

    SAMPLE=${R1%_1.fastq}

    R2=${SAMPLE}_2.fastq



    echo "Processing sample: $SAMPLE"



    trim_galore --paired --cores 8 --quality 20 --length 20 --fastqc --fastqc_args "--outdir $OUT_DIR" --output_dir "$OUT_DIR" "$RAW_DIR/$R1" "$RAW_DIR/$R2"



done




# ------------------------------

# 4. Done

# ------------------------------

echo "All samples processed successfully."
