#!/bin/bash
#SBATCH --job-name=BWA_G5-test
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

# ------------------------------
# 2. Define directories
# ------------------------------

IN_DIR="/projects/intro2gds/I2GDS2025/G5_MG_AMR/02_after_deduper"

OUT_DIR="/projects/intro2gds/I2GDS2025/G5_MG_AMR/03_BWA"

# ------------------------------
# 3. Run BWA for each sample
# ------------------------------

cd "$IN_DIR"

for R1 in *_R1.fastq.gz
do

    SAMPLE=${R1%_R1.fastq.gz}
    R2=${SAMPLE}_R2.fastq.gz

    echo "Processing sample: $SAMPLE"

     bwa mem /projects/intro2gds/I2GDS2025/G5_MG_AMR/03a_Human_ref/hg19.fa "$IN_DIR/$R1" "$IN_DIR/$R2" | \
	     samtools fastq -t -T BX -f 4 \
	     -1 "$OUT_DIR/${SAMPLE}_1.fq.gz" \
	     -2 "$OUT_DIR/${SAMPLE}_2.fq.gz" \
	     -s "$OUT_DIR/${SAMPLE}_singletons.fq.gz"  
done
