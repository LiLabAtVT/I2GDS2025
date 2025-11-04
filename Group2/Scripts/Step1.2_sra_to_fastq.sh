#!/bin/bash
# -------------------------------------------
# sra_to_fastq.sh
#SBATCH --account=introtogds
#SBATCH --output=/projects/intro2gds/I2GDS2025/G2_PlantDisease/Jingjing/script/logs/sra_to_fastq.%j.out
#SBATCH --error=/projects/intro2gds/I2GDS2025/G2_PlantDisease/Jingjing/script/logs/sra_to_fastq.%j.err
#SBATCH --job-name=sra_to_fastq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mail-user=jingjingy@vt.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=200GB
#SBATCH --cpus-per-task=8
# -------------------------------------------

SRA_FILE="/projects/intro2gds/I2GDS2025/G2_PlantDisease/Jingjing/RawData/SRR29108932"
FASTQ_DIR="/projects/intro2gds/I2GDS2025/G2_PlantDisease/Jingjing/RawData/fastq_historical"
mkdir -p "$FASTQ_DIR"
THREADS=8

echo "Converting $SRA_FILE to FASTQ ..."

# 1. Add SRA Toolkit to PATH
export PATH="/projects/intro2gds/I2GDS2025/G2_PlantDisease/Jingjing/sratoolkit.3.2.1-ubuntu64/bin:$PATH"

# 2. Convert SRA to FASTQ
fasterq-dump "$SRA_FILE" \
  --split-files \
  --threads "$THREADS" \
  -O "$FASTQ_DIR"

# 3. Compress FASTQ
echo "Compressing FASTQ files ..."
gzip -f "$FASTQ_DIR"/*.fastq

echo "Done! Compressed FASTQ files are in $FASTQ_DIR"