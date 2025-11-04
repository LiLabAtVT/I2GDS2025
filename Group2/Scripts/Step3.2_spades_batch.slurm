#!/bin/bash
#SBATCH --job-name=spades_batch
#SBATCH --output=logs/spades_%x_%j.log
#SBATCH --error=logs/spades_%x_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --account=introtogds


# Load SPAdes module
module load spades

# Stop on first error
set -e

# Directories
READS_DIR="data/Modern_strain"
OUT_BASE="Carter/data/assemblies"

# Create output + log directories
mkdir -p "${OUT_BASE}" logs

# Loop over all *_1.fastq.gz files
for R1 in ${READS_DIR}/*_1.fastq.gz; do
    SAMPLE=$(basename "${R1}" _1.fastq.gz)
    R2="${READS_DIR}/${SAMPLE}_2.fastq.gz"
    OUT_DIR="${OUT_BASE}/${SAMPLE}_spades"

    echo "======================================"
    echo "Starting SPAdes for ${SAMPLE}"
    echo "Output -> ${OUT_DIR}"
    echo "======================================"

    mkdir -p "${OUT_DIR}"

    spades.py \
        -1 "${R1}" \
        -2 "${R2}" \
        -o "${OUT_DIR}" \
        -t ${SLURM_CPUS_PER_TASK} \
        -m 32 \
        --isolate \
        --cov-cutoff auto

    echo "SPAdes assembly finished for ${SAMPLE}"
    echo "Main contigs file: ${OUT_DIR}/contigs.fasta"
done

echo "======================================"
echo "All assemblies completed!"
echo "======================================"
