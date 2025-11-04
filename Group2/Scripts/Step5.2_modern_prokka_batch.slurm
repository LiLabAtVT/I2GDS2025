#!/bin/bash
#SBATCH --job-name=prokka_batch
#SBATCH --output=logs/prokka_%x_%j.out
#SBATCH --error=logs/prokka_%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=normal_q
#SBATCH --account=introtogds

set -euo pipefail

# --- Load environment ---
module load Miniconda3/24.7.1-0
source activate ~/.conda/envs/prokka_env

# --- Set paths ---
ASSEMBLY_DIR="all_bins"       # change to correct location relative to script
OUT_BASE="data/annotations"

mkdir -p logs
mkdir -p "${OUT_BASE}"

# --- Change directory to where the job was submitted ---
cd "${SLURM_SUBMIT_DIR:-$PWD}"

# Enable nullglob so loop doesn’t fail if no matches
shopt -s nullglob

# --- Loop over all FASTA assemblies ---
for ASSEMBLY in ${ASSEMBLY_DIR}/SRR*.fasta; do
    SAMPLE=$(basename "${ASSEMBLY}" .fasta)
    OUT_DIR="${OUT_BASE}/${SAMPLE}_prokka"

    echo "==========================================="
    echo "Starting Prokka for sample: ${SAMPLE}"
    echo "Input file: ${ASSEMBLY}"
    echo "Output dir: ${OUT_DIR}"
    echo "==========================================="

    prokka \
        --outdir "${OUT_DIR}" \
        --prefix "${SAMPLE}" \
        --cpus ${SLURM_CPUS_PER_TASK:-8} \
        --force \
        --locustag "${SAMPLE}" \
        "${ASSEMBLY}"

    echo "Finished Prokka for ${SAMPLE}"
    echo
done

echo  "All Prokka annotations completed successfully!"
