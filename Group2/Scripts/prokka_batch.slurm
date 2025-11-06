#!/bin/bash
#SBATCH --job-name=prokka_batch
#SBATCH --output=logs/prokka_%x_%j.out
#SBATCH --error=logs/prokka_%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=normal_q
#SBATCH --account=introtogds
#SBATCH --mail-type=END,FAIL

set -euo pipefail

echo "==== Prokka batch job started at $(date) ===="

# --- Load Conda environment ---
module load Miniconda3/24.7.1-0
source $(conda info --base)/etc/profile.d/conda.sh
conda activate group2_env

# --- Directories ---
ASSEMBLY_DIR="./assemblies"
OUT_BASE="./annotations"

mkdir -p "${OUT_BASE}" logs

cd "${SLURM_SUBMIT_DIR:-$PWD}"
shopt -s nullglob

# --- Loop over assemblies ---
for ASSEMBLY in ${ASSEMBLY_DIR}/*/contigs.fasta; do
    RAW_SAMPLE=$(basename "$(dirname "${ASSEMBLY}")")
    SAMPLE="${RAW_SAMPLE/_spades/_prokka}"
    OUT_DIR="${OUT_BASE}/${SAMPLE}"

    echo "==========================================="
    echo "Starting Prokka for sample: ${SAMPLE}"
    echo "Input file: ${ASSEMBLY}"
    echo "Output dir: ${OUT_DIR}"
    echo "==========================================="

    prokka \
        --outdir "${OUT_DIR}" \
        --prefix "${SAMPLE}" \
        --cpus "${SLURM_CPUS_PER_TASK:-8}" \
        --force \
        --locustag "${SAMPLE}" \
        "${ASSEMBLY}"

    echo "Finished Prokka for ${SAMPLE}"
    echo
done

echo "==== All Prokka annotations completed successfully at $(date)! ===="
