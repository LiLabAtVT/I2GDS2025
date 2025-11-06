#!/bin/bash
#SBATCH --job-name=checkm_batch
#SBATCH --output=logs/checkm_%x_%j.out
#SBATCH --error=logs/checkm_%x_%j.err
#SBATCH --account=introtogds
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --mail-type=END,FAIL

set -euo pipefail

echo "==== CheckM batch job started at $(date) ===="

# ------------------------------
# 1. Activate CheckM environment
# ------------------------------
module load Miniconda3/24.7.1-0
source $(conda info --base)/etc/profile.d/conda.sh
conda activate group2_env

# ------------------------------
# 2. Set directories
# ------------------------------
SPADES_DIR="assemblies"          # location of SPAdes outputs
CHECKM_OUT="checkm_results"      # where to save CheckM outputs
SUMMARY_FILE="${CHECKM_OUT}/checkm_summary.csv"

mkdir -p "${CHECKM_OUT}" logs

# ------------------------------
# 3. Loop through each assembly
# ------------------------------
for strain_dir in "${SPADES_DIR}"/*_spades; do
    if [ -d "${strain_dir}" ]; then
        sample=$(basename "${strain_dir}" _spades)
        fasta="${strain_dir}/contigs.fasta"

        if [ ! -f "${fasta}" ]; then
            echo "No contigs.fasta found for ${sample}, skipping..."
            continue
        fi

        OUTDIR="${CHECKM_OUT}/${sample}"
        mkdir -p "${OUTDIR}"

        echo "Running CheckM for sample: ${sample} ..."
        checkm lineage_wf \
            -x fasta \
            --reduced_tree \
            "${strain_dir}" \
            "${OUTDIR}" \
            --threads ${SLURM_CPUS_PER_TASK}
    fi
done

# ------------------------------
# 4. Summarize results
# ------------------------------
echo "Generating summary file ..."

checkm qa \
    -o 2 \
    -f "${SUMMARY_FILE}" \
    --tab_table \
    "${CHECKM_OUT}"/*/lineage.ms \
    "${CHECKM_OUT}"

echo "All CheckM analyses completed at $(date)"
echo "Summary saved to: ${SUMMARY_FILE}"
