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
# 4. Generate single summary table
# ------------------------------
echo "Generating summary file ..."
echo "Sample,Completeness,Contamination" > "${SUMMARY_FILE}"

for OUTDIR in ${CHECKM_OUT}/*; do
    if [ -f "${OUTDIR}/storage/bin_stats_ext.tsv" ]; then
        sample=$(basename "${OUTDIR}")
        # Extract completeness and contamination from CheckM output
        awk -v s="${sample}" -F'\t' '/completeness/ {
            gsub(/[^0-9.]/, "", $2); comp=$2
        } /contamination/ {
            gsub(/[^0-9.]/, "", $2); cont=$2
        } END {
            if (comp != "" && cont != "") print s","comp","cont
        }' "${OUTDIR}/storage/bin_stats_ext.tsv" >> "${SUMMARY_FILE}"
    fi
done

echo "All CheckM analyses completed at $(date)"
echo "Summary saved to: ${SUMMARY_FILE}"
