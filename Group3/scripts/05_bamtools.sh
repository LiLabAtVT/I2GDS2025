#!/bin/bash
#SBATCH --job-name=bamtools_merge
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=05_bamtools/bamtools_%A_%a.out
#SBATCH --error=05_bamtools/bamtools_%A_%a.err
# Submit with: sbatch --array=1-1000 05_bamtools.sh

set -euo pipefail

# -------------------------------
# Basic paths
# -------------------------------
INDIR="04_bismark_align"
OUTDIR="05_bamtools"
ID_FILE="03_idemp/cell_ids.txt"

mkdir -p "${OUTDIR}"

ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ID_FILE}")

R1_BAM="${INDIR}/sub_R1_extracted_val_1.fq_${ID}_bismark_bt2.bam"
R2_BAM="${INDIR}/sub_R2_extracted_val_2.fq_${ID}_bismark_bt2.bam"
MERGED_BAM="${OUTDIR}/${ID}.merged.bam"

# -------------------------------
# Check input existence
# -------------------------------
if [[ ! -f "${R1_BAM}" || ! -f "${R2_BAM}" ]]; then
  echo "[ERR] Missing BAM(s) for cell ${ID}"
  exit 1
fi

echo "[$(date)] Merging BAMs for cell ${ID}"
echo "R1: ${R1_BAM}"
echo "R2: ${R2_BAM}"

# -------------------------------
# Merge with bamtools
# -------------------------------
bamtools merge \
  -in "${R1_BAM}" \
  -in "${R2_BAM}" \
  -out "${MERGED_BAM}"

echo "[$(date)] Finished merge for ${ID} -> ${MERGED_BAM}"
