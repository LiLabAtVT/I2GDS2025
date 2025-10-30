#!/bin/bash
#SBATCH --job-name=bismark_SE
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=04_bismark_align/bismark_%A_%a.out
#SBATCH --error=04_bismark_align/bismark_%A_%a.err
# Submit with: sbatch --array=1-1000 04_bismark_align.sh

set -euo pipefail

# -------------------------------
# Basic settings
# -------------------------------
THREADS="${SLURM_CPUS_PER_TASK:-8}"
GENOME_DIR="/projects/intro2gds/I2GDS2025/G3_SingleCell/Rui/reference/Bisulfite_Genome_human"
INDIR="03_idemp"
ALIGN_DIR="04_bismark_align"
DONE_DIR=".bmk_done"

mkdir -p "${ALIGN_DIR}" "${DONE_DIR}"

ID_FILE="${INDIR}/cell_ids.txt"
ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ID_FILE}")

R1="${INDIR}/sub_R1_extracted_val_1.fq_${ID}.fastq.gz"
R2="${INDIR}/sub_R2_extracted_val_2.fq_${ID}.fastq.gz"

[[ -f "${R1}" ]] || { echo "[ERR] Missing ${R1}"; exit 1; }
[[ -f "${R2}" ]] || { echo "[ERR] Missing ${R2}"; exit 1; }

echo "[$(date)] Start Drop-BS alignment for cell ${ID}"
echo "  Genome: ${GENOME_DIR}"
echo "  Threads: ${THREADS}"
echo "  R1 (PBAT): ${R1}"
echo "  R2 (default): ${R2}"

# -------------------------------
# (1) R2: default single-end alignment
# -------------------------------
bismark \
  --genome "${GENOME_DIR}" \
  -p "${THREADS}" \
  -o "${ALIGN_DIR}" \
  --temp_dir "${ALIGN_DIR}/tmp_${ID}_R2" \
  --gzip \
  "${R2}"

# -------------------------------
# (2) R1: PBAT single-end alignment (⚠️ no --gzip)
# -------------------------------
bismark \
  --genome "${GENOME_DIR}" \
  --pbat \
  -p "${THREADS}" \
  -o "${ALIGN_DIR}" \
  --temp_dir "${ALIGN_DIR}/tmp_${ID}_R1" \
  "${R1}"

# -------------------------------
# (3) Extract basic metrics only (mapping efficiency & CpG methylation)
#     Writes a simple TSV: 04_bismark_align/basic_metrics.tsv
# -------------------------------
R1_REPORT="${ALIGN_DIR}/sub_R1_extracted_val_1.fq_${ID}_bismark_bt2_SE_report.txt"
R2_REPORT="${ALIGN_DIR}/sub_R2_extracted_val_2.fq_${ID}_bismark_bt2_SE_report.txt"
METRICS_TSV="${ALIGN_DIR}/basic_metrics.tsv"

# Create header once
if [[ ! -s "${METRICS_TSV}" ]]; then
  echo -e "cell_id\tread_set\tmapping_efficiency_pct\tcpg_methylation_pct" > "${METRICS_TSV}"
fi

parse_and_append () {
  local report="$1"
  local label="$2"
  if [[ -f "${report}" ]]; then
    local map cpg
    map=$(awk -F': *' '/^Mapping efficiency/ {gsub("%","",$2); print $2}' "${report}" | head -n1)
    cpg=$(awk -F': *' '/^C methylated in CpG context/ {gsub("%","",$2); print $2}' "${report}" | head -n1)
    [[ -n "${map:-}" && -n "${cpg:-}" ]] && echo -e "${ID}\t${label}\t${map}\t${cpg}" >> "${METRICS_TSV}"
  fi
}

parse_and_append "${R2_REPORT}" "R2_default_SE"
parse_and_append "${R1_REPORT}" "R1_PBAT_SE"

# -------------------------------
# (4) Finish
# -------------------------------
touch "${DONE_DIR}/${ID}.done"
echo "[$(date)] Finished cell ${ID}"
