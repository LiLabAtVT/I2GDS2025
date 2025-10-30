#!/bin/bash
#SBATCH --job-name=umi_whitelist_extract
#SBATCH --account=introtogds
#SBATCH --output=umi_whitelist_extract_%j.log    
#SBATCH --error=umi_whitelist_extract_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal_q

set -euo pipefail

# === Config ===
OUTDIR="01_umi_tools"
BC_PATTERN='(?P<discard_1>TAAGTAGAAGATGGTATATGAGAT){s<=4}(?P<cell_1>.{15})(?P<umi_1>.{0})'
SUBSET_READS=500000000
SET_CELL_NUMBER=1000     
ERR_CORR_THRESH=1

# === Locate input files ===
R1="sub_R1.fastq"
R2="sub_R2.fastq"

[[ -f "$R1" ]] || { echo "Missing $R1"; exit 1; }
[[ -f "$R2" ]] || { echo "Missing $R2"; exit 1; }

# === Prepare outdir ===
mkdir -p "$OUTDIR"

# === 1) whitelist ===
srun umi_tools whitelist \
  --subset-reads="${SUBSET_READS}" \
  --stdin "${R2}" \
  --error-correct-threshold="${ERR_CORR_THRESH}" \
  --extract-method=regex \
  --bc-pattern="${BC_PATTERN}" \
  --set-cell-number="${SET_CELL_NUMBER}" \
  --log2stderr \
  > "${OUTDIR}/whitelist_1bp_${SET_CELL_NUMBER}_allreads.txt" \
  2> "${OUTDIR}/whitelist.log"

# === 2) extract reads ===
srun umi_tools extract \
  --extract-method=regex \
  --bc-pattern="${BC_PATTERN}" \
  --error-correct-cell 1 \
  --stdin "${R2}" \
  --stdout "${OUTDIR}/sub_R2_extracted.fastq" \
  --read2-in "${R1}" \
  --read2-out="${OUTDIR}/sub_R1_extracted.fastq" \
  --filter-cell-barcode \
  --whitelist "${OUTDIR}/whitelist_1bp_${SET_CELL_NUMBER}_allreads.txt" \
  2> "${OUTDIR}/extract.log"

echo "Done. Outputs written to ${OUTDIR}/"
