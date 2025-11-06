#!/bin/bash
#SBATCH --job-name=picard_dedup
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=picard_dedup_%j.out
#SBATCH --error=picard_dedup_%j.err

INDIR="/projects/intro2gds/I2GDS2025/TestData_LinuxPeerEval/G3_TestData/05_bamtools"
OUTDIR="/projects/intro2gds/I2GDS2025/TestData_LinuxPeerEval/G3_TestData/06_picard"

mkdir -p "$OUTDIR"

for BAM in ${INDIR}/*_sorted.bam
do
  ID=$(basename "$BAM" _sorted.bam)
  OUT_BAM="${OUTDIR}/${ID}_dedup.bam"
  METRICS="${OUTDIR}/${ID}_metrics.txt"

  echo "[$(date)] Processing ${ID} ..."

  picard MarkDuplicates \
    I="$BAM" \
    O="$OUT_BAM" \
    M="$METRICS" \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=SILENT

  echo "[$(date)] Done ${ID} -> ${OUT_BAM}"
done
