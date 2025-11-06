#!/bin/bash
#SBATCH --job-name=merge_sort_all_cells
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=merge_sort_all_cells_%j.out
#SBATCH --error=merge_sort_all_cells_%j.err

INDIR="/projects/intro2gds/I2GDS2025/TestData_LinuxPeerEval/G3_TestData/04_bismark_align"
OUTDIR="/projects/intro2gds/I2GDS2025/TestData_LinuxPeerEval/G3_TestData/05_bamtools"

mkdir -p "$OUTDIR"

# Merge + Sort
for R1 in ${INDIR}/sub_R1_extracted_val_1.fq_*_bismark_bt2.bam
do
  ID=$(basename "$R1" | sed 's/.*\.fq_//; s/_bismark_bt2.bam//')
  R2="${INDIR}/sub_R2_extracted_val_2.fq_${ID}_bismark_bt2.bam"
  MERGED="${OUTDIR}/${ID}.merged.bam"
  SORTED="${OUTDIR}/${ID}_sorted.bam"
  
  # Merge R1 + R2
  bamtools merge -in "$R1" -in "$R2" -out "$MERGED"

  # Sort merged BAM
  samtools sort "$MERGED" -o "$SORTED"

  # Delete unsorted intermediate file to save space
  rm -f "$MERGED"

  echo "[$(date)] Done ${ID} -> ${SORTED}"
done
