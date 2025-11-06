#!/bin/bash
#SBATCH --job-name=bismark_align
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=bismark_align_%j.out
#SBATCH --error=bismark_align_%j.err


GENOME_DIR="/projects/intro2gds/I2GDS2025/G3_SingleCell/Rui/reference/Bisulfite_Genome_human"
INDIR="03_idemp"
OUTDIR="04_bismark_align"
THREADS="${SLURM_CPUS_PER_TASK:-8}"


mkdir -p "$OUTDIR"


# ---- Align R1 reads (PBAT mode) ----
for f in ${INDIR}/*_R1_*.fastq.gz; do
  base=$(basename "$f" .fastq.gz)
  echo "[$(date)] Aligning (PBAT) $base ..."
  bismark --genome "$GENOME_DIR" --pbat -p "$THREADS" \
          -o "$OUTDIR" "$f"
done

# ---- Align R2 reads (default mode) ----
for f in ${INDIR}/*_R2_*.fastq.gz; do
  base=$(basename "$f" .fastq.gz)
  echo "[$(date)] Aligning (default) $base ..."
  bismark --genome "$GENOME_DIR" -p "$THREADS" \
          -o "$OUTDIR" "$f"
done
