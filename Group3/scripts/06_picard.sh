#!/bin/bash
#SBATCH --job-name=picard_dedup
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=06_picard/picard_%A_%a.log
#SBATCH --error=06_picard/picard_%A_%a.err
#SBATCH --array=1-1000          

module load picard
module load samtools

IN_DIR="/projects/intro2gds/I2GDS2025/G3_SingleCell/BisulfiteData/SRR19391270/subset/05_bamtools"
OUT_DIR="/projects/intro2gds/I2GDS2025/G3_SingleCell/BisulfiteData/SRR19391270/subset/06_picard"

mkdir -p "$OUT_DIR"


BAM_FILE=$(ls $IN_DIR/*.merged.bam | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [ -z "$BAM_FILE" ]; then
  echo "No input file for task $SLURM_ARRAY_TASK_ID"
  exit 0
fi

BASE=$(basename "$BAM_FILE" .merged.bam)


samtools sort -@4 -o $OUT_DIR/${BASE}.sorted.bam "$BAM_FILE"


picard MarkDuplicates \
    I=$OUT_DIR/${BASE}.sorted.bam \
    O=$OUT_DIR/${BASE}.dedup.bam \
    M=$OUT_DIR/${BASE}.dedup.metrics.txt \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=SILENT


samtools index $OUT_DIR/${BASE}.dedup.bam

echo "✅ Done: $BASE"
