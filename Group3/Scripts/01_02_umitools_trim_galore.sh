#!/bin/bash
#SBATCH --job-name=extract_trim
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=extract_trim_%j.out
#SBATCH --error=extract_trim_%j.err


THREADS="${SLURM_CPUS_PER_TASK:-8}"
OUTDIR="extract_trim_out"
mkdir -p "$OUTDIR"

sub_R1="sub_R1.fastq"
sub_R2="sub_R2.fastq"

# -------------------------------
# Step 1: whitelist
# -------------------------------
umi_tools whitelist \
  --stdin "$sub_R2" \
  --extract-method=regex \
  --bc-pattern='(?P<discard_1>TAAGTAGAAGATGGTATATGAGAT){s<=4}(?P<cell_1>.{15})(?P<umi_1>.{0})' \
  --subset-reads=10000000 \
  --error-correct-threshold=1 \
  --set-cell-number=100 \
  --log2stderr > "${OUTDIR}/whitelist_top100.txt"

# -------------------------------
# Step 2: extract
# -------------------------------
umi_tools extract \
  --extract-method=regex \
  --bc-pattern='(?P<discard_1>TAAGTAGAAGATGGTATATGAGAT){s<=4}(?P<cell_1>.{15})(?P<umi_1>.{0})' \
  --error-correct-cell \
  --stdin "$sub_R2" --stdout "${OUTDIR}/${sub_R2%.fastq}_extracted.fastq" \
  --read2-in "$sub_R1" --read2-out "${OUTDIR}/${sub_R1%.fastq}_extracted.fastq" \
  --filter-cell-barcode --whitelist="${OUTDIR}/whitelist_top100.txt"

# -------------------------------
# Step 3: Trim Galore
# -------------------------------
trim_galore --paired --quality 20 --length 20 --cores "$THREADS" \
  --output_dir "$OUTDIR" \
  "${OUTDIR}/${sub_R1%.fastq}_extracted.fastq" \
  "${OUTDIR}/${sub_R2%.fastq}_extracted.fastq"
