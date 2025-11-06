#!/bin/bash
#SBATCH --job-name=idemp
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=idemp_%j.out
#SBATCH --error=idemp_%j.err

OUTDIR="03_idemp"
mkdir -p "$OUTDIR"

# input
WL="extract_trim_out/whitelist_top100.txt"
R1="extract_trim_out/sub_R1_extracted_val_1.fq"
R2="extract_trim_out/sub_R2_extracted_val_2.fq"

# 1. barcode_r.txt
awk '{print $1, NR}' "$WL" > "$OUTDIR/barcode_r.txt"

# 2. index.fastq
awk 'NR%4==1{
  if (match($0, /_([ACGTN]+)_/, m)) bc=m[1]; else bc="N";
  q=""; for(i=1;i<=length(bc);i++) q=q"F";
  print "@idx" ++c; print bc; print "+"; print q;
}' "$R1" > "$OUTDIR/index.fastq"

# 3. run idemp
idemp \
  -b "$OUTDIR/barcode_r.txt" \
  -I1 "$OUTDIR/index.fastq" \
  -R1 "$R1" \
  -R2 "$R2" \
  -m 0 \
  -o "$OUTDIR"
