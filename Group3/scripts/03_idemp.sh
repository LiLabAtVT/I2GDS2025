#!/bin/bash
#SBATCH --job-name=idemp_split
#SBATCH --account=introtogds
#SBATCH --output=03_idemp/idemp_split_%j.log
#SBATCH --error=03_idemp/idemp_split_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal_q

OUTDIR="03_idemp"
mkdir -p "$OUTDIR"
WL="01_umi_tools/whitelist_1bp_1000_allreads.txt"
R1="02_Trim_Galore/sub_R1_extracted_val_1.fq"
R2="02_Trim_Galore/sub_R2_extracted_val_2.fq"

# 1)  barcode_r.txt   
awk '{print $1, NR}' "$WL" > "$OUTDIR/barcode_r.txt"

# 2)  generate index.fastq
awk 'NR%4==1{
  if (match($0, /_([ACGTN]+)_/, m)) {
    bc=m[1];
  } else {
    bc="N";
  }
  q=""; for(i=1;i<=length(bc);i++) q=q"F";
  print "@idx" ++c; print bc; print "+"; print q;
}' "$R1" > "$OUTDIR/index.fastq"

# 3) run idemp
idemp \
  -b "$OUTDIR/barcode_r.txt" \
  -I1 "$OUTDIR/index.fastq" \
  -R1 "$R1" \
  -R2 "$R2" \
  -m 0 \
  -o "$OUTDIR"
