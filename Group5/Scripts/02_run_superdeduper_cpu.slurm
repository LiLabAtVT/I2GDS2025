#!/bin/bash                          
#SBATCH --job-name=sdeduper_cpu       
#SBATCH --account=introtogds #Make sure this is an account you have access to      
#SBATCH --partition=normal_q          
#SBATCH --time=06:00:00               
#SBATCH --cpus-per-task=8             
#SBATCH --mem=32G                        
#SBATCH --output=sdeduper.%j.out  
#SBATCH --error=sduduper.%j.err  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=userid@vt.edu #insert your email here 

shopt -s nullglob                     # If no files match a glob, expand to empty (not the pattern)

# ---- Paths ----
INDIR="insert/your/input/directory/filepath/here" #make sure you change this to your output directory from TrimGalore     
OUTDIR="insert/your/output/directory/filepath/here" #make sure you change this to your preferred directory 

# ---- Activate software environment ----
source /projects/intro2gds/I2GDS2025/tools/miniconda3/etc/profile.d/conda.sh  # Load conda functions
conda activate htstream12            # Activate env that provides hts_SuperDeduper

# ---- Collect all R1 files (Trim Galore naming: *_1_val_1.fq) ----
R1S=( "$INDIR"/*_1_val_1.fq )        # Array of R1 fastq files to iterate over
[[ ${#R1S[@]} -eq 0 ]] && echo "[INFO] No *_1_val_1.fq under $INDIR; nothing to do."  # Inform when empty set

# ---- Main loop over samples ----
for R1 in "${R1S[@]}"; do            # Iterate through every R1 file found
  base=$(basename "$R1")             # Get filename only (strip directory)
  sample=${base%_1_val_1.fq}         # Strip suffix to get sample ID (e.g., SRRxxxx)

  R2="$INDIR/${sample}_2_val_2.fq"   # Expected R2 mate based on Trim Galore naming
  if [[ ! -f "$R2" ]]; then          # If paired R2 is missing,
    echo "[WARN] Missing R2 for $sample; skip."  # warn and skip this sample
    continue
  fi

  prefix="$OUTDIR/${sample}"         # Output prefix used by HTStream tools
  echo "[INFO] $sample -> hts_SuperDeduper"   # Progress message

  # Run deduplication; -1/-2: paired inputs; -f: output prefix; -F: overwrite if exists
  hts_SuperDeduper -1 "$R1" -2 "$R2" -f "$prefix" -F \
    || { echo "[WARN] hts_SuperDeduper failed for $sample; continue."; continue; }

  # Quick sanity check that expected gzipped outputs exist and are non-empty
  if [[ -s "${prefix}_R1.fastq.gz" && -s "${prefix}_R2.fastq.gz" ]]; then
    echo "[OK] $sample done."
  else
    echo "[WARN] Outputs missing/empty for $sample."
  fi
done

echo "[DONE] Pipeline finished."      # Final message when loop completes
