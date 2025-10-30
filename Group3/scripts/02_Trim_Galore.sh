#!/bin/bash
#SBATCH --job-name=Trim_Galore
#SBATCH --account=introtogds
#SBATCH --output=Trim_Galore_%j.log
#SBATCH --error=Trim_Galore_%j.err
#SBATCH --time=12:00:00          
#SBATCH --mem=32G                
#SBATCH --cpus-per-task=4        
#SBATCH --partition=normal_q    

# Creat output directory
mkdir -p 02_Trim_Galore
OUTDIR="02_Trim_Galore"

# define input directory
R1="01_umi_tools/sub_R1_extracted.fastq"
R2="01_umi_tools/sub_R2_extracted.fastq"

# Running Trim Galore
echo "Running Trim Galore! on ${R1} and ${R2} ..."
trim_galore \
  --paired \
  --clip_R1 15 \
  --clip_R2 40 \
  --cores 4 \
  --output_dir "$OUTDIR" \
  "$R1" "$R2"

echo "All trimming steps finished at $(date)"
