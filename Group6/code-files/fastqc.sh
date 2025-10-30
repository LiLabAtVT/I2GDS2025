#!/bin/bash

#SBATCH --account=prudenlab  ## change to desired account 

#SBATCH --partition=normal_q   ## change to desired partition 

#SBATCH -t 3:00:00

#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --mail-user=aaamber@vt.edu  ## change to desired email address 

mkdir -p /projects/intro2gds/I2GDS2025/G6_AMR_ARG/Amber/fastqc_results

module load Miniconda3

source activate fastqc 


cd /projects/intro2gds/I2GDS2025/G6_AMR_ARG/Amber/RefPaperSubset




for file in *.fastq.gz; do

    fastqc "$file" -o /projects/intro2gds/I2GDS2025/G6_AMR_ARG/Amber/fastqc_results/



done
