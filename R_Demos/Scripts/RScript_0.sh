#!/bin/bash
#SBATCH -A introtogds
#SBATCH -p normal_q
#SBATCH --cpus-per-task=8
#SBATCH -J Razan_Integration_New_Nina
#SBATCH -o 	CellMarkers_DEG_Analysis_%j.out

echo "date: `date`"
echo "hostname: $HOSTNAME"; echo
echo -e "\nChecking job details for CPU_IDs..."
scontrol show job --details $SLURM_JOB_ID

#loading R
#module load site/owl-genoa/easybuild/setup
module load R/4.4.2-gfbf-2024a 


# If you want to install packages
#export R_LIBS_USER="/home/arazan/R/x86_64-pc-linux-gnu-library/4.4"



#To run your R script
Rscript CellMarkers_DEG_Analysis.R