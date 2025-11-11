#!/bin/bash
#SBATCH -A introtogds
#SBATCH -p normal_q
#SBATCH --cpus-per-task=8
#SBATCH -J scRNA_Seq
#SBATCH -o 	Integration_%j.out #Replace this with the name of the job name  you want to run.

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
Rscript Integration.R #Replace this with the name of the script you want to run.
