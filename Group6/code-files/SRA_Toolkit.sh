#!/bin/bash

#SBATCH --account=metagen   ## change to desired account 

#SBATCH --partition=normal_q  ## change to desired partition 

#SBATCH -t 12:00:00

#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --mail-user=aaamber@vt.edu  ## change to desired email address 


mkdir -p fastq_files 

# Add your SRA Toolkit bin directory to PATH

export PATH=/projects/intro2gds/I2GDS2025/G6_AMR_ARG/Amber/SRA_Demo/sratoolkit.3.2.1-alma_linux64/bin:$PATH

for acc in $(cat SRR_File.txt); do

    prefetch "$acc"

    fasterq-dump "$acc" --split-files --threads 4 -O fastq_files

done
