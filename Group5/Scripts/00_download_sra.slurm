#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH --account=introtogds
#SBATCH --output=sra_download.%j.out
#SBATCH --error=sra_download.%j.err
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hibyeon@vt.edu
# ---------------------------------
#1. Load a module
module load SRA-Toolkit

#2. Move to the working directory (if already there, this step is not necessary).

#cd /projects/intro2gds/I2GDS2025/G5_MG_AMR/raw_data

#3. Download SRA data (prefetch)

while read SRR; do
    echo "Downloading $SRR ..."
    prefetch $SRR
done < SRR_Acc_List.txt

#4. Convert .sra to fastq format (fasterq-dump)

for SRR in $(cat SRR_Acc_List.txt); do
    echo "Converting $SRR to FASTQ ..."
    fasterq-dump --split-files --threads $SLURM_CPUS_PER_TASK --outdir . $SRR
done

#5. Compress the data (tar.gz)

echo "Creating tar.gz archive of all FASTQ files ..."
tar -czvf PRJNA678454_fastq.tar.gz SRR*_*.fastq
echo "All downloads, conversions, and compression completed!"
