#!/bin/bash

#SBATCH --account=prudenlab
#SBATCH --partition=normal_q
#SBATCH --mem=16G
#SBATCH -t 30:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=dceglio@vt.edu

module load Miniconda3

source activate trimmomatic

cd /projects/intro2gds/I2GDS2025/G6_AMR_ARG/SampleSubset

OUTPUT_DIR=/projects/intro2gds/I2GDS2025/G6_AMR_ARG/Daniel/Subset/trimmomatic_output
ADAPTERS=TruSeq3-PE.fa

sample=$(ls *.fastq* | awk -F/ '{gsub(/_[12].fastq.gz/, "", $NF); print $NF}' | sort | uniq)

for samples in $sample; do
	echo ${samples}
	trimmomatic PE -Xmx16G -threads 8 -phred33 \
		 "${samples}"_1.fastq.gz "${samples}"_2.fastq.gz \
		 "${OUTPUT_DIR}"/"${samples}"_1_paired.fastq.gz "${OUTPUT_DIR}"/"${samples}"_1_unpaired.fastq.gz \
		 "${OUTPUT_DIR}"/"${samples}"_2_paired.fastq.gz "${OUTPUT_DIR}"/"${samples}"_2_unpaired.fastq.gz \
		 ILLUMINACLIP:${ADAPTERS}:2:30:10 \
 		 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

		echo "finished ${samples}"
done




