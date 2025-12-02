#!/bin/bash
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yourPID@vt.edu
#SBATCH --cpus-per-task=16

module load Miniconda3

source activate diamond_and_vsearch


RPOB_DB='/projects/ciwars/databases/RpoB.dmnd' #replace with location of your downloaded RpoB.dmnd file


cd #insert path to directory where all vsearch output files are located

samples=$(ls *_clean_merged.fastq | awk '{gsub(/_clean_merged.fastq/,"",$0); print}' | sort | uniq)

for sample in $samples
do
	diamond blastx -e 1e-10 --id 40 -k 1 --threads 10 -d $RPOB_DB -q ${sample}_clean_merged.fastq -o /projects/intro2gds/I2GDS2025/G6_AMR_ARG/Daniel/Subset/diamond_output/${sample}_rpob.csv --outfmt 6
done

