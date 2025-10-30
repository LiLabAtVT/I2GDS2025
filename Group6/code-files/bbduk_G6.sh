#!/bin/bash

#SBATCH --account=prudenlab
#SBATCH --partition=normal_q
#SBATCH --mem=200G
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=dceglio@vt.edu

module load Miniconda3

source activate bbduk

cd /projects/intro2gds/I2GDS2025/G6_AMR_ARG/Daniel/Subset/trimmomatic_output

sample=$(ls *_[12]_paired.fastq.gz | awk -F/ '{gsub(/_[12]_paired.fastq.gz/, "", $NF); print $NF}' | sort | uniq)

for one_sample in $sample; do
    echo ${one_sample}

    bbduk.sh preallocate=t -Xmx180G \
        minlength=51 \
        in="${one_sample}"_1_paired.fastq.gz \
        in2="${one_sample}"_2_paired.fastq.gz \
        out=/projects/intro2gds/I2GDS2025/G6_AMR_ARG/Daniel/Subset/bbduk_output/"${one_sample}"_1_decontam.fastq.gz \
        out2=/projects/intro2gds/I2GDS2025/G6_AMR_ARG/Daniel/Subset/bbduk_output/"${one_sample}"_2_decontam.fastq.gz \
        stats=/projects/intro2gds/I2GDS2025/G6_AMR_ARG/Daniel/Subset/bbduk_output/"${one_sample}"_decontamination_stats.txt \
        ref='/projects/apruden_lab/thomasbyrne/metagenomics_workshop/merged_ref_5081444234059102403.fa.gz'


done
