#!/bin/bash
#SBATCH -J QA/QC
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --time=0-20:00:00
#SBATCH --mem=20G
#SBATCH --mail-user=peterfs@vt.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2025.7

qiime dada2 denoise-paired \
--i-demultiplexed-seqs /home/peterfs/practice/QIIME/data_test.qza \
--p-trunc-len-f 0 \
--p-trunc-len-r 200 \
--o-representative-sequences /home/peterfs/practice/QIIME/sequences.qza \
--o-table /home/peterfs/practice/QIIME/table.qza \
--o-denoising-stats /home/peterfs/practice/QIIME/denoising_stats.qza
