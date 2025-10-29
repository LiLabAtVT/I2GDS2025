#!/bin/bash
#SBATCH -J Denoise
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --time=0-20:00:00
#SBATCH --mem=40G
#SBATCH --mail-user=yourusername@vt.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2025.7

qiime dada2 denoise-paired \
--i-demultiplexed-seqs /path/to/your/directory/data.qza \
--p-trunc-len-f 0 \
--p-trunc-len-r 200 \
--o-representative-sequences /path/to/your/directory/sequences.qza \
--o-table /path/to/your/directory/table.qza \
--o-denoising-stats /path/to/your/directory/denoising_stats.qza
