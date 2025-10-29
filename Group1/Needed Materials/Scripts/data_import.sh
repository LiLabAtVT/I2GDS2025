#!/bin/bash
#SBATCH -J QA/QC
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --time=0-08:00:00
#SBATCH --mem=2G
#SBATCH --mail-user=yourusername@vt.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2025.7

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path /path/to/your/directory/bacteria_manifest.tsv \
--output-path /path/to/your/directory/data.qza \
--input-format PairedEndFastqManifestPhred33V2 \

qiime demux summarize \
--i-data /path/to/your/directory/data.qza \
--o-visualization /path/to/your/directory/qualityplot.qzv\
