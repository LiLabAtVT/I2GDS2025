#!/bin/bash
#SBATCH -J QA/QC
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --time=0-08:00:00
#SBATCH --mem=2G
#SBATCH --mail-user=peterfs@vt.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2025.7

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path /home/peterfs/practice/files/bacteria_manifest.tsv \
--output-path /home/peterfs/practice/QIIME/data.qza \
--input-format PairedEndFastqManifestPhred33V2 \

qiime demux summarize \
--i-data /home/peterfs/practice/QIIME/data.qza \
--o-visualization /home/peterfs/practice/QIIME/qualityplot.qzv\