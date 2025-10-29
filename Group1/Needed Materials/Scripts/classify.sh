#!/bin/bash
#SBATCH -J Classifying
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

# Step 1: Cluster features de novo
qiime vsearch cluster-features-de-novo \
--i-sequences /home/peterfs/practice/QIIME/sequences.qza \
--i-table /home/peterfs/practice/QIIME/table.qza \
--p-perc-identity 0.97 \
--o-clustered-table /home/peterfs/practice/QIIME/table_97.qza \
--o-clustered-sequences /home/peterfs/practice/QIIME/rep_seqs_97.qza

# Step 2: Classify sequences
qiime feature-classifier classify-sklearn \
--i-reads /home/peterfs/practice/QIIME/rep_seqs_97.qza \
--i-classifier /home/peterfs/practice/files/silva-138-99-nb-classifier.qza \
--o-classification /home/peterfs/practice/QIIME/updated_taxonomy.qza \

# Step 3: Tabulate metadata (Optional)
qiime metadata tabulate \
--m-input-file /home/peterfs/practice/QIIME/updated_taxonomy.qza \
--o-visualization /home/peterfs/practice/QIIME/taxa-meta.qzv

# Step 4: Filter table by taxonomy
qiime taxa filter-table \
--i-table /home/peterfs/practice/QIIME/table_97.qza \
--i-taxonomy /home/peterfs/practice/QIIME/updated_taxonomy.qza \
--p-exclude mitochondria,chloroplast \
--o-filtered-table /home/peterfs/practice/QIIME/tax-class-filter-table.qza

# OPTIONAL: Remove singletons
qiime feature-table filter-features \
--i-table /home/peterfs/practice/QIIME/tax-class-filter-table.qza \
--p-min-frequency 2 \
--o-filtered-table /home/peterfs/practice/QIIME/feature-frequency-filtered-table.qza

# Step 5: Generate bar plots
qiime taxa barplot \
--i-table /home/peterfs/practice/QIIME/feature-frequency-filtered-table.qza \
--i-taxonomy /home/peterfs/practice/QIIME/updated_taxonomy.qza \
--m-metadata-file /home/peterfs/practice/files/bacteria_manifest.tsv \
--o-visualization /home/peterfs/practice/QIIME/taxa-barplot.qzv
