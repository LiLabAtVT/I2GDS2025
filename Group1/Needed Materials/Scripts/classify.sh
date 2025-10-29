#!/bin/bash
#SBATCH -J Classifying
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --time=0-20:00:00
#SBATCH --mem=20G
#SBATCH --mail-user=yourusername@vt.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2025.7

# Step 1: Cluster features de novo
qiime vsearch cluster-features-de-novo \
--i-sequences /path/to/your/directory/sequences.qza \
--i-table /path/to/your/directory/table.qza \
--p-perc-identity 0.97 \
--o-clustered-table /path/to/your/directory/table_97.qza \
--o-clustered-sequences /path/to/your/directory/rep_seqs_97.qza

# Step 2: Classify sequences
qiime feature-classifier classify-sklearn \
--i-reads /path/to/your/directory/rep_seqs_97.qza \
--i-classifier /projects/intro2gds/I2GDS2025/TestData_LinuxPeerEval/G1_testdata/silva-138-99-nb-classifier.qza \
--o-classification /path/to/your/directory/updated_taxonomy.qza \

# Step 3: Tabulate metadata (Optional)
qiime metadata tabulate \
--m-input-file /path/to/your/directory/updated_taxonomy.qza \
--o-visualization /path/to/your/directory/taxa-meta.qzv

# Step 4: Filter table by taxonomy
qiime taxa filter-table \
--i-table /path/to/your/directory/table_97.qza \
--i-taxonomy /path/to/your/directory/updated_taxonomy.qza \
--p-exclude mitochondria,chloroplast \
--o-filtered-table /path/to/your/directory/tax-class-filter-table.qza

# Step 5: Remove singletons
qiime feature-table filter-features \
--i-table /path/to/your/directory/tax-class-filter-table.qza \
--p-min-frequency 2 \
--o-filtered-table /path/to/your/directory/feature-frequency-filtered-table.qza

# Step 6: Generate bar plots
qiime taxa barplot \
--i-table /path/to/your/directory/feature-frequency-filtered-table.qza \
--i-taxonomy /path/to/your/directory/updated_taxonomy.qza \
--m-metadata-file /path/to/your/directory/bacteria_manifest.tsv \
--o-visualization /path/to/your/directory/taxa-barplot.qzv
