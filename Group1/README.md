# Introduction to Group 1 16S Amplicon Pipeline
Our article examines the relative abundance, taxonomic profiles, and community structure of bacterial and fungal communities associated with parsley (Petroselinum crispum) and celery (Apium graveolens) roots via monocropping and intercropping systems.  The study aims to provide a baseline understanding of how intercropping influences rhizosphere microbial dynamics.

Article: https://bmcgenomdata.biomedcentral.com/articles/10.1186/s12863-025-01351-0#additional-information

This pipeline utilizes QIIME2's (Quantitative Insights Into Microbial Ecology 2) "Amplicon" Distribution, a suite of plug-ins that provide broad analytic functionality to support microbiome marker gene analysis from raw sequencing data to pubication-quality visualizations and statistics.
## Installing QIIME2 into Conda environment
Installing the base distribution's conda environment: QIIME 2 recommends creating a new environment specifically for the QIIME 2 distribution, as they are many dependencies that you wouldn't want in an existing environment
```{r}
# Installation for Linux/Windows users
conda env create \
  --name qiime2-amplicon-2025.7 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.7/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml
```

```{r}
# Installation for MacOS users (Apple Silicon)
CONDA_SUBDIR=osx-64 conda env create \
  --name qiime2-amplicon-2025.7 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.7/amplicon/released/qiime2-amplicon-macos-latest-conda.yml
conda activate qiime2-amplicon-2025.7
conda config --env --set subdir osx-64
```
Verify your installation:
```{r}
# Test your installation
conda deactivate
conda activate qiime2-amplicon-2025.7
qiime info
```

Go ahead and activate your QIIME version (after verification)
```{linux}
conda activate qiime2-2025.7 # reference your specific version of QIIME
OR
source activate qiime2-amplicon-2025.7
```

## Test Data
The test data is labeled "Group1_TestData" in the following ARC Directory: /projects/intro2gds/I2GDS2025/TestData_LinuxPeerEval
It consists of **12 bacterial samples** (forward & reverse reads) in FASTQ format (.fastq.gz), which have been deposited in the Sequence Read Archive (SRA) of the National Center for Biotechnology Information (NCBI) under the Bioproject Accession numbers; **SRP540554 (16S rRNA)** and **SRP540675 (ITS)**.

We downloaded the provided script to retrieve each samples' reads: ena-file-download-read_run-SRP540554-fastq_ftp-20250915-1609.sh
This tutorial will focus on the bacterial samples (Link: http)://identifiers.org/insdc.sra:SRP540554)

Here are the provided materials for this tutorial before we move to QIIME2:
- Manifest file ('Needed Materials' folder)
- Test data ('Needed Materials' folder)
- Classifier (silva-138-99-nb-classifier.qza)
- Bash scripts ('Needed Materials' folder)
# QIIME2 Pipeline
The goal is to generate high-resolution OTUs (operational taxonomic unit). OTUs are employed to classify groups of similar sequences at 97% similarity.
We have provided our manifest file (TSV format) into our test data; it will be used throughout this pipeline. The name is "bacteria_manifest.tsv".

### REMEMBER!
Please remember to change the directory paths to your own specific path!

## Data Import & Quality Information
 This visualization provides a summary of sequence counts per sample and plots of sequence quality at each position. The output file will be a .qzv file, so it can be viewed in QIIME2 View (//https://view.qiime2.org/)

 
 
 Use this script under 'Needed Materials': data_import.sh
 
 ```{linux}
#!/bin/bash
#SBATCH -J Job Name
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
--input-path /home/peterfs/practice/files/bacteria_manifest.tsv \
--output-path /home/peterfs/practice/QIIME/data.qza \
--input-format PairedEndFastqManifestPhred33V2 \

qiime demux summarize \
--i-data /home/peterfs/practice/QIIME/data.qza \
--o-visualization /home/peterfs/practice/QIIME/qualityplot.qzv\
```

## Denoising via DADA2
The DADA2 (Divisive Amplicon Denoising Algorithm 2) will do the following steps:
- Quality filtering, trimming, error correction, dereplication, and chimera removal 
- Demultiplexing and denoising of raw sequence datasets in FASTQ format

 Use this script under 'Needed Materials': denoise.sh
 
```{linux}
#!/bin/bash
#SBATCH -J Job Name
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
```

## Clustering and Taxonomy Information
Taxonomy will be assigned based on your classifier. Here, we utilized the SILVA 138 classifier (available in Group 1's 'Needed Materials' folder.  In this final step, a taxa bar plot will be created to visualize community structure of your samples. 

 Use this script under 'Needed Materials': classify.sh 
```{linux}
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
```
Visualize the barplot on QIIME2 View (//https://view.qiime2.org/)

## Exporting Taxonomic Information
After visualization, feel free to export the OTU abundance file (TSV format) in the level you want for downstream analysis (e.g. Phylum = Level 2)
