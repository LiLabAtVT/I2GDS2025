# Introduction to Group 1 16S Amplicon Pipeline
Our article examines the relative abundance, taxonomic profiles, and community structure of bacterial and fungal communities associated with parsley (*Petroselinum crispum*) and celery (*Apium graveolens*) roots via monocropping and intercropping systems.  The study aims to provide a baseline understanding of how intercropping influences rhizosphere microbial dynamics.

Article: https://bmcgenomdata.biomedcentral.com/articles/10.1186/s12863-025-01351-0#additional-information

This pipeline utilizes QIIME2's (Quantitative Insights Into Microbial Ecology 2) "Amplicon" Distribution, a suite of plug-ins that provide broad analytic functionality to support microbiome marker gene analysis from raw sequencing data to pubication-quality visualizations and statistics.
We expect to reproduce taxonomic composition tables for both bacteria and fungi.

Figure from article we plan to reproduce:

<img width="350" height="200" alt="image" src="https://github.com/user-attachments/assets/8210a433-4b7a-4c40-aecd-a499ed9a6afc" />

Babalola et al. 2025

## Installing QIIME2 into Conda environment
Installing the base distribution's conda environment: QIIME 2 recommends creating a new environment specifically for the QIIME 2 distribution, as they are many dependencies that you wouldn't want in an existing environment.
```{linux}
# Installation for Linux/Windows users
# Load Miniconda3 (or other verison of conda) 
module load Miniconda3/24.7.1-0
conda env create \
  --name qiime2-amplicon-2025.7 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.7/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml
```

```{linux}
# Installation for MacOS users (Apple Silicon)
CONDA_SUBDIR=osx-64 conda env create \
  --name qiime2-amplicon-2025.7 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.7/amplicon/released/qiime2-amplicon-macos-latest-conda.yml
conda activate qiime2-amplicon-2025.7
conda config --env --set subdir osx-64
```
Verify your installation:
```{linux}
# Test your installation
conda deactivate
conda activate qiime2-amplicon-2025.7
qiime info
```

Activate your QIIME2 version (after verification)
```{linux}
conda activate qiime2-2025.7 # reference your specific version of QIIME
OR
source activate qiime2-amplicon-2025.7
```

## Test Data
The test data is labeled "Group1_TestData" in the following ARC Directory: /projects/intro2gds/I2GDS2025/TestData_LinuxPeerEval
It consists of **12 bacterial samples** (forward & reverse reads) in FASTQ format (.fastq.gz), which have been deposited in the Sequence Read Archive (SRA) of the National Center for Biotechnology Information (NCBI) under the Bioproject Accession numbers; **SRP540554 (16S rRNA)** and **SRP540675 (ITS)**.

We downloaded the provided script to retrieve each samples' reads: ena-file-download-read_run-SRP540554-fastq_ftp-20250915-1609.sh
This tutorial will focus on the bacterial samples (Link: http://identifiers.org/insdc.sra:SRP540554).

Here are the provided materials for this tutorial before we move to QIIME2:
- Manifest file ('Needed Materials' Folder)
- Test data ('Group1_testdata' ARC Directory)
- Classifier ('Group1_testdata' ARC Directory)
- Bash scripts ('Needed Materials' Folder)
# QIIME2 Pipeline
The goal is to generate high-resolution OTUs (operational taxonomic unit). OTUs are employed to classify groups of similar sequences at 97% similarity.
We have provided our manifest file (TSV format) into our test data; it will be used throughout this pipeline. The name is "bacteria_manifest.tsv".

## Visual Pipeline
<img width="300" height="700" alt="16S Visual Pipeline" src="https://github.com/user-attachments/assets/be07b3e2-4911-48c2-8fda-84a66b262f2d" />

## REMEMBER!
Please remember to change the directory paths to your own specific path!

## Data Import & Quality Information
 This visualization provides a summary of sequence counts per sample and plots of sequence quality at each position. The output file will be a .qzv file, so it can be viewed in QIIME2 View (https://view.qiime2.org/)

 Use this script under 'Needed Materials': data_import.sh

INPUT: "bacteria_manifest.tsv"

OUTPUT: "qualityplot.qzv"

 <img width="800" height="450" alt="image" src="https://github.com/user-attachments/assets/9bfb7530-2794-40fd-927b-670a64fb307a" />

 ```{linux}
#!/bin/bash
#SBATCH -J Data Import Quality Plot
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
```

## Denoising via DADA2
The DADA2 (Divisive Amplicon Denoising Algorithm 2) will do the following steps:
- Quality filtering, trimming, error correction, dereplication, and chimera removal 
- Demultiplexing and denoising of raw sequence datasets in FASTQ format

 Use this script under 'Needed Materials': denoise.sh
 
INPUT: "data.qza"

OUTPUT: "sequences.qza", "table.qza", "denoising_stats.qza", and "denoising_stats.qzv"
 
```{linux}
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
```
Generate a .qzv file of denoising_stats to visualize results in QIIME2 View 
```{linux}
qiime metadata tabulate \
  --m-input-file /path/to/your/directory/denoising_stats.qza \
  --o-visualization /path/to/your/directory/denoising_stats.qzv
```

## Clustering and Taxonomy Information
Taxonomy will be assigned based on your classifier. Here, we utilized the SILVA 138 classifier (available in Group 1's 'Needed Materials' folder.  In this final step, a taxa bar plot will be created to visualize community structure of your samples. 

 Use this script under 'Needed Materials': classify.sh 
 
 INPUT:
 
 - Step 1: "sequences.qza" and "table.qza"
 - Step 2: "rep_seqs_97.qza" and "silva-138-99-nb-classifier.qza"
 - Step 3: "updated_taxonomy.qza"
 - Step 4: "table_97.qza"  and "updated_taxonomy.qza"
 - Step 5: "tax-class-filter-table.qza"
 - Step 6: "feature-frequency-filtered-table.qza", "bacteria_manifest.tsv", and "updated_taxonomy.qza"
   
OUTPUT:

 - Step 1: "table_97.qza" and "rep_seqs_97.qza"
 - Step 2: "updated_taxonomy.qza"
 - Step 3: "taxa-meta.qzv"
 - Step 4: "tax-class-filter-table.qza"
 - Step 5: "feature-frequency-filtered-table.qza"
 - Step 6: "taxa-barplot.qzv"
   
```{linux}
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
```

Expected output:

<img width="550" height="350" alt="image" src="https://github.com/user-attachments/assets/b3d36d1d-7a41-4167-978a-80bf6419be70" />

Visualize the barplot on QIIME2 View (//https://view.qiime2.org/)
## Exporting Taxonomic Information
After visualization, export the OTU abundance file (TSV format) in the level you want for downstream analysis (e.g. Phylum = Level 2)

```{linux}
qiime tools export --input-path feature-frequency-filtered-table.qza --output-path exported
qiime tools export --input-path updated_taxonomy.qza --output-path exported
```
```{linux}
# move into new directory "exported" 
cp taxonomy.tsv biom-taxonomy.tsv
```
```{linux}
# Change the first line of biom-taxonomy.tsv (i.e. the header) to this:
#    #OTUID taxonomy confidence
# (use file editor in ARC dashboard)
biom add-metadata \
-i feature-table.biom \
-o table-with-taxonomy.biom \
--observation-metadata-fp biom-taxonomy.tsv \
--sc-separated taxonomy
```
```{linux}
qiime taxa collapse \
 --i-table /path/to/your/directory/feature-frequency-filtered-table.qza \
 --i-taxonomy /path/to/your/directory/updated_taxonomy.qza \
 --p-level 2 \
 --o-collapsed-table level2-table.qza
```
```{linux}
qiime tools export \
 --input-path level2-table.qza \
 --output-path exported_table
```
```{linux}
# move into new directory "exported_table"
biom convert \
 -i /path/to/your/directory/exported_table/feature-table.biom \
 -o /path/to/your/directory/exported_table/level2-table.tsv \
 --to-tsv
```

# References
### QIIME2: ​

https://doi.org/10.1038/s41587-019-0209-9​

### QIIME2 Code: 

Riddley, M. "Mia's QIIME2 Workflow" [Mia's QIIME2 Workflow.pdf](https://github.com/user-attachments/files/23196285/Mia.s.QIIME2.Workflow.pdf)


### SILVA:​

Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO (2013) The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucl. Acids Res. 41 (D1): D590-D596.​
Yilmaz P, Parfrey LW, Yarza P, Gerken J, Pruesse E, Quast C, Schweer T, Peplies J, Ludwig W, Glöckner FO (2014) The SILVA and "All-species Living Tree Project (LTP)" taxonomic frameworks. Nucl. Acids Res. 42:D643-D648​
Glöckner FO, Yilmaz P, Quast C, Gerken J, Beccati A, Ciuprina A, Bruns G, Yarza P, Peplies J, Westram R, Ludwig W (2017) 25 years of serving the community with ribosomal RNA gene reference databases and tools. J. Biotechnol.
