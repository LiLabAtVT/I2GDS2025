## Introduction
This project replicates a bioinformatics pipeline for analyzing a collection of Xylella fastidiosa (Xf) genomes obtained from both pure cultures and century-old herbarium specimens.

Reference: Century-old herbarium specimen provides insights into Pierce’s disease of grapevines emergence in the Americas
https://doi.org/10.1016/j.cub.2024.11.029

## Overview of workflow

| Step                   | Description                                                                                                        | Tool         |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------ | ------------ |
| **1. Data download**   | Retrieve metagenomic reads from herbarium specimens, raw reads of 44 modern *Xf* strains, and reference genome(s). | **SRAtools** |
| **2. Mapping**         | Map herbarium metagenomic reads to the *Xf* reference genome.                                                      | **Bowtie2**  |
| **3. Assembly**        | Assemble mapped reads (herbarium) and raw reads (modern strains).                                                  | **SPAdes**   |
| **4. Quality control** | Assess completeness and contamination of assembled genomes.                                                        | **CheckM**   |
| **5. Annotation**      | Annotate genes in assembled genomes.                                                                               | **Prokka**   |

1. **Data download** - **SRAtools**
   - Metagenomic reads from herbarium specimen(s)
   - Raw reads of 44 modern Xf strains
   - Reference genome(s) of *X. fastidiosa*

2. **Mapping Xf in metagenomes** - **Bowtie2**
   - Use **Bowtie2** to map herbarium metagenomic reads to the Xf reference genome

3. **Assembly** - **SPAdes**
     - 2.1 Assemble mapped reads from herbarium specimen
     - 2.2 Assemble raw reads of modern strains

4. **Quality control** - **CheckM**
   - QC of assembled genomes/scaffolds 

5. **Annotation** - **Prokka** for gene annotation

## Environment setup

## Data download - SRAtools
1.1 Metagenome dataset (PRJNA1114123) and reference sequence genome (GCA_000007245.1) was directly downloaded from NCBI

1.2 The accession list of 44 modern strains was download from NCBI and then using SRAtools for downloading. 
```
#!/bin/bash
# -------------------------------------------
# Download_SRR.sh
#SBATCH --account=introtogds
#SBATCH --job-name=download_SRR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mail-user=jingjingy@vt.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=200GB
#SBATCH --cpus-per-task=4
# -------------------------------------------

set -euo pipefail
echo "Job started at $(date)"

# set path to local sratoolkit
export PATH=/projects/intro2gds/I2GDS2025/G2_PlantDisease/Jingjing/sratoolkit.3.2.1-ubuntu64/bin:$PATH

# 2. set working directory 
OUTDIR=/projects/intro2gds/I2GDS2025/G2_PlantDisease/Jingjing/RawData
LIST=${OUTDIR}/sra_list.tx
cd "$OUTDIR"

# 3. Batch download
while read ACC; do
    echo "=== Processing $ACC ==="
    prefetch --output-directory "$OUTDIR" "$ACC"
    fasterq-dump --threads $SLURM_CPUS_PER_TASK --outdir "$OUTDIR" "$ACC"
    gzip "${OUTDIR}/${ACC}"*.fastq
    echo "=== $ACC done ==="
done < "$LIST"

echo "✅ All downloads finished at $(date)"
```


1.3 As downloaded files were in ".sra" format, command "fasterq-dump" was used for transfering SRA to FASTQ 
```
#!/bin/bash
# -------------------------------------------
# sra_to_fastq.sh
#SBATCH --account=introtogds
#SBATCH --output=/projects/intro2gds/I2GDS2025/G2_PlantDisease/Jingjing/script/logs/sra_to_fastq.%j.out
#SBATCH --error=/projects/intro2gds/I2GDS2025/G2_PlantDisease/Jingjing/script/logs/sra_to_fastq.%j.err
#SBATCH --job-name=sra_to_fastq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mail-user=jingjingy@vt.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=200GB
#SBATCH --cpus-per-task=4
# -------------------------------------------

SRA_FILE="/projects/intro2gds/I2GDS2025/G2_PlantDisease/Jingjing/RawData/SRR29108932"

FASTQ_DIR="/projects/intro2gds/I2GDS2025/G2_PlantDisease/Jingjing/RawData/fastq_historical"
mkdir -p "$FASTQ_DIR"

THREADS=8

echo "Converting $SRA_FILE to FASTQ ..."

export PATH="/projects/intro2gds/I2GDS2025/G2_PlantDisease/Jingjing/sratoolkit.3.2.1-ubuntu64/bin:$PATH"

#  use fasterq-dump transfer SRA to FASTQ
fasterq-dump "$SRA_FILE" \
  --split-files \
  --threads "$THREADS" \
  -O "$FASTQ_DIR"

#  gzip fastq
echo "Compressing FASTQ files ..."
gzip -f "$FASTQ_DIR"/*.fastq

echo "Done! Compressed FASTQ files are in $FASTQ_DIR"
```


## Mapping Xf genome in herbarium sample - Bowtie2
As Xf is a bacterial pathogen colonizing in xylem vessel, it is possible to sequece and assemble its historical genome from historical material. Reads were mapped against the reference geonme sequence using Bowtie2.




## SPAdes


## CheckM


## Prokka
