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


## Environment setup

## Data download - SRAtools
1.1 Metagenome dataset (PRJNA1114123) and reference sequence genome (X. fastidiosa Temecula1: GCA_000007245.1) were both downloaded directly from NCBI

1.2 Modern isolates: Accession list of 44 Xf strains was retrieved from NCBI and downloaded using SRAtools. 
Note: Make sure sra_list.txt contains one SRA accession per line (e.g., SRR12345678).
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

# 1. Add SRA Toolkit to PATH
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
#SBATCH --cpus-per-task=8
# -------------------------------------------

SRA_FILE="/projects/intro2gds/I2GDS2025/G2_PlantDisease/Jingjing/RawData/SRR29108932"
FASTQ_DIR="/projects/intro2gds/I2GDS2025/G2_PlantDisease/Jingjing/RawData/fastq_historical"
mkdir -p "$FASTQ_DIR"
THREADS=8

echo "Converting $SRA_FILE to FASTQ ..."

# 1. Add SRA Toolkit to PATH
export PATH="/projects/intro2gds/I2GDS2025/G2_PlantDisease/Jingjing/sratoolkit.3.2.1-ubuntu64/bin:$PATH"

# 2. Convert SRA to FASTQ
fasterq-dump "$SRA_FILE" \
  --split-files \
  --threads "$THREADS" \
  -O "$FASTQ_DIR"

# 3. Compress FASTQ
echo "Compressing FASTQ files ..."
gzip -f "$FASTQ_DIR"/*.fastq

echo "Done! Compressed FASTQ files are in $FASTQ_DIR"
```

## Mapping - Bowtie2
As Xylella fastidiosa is a xylem-limited bacterial pathogen, it is possible to reconstruct its historical genome from preserved herbarium tissue.
Metagenomic reads from historical samples were mapped against the Xf reference genome using Bowtie2.




## Genome assembly - SPAdes
Raw Illumina paired-end reads were assembled de novo using SPAdes v4.1.0 (Bankevich et al., 2012). Each isolate’s forward (*_1.fastq.gz) and reverse (*_2.fastq.gz) reads are assembled independently with the --only-assembler flag to disable read error correction. The resulting contigs output to data/assemblies/, and the primary assembly file (scaffolds.fasta) for each isolate will be used in downstream annotation.

SPAdes was chosen for its balance of accuracy and computational efficiency in assembling bacterial genomes from Illumina short reads, providing robust contigs/scaffolds suitable for subsequent annotation (Prokka), ortholog detection, and phylogenetic analysis.

It will be best to submit this step as a slurm job:
``` 


## Quality control - CheckM


## Genome annotation - Prokka
Assembled and quality-controlled contigs were annotated using Prokka v1.14.6 (Seemann, 2014), a rapid annotation pipeline designed for prokaryotic genomes. Each assembly (scaffolds.fasta) from the quality-controlled assembly folder all_bins/ was annotated independently in parallel using 8 CPU threads. The output for each genome was written to data/annotations/, generating standard annotation files including GFF3, GenBank, and FAA (protein) files. 

Prokka was selected for its speed, consistency, and compatibility with downstream comparative genomics workflows. Prokka identifies and functionally annotates genes, rRNAs, tRNAs, and other genomic features using curated databases such as UniProt, RefSeq, and Pfam. It produces high-quality, standardized annotations that enable reliable gene-based analyses.






