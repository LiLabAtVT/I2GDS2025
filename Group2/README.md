## Introduction
This project replicates a bioinformatics pipeline for analyzing a collection of Xylella fastidiosa (Xf) genomes obtained from both pure cultures and century-old herbarium specimens.

Reference: Century-old herbarium specimen provides insights into Pierce’s disease of grapevines emergence in the Americas
https://doi.org/10.1016/j.cub.2024.11.029

## Overview of workflow

| Step                   | Description                                                                                                        | Tool         |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------ | ------------ |
| **1. Data download**   | Retrieve metagenomic reads from herbarium specimens, raw reads of 44 modern *Xf* strains, and reference genome(s). | **SRAtools** |
| **2. Mapping**         | Map herbarium metagenomic reads to the *Xf* reference genome.                                                      | **Bowtie2 and SAMtools**  |
| **3. Assembly**        | Assemble mapped reads (herbarium) and raw reads (modern strains).                                                  | **SPAdes**   |
| **4. Quality control** | Assess completeness and contamination of assembled genomes.                                                        | **CheckM**   |
| **5. Annotation**      | Annotate genes in assembled genomes.                                                                               | **Prokka**   |


## Environment setup

## 1. Data download - SRAtools
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

echo "All downloads finished at $(date)"
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

## 2. Mapping - Bowtie2 and SAMtools
Mapping metagenomic reads from herbarium specimens to the Xylella fastidiosa reference genome allows selective recovery of pathogen-derived sequences from mixed plant and microbial DNA.
This step confirms the presence of X. fastidiosa in historical material, enables genome reconstruction, and provides quality metrics (e.g., alignment rate and coverage) essential for subsequent assembly and evolutionary analyses.

Bowtie2 is used to align metagenomic reads from herbarium specimens to the X. fastidiosa reference genome. Since X. fastidiosa is a xylem-limited bacterial pathogen, mapping reads to its genome allows the identification of ancient pathogen sequences preserved in historical plant tissue.

SAMtools was used to convert, sort, and manage the Bowtie2 alignment files.
It converts large SAM files into compressed BAM format, sorts alignments by genomic coordinates, and enables generation of mapping statistics and indexing for efficient downstream analysis.

```
#!/bin/bash
# -------------------------------------------
# map_metagenomes.sh
# Purpose: Map metagenomic reads from herbarium specimens 
#          to Xylella fastidiosa reference genome using Bowtie2 and SAMtools
# -------------------------------------------

#SBATCH --account=introtogds
#SBATCH --job-name=map_metagenomes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mail-user=jingjingy@vt.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=200GB
#SBATCH --cpus-per-task=8
# -------------------------------------------

echo "Job started at $(date)"
echo "Running on node: $(hostname)"

# ------------------------------
# 1️⃣ Set working directory
# ------------------------------
cd /projects/intro2gds/I2GDS2025/G2_PlantDisease/Jingjing || exit
echo "Working directory: $(pwd)"

# ------------------------------
# 2️⃣ Load required modules
# ------------------------------
module load Bowtie2
module load SAMtools

# ------------------------------
# 3️⃣ Define input/output paths
# ------------------------------
REF_FASTA=./RawData/Ref_genome/GCF_000007245.1_ASM724v1_genomic.fna    # Reference genome
METAGENOME_DIR=./RawData/Historical_fastq                              # Input metagenomic FASTQ files
OUTPUT_DIR=./mapping_results                                           # Output BAM directory
INDEX_PREFIX=./RawData/Ref_genome/xf_ref                               # Bowtie2 index prefix

mkdir -p "$OUTPUT_DIR"

# ------------------------------
# 4️⃣ Build Bowtie2 index (if not exists)
# ------------------------------
if [ ! -f ${INDEX_PREFIX}.1.bt2 ]; then
    echo "Building Bowtie2 index..."
    bowtie2-build "$REF_FASTA" "$INDEX_PREFIX"
else
    echo "Index already exists. Skipping build step."
fi

# ------------------------------
# 5️⃣ Map metagenomic reads
# ------------------------------
for f1 in "$METAGENOME_DIR"/*_1.fastq.gz; do
    fname=$(basename "$f1" _1.fastq.gz)
    f2="${METAGENOME_DIR}/${fname}_2.fastq.gz"

    # Check for paired-end files
    if [ ! -f "$f2" ]; then
        echo "Paired file for $f1 not found! Skipping..."
        continue
    fi

    echo "Mapping sample: $fname"

    # Bowtie2 alignment (output SAM)
    bowtie2 -x "$INDEX_PREFIX" -1 "$f1" -2 "$f2" \
        --very-sensitive -p "$SLURM_CPUS_PER_TASK" \
        -S "$OUTPUT_DIR/${fname}.sam"

    # Convert SAM → sorted BAM
    samtools view -bS "$OUTPUT_DIR/${fname}.sam" \
        | samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$OUTPUT_DIR/${fname}.sorted.bam"

    # Index BAM
    samtools index "$OUTPUT_DIR/${fname}.sorted.bam"

    # Remove SAM to save space
    rm "$OUTPUT_DIR/${fname}.sam"

    echo "$fname mapping complete. BAM: ${fname}.sorted.bam"
done

echo "All paired-end metagenomes successfully mapped!"
echo "Job finished at $(date)"

```

## 3. Genome assembly - SPAdes
Raw Illumina paired-end reads were assembled de novo using SPAdes v4.1.0 (Bankevich et al., 2012). Each isolate’s forward (*_1.fastq.gz) and reverse (*_2.fastq.gz) reads are assembled independently with the --only-assembler flag to disable read error correction. The resulting contigs output to data/assemblies/, and the primary assembly file (scaffolds.fasta) for each isolate will be used in downstream annotation.

SPAdes was chosen for its balance of accuracy and computational efficiency in assembling bacterial genomes from Illumina short reads, providing robust contigs/scaffolds suitable for subsequent annotation (Prokka), ortholog detection, and phylogenetic analysis.

It will be best to submit this step as a slurm job:
``` 
#!/bin/bash
#SBATCH --job-name=spades_batch
#SBATCH --output=logs/spades_%x_%j.log
#SBATCH --error=logs/spades_%x_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --account=introtogds


# Load SPAdes module
module load spades

# Stop on first error
set -e

# Directories
READS_DIR="data/Modern_strain"
OUT_BASE="Carter/data/assemblies"

# Create output + log directories
mkdir -p "${OUT_BASE}" logs

# Loop over all *_1.fastq.gz files
for R1 in ${READS_DIR}/*_1.fastq.gz; do
    SAMPLE=$(basename "${R1}" _1.fastq.gz)
    R2="${READS_DIR}/${SAMPLE}_2.fastq.gz"
    OUT_DIR="${OUT_BASE}/${SAMPLE}_spades"

    echo "======================================"
    echo "Starting SPAdes for ${SAMPLE}"
    echo "Output -> ${OUT_DIR}"
    echo "======================================"

    mkdir -p "${OUT_DIR}"

    spades.py \
        -1 "${R1}" \
        -2 "${R2}" \
        -o "${OUT_DIR}" \
        -t ${SLURM_CPUS_PER_TASK} \
        -m 32 \
        --isolate \
        --cov-cutoff auto

    echo "SPAdes assembly finished for ${SAMPLE}"
    echo "Main contigs file: ${OUT_DIR}/contigs.fasta"
done

echo "======================================"
echo "All assemblies completed!"
echo "======================================"
```


## 4. Quality control - CheckM
```

```

## 5. Genome annotation - Prokka
Assembled and quality-controlled contigs were annotated using Prokka v1.14.6 (Seemann, 2014), a rapid annotation pipeline designed for prokaryotic genomes. Each assembly (scaffolds.fasta) from the quality-controlled assembly folder all_bins/ was annotated independently in parallel using 8 CPU threads. The output for each genome was written to data/annotations/, generating standard annotation files including GFF3, GenBank, and FAA (protein) files. 

Prokka was selected for its speed, consistency, and compatibility with downstream comparative genomics workflows. Prokka identifies and functionally annotates genes, rRNAs, tRNAs, and other genomic features using curated databases such as UniProt, RefSeq, and Pfam. It produces high-quality, standardized annotations that enable reliable gene-based analyses.

```
#!/bin/bash
#SBATCH --job-name=prokka_batch
#SBATCH --output=logs/prokka_%x_%j.out
#SBATCH --error=logs/prokka_%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=normal_q
#SBATCH --account=introtogds

set -euo pipefail

# --- Load environment ---
module load Miniconda3/24.7.1-0
source activate ~/.conda/envs/prokka_env

# --- Set paths ---
ASSEMBLY_DIR="all_bins"       # change to correct location relative to script
OUT_BASE="data/annotations"

mkdir -p logs
mkdir -p "${OUT_BASE}"

# --- Change directory to where the job was submitted ---
cd "${SLURM_SUBMIT_DIR:-$PWD}"

# Enable nullglob so loop doesn’t fail if no matches
shopt -s nullglob

# --- Loop over all FASTA assemblies ---
for ASSEMBLY in ${ASSEMBLY_DIR}/SRR*.fasta; do
    SAMPLE=$(basename "${ASSEMBLY}" .fasta)
    OUT_DIR="${OUT_BASE}/${SAMPLE}_prokka"

    echo "==========================================="
    echo "Starting Prokka for sample: ${SAMPLE}"
    echo "Input file: ${ASSEMBLY}"
    echo "Output dir: ${OUT_DIR}"
    echo "==========================================="

    prokka \
        --outdir "${OUT_DIR}" \
        --prefix "${SAMPLE}" \
        --cpus ${SLURM_CPUS_PER_TASK:-8} \
        --force \
        --locustag "${SAMPLE}" \
        "${ASSEMBLY}"

    echo "Finished Prokka for ${SAMPLE}"
    echo
done

echo  "All Prokka annotations completed successfully!"
```




