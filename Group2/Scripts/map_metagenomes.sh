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
