#!/bin/bash
# -------------------------------------------
#bam_to_spades.sh
#SBATCH --account=introtogds
#SBATCH --job-name=download_SRR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mail-user=jingjingy@vt.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=200GB
#SBATCH --cpus-per-task=8
# -------------------------------------------

echo "Job started at $(date)"

# ------------------------------
# Set working directory
# ------------------------------
cd /projects/intro2gds/I2GDS2025/G2_PlantDisease/Jingjing/mapping_results

# ------------------------------
# enter parameters
# ------------------------------
BAM="SRR29108932.sorted.bam"            # input BAM file
OUTDIR="/projects/intro2gds/I2GDS2025/G2_PlantDisease/Jingjing/spades_output"         # SPAdes  output
THREADS=8
MEMORY=200


# ------------------------------
# Load modules
# ------------------------------
module load SAMtools
module load SPAdes

BASENAME="Historial_Xf"            # set prefix of output
FASTQ1="${BASENAME}_1.fastq"       # paired-end reads 1
FASTQ2="${BASENAME}_2.fastq"       # paired-end reads 2
                                                              
mkdir -p "$OUTDIR"

# ------------------------------
# Extract paired-end reads from BAM
# ------------------------------
echo "Extracting paired-end reads from $BAM ..."
samtools fastq -1 "$FASTQ1" -2 "$FASTQ2" -0 /dev/null -s /dev/null -n "$BAM"


# ------------------------------
# Assemble with SPAdes
# ------------------------------
echo "Running SPAdes assembly..."
spades.py --isolate \
    -1 "$FASTQ1" -2 "$FASTQ2" \
    -o "$OUTDIR" \
    --threads $THREADS \
    --memory $MEMORY

echo "Assembly finished. Output in $OUTDIR/"
