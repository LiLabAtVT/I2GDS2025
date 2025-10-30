#!/bin/bash

#SBATCH -t 70:00:00
#SBATCH -p normal_q
#SBATCH -A introtogds
#SBATCH --mail-type=ALL
#SBATCH --mail-user=###vt-pid@vt.edu # Change to whichever email you would like to receive job updates
#SBATCH --cpus-per-task=4
#SBATCH --mem=200GB
#SBATCH --output=bwa_%j.out
#SBATCH --error=bwa_%j.err

#Path to main folder (where the github folders are located)
cd /projects/intro2gds/I2GDS2025/G4_Viruses/github/

#Set Conda Environment
source ~/.bashrc
conda activate g4_viruses

#Create an input and output directory for BWA samples, set the thread count, set reference database directory, and create a log
REF="/projects/intro2gds/I2GDS2025/G4_Viruses/databases/bwa/human_ref.fna"
INPUT_DIR="outputs/trimmed_outputs"
OUTPUT_DIR="outputs/bwa_outputs"
LOG_DIR="logs"
THREADS=16

LOGFILE="$LOG_DIR/bwa_${SLURM_JOB_ID}.log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"; }

log "Starting BWA filtering on $(hostname)"
log "Reference: $REF"
mkdir -p "$OUTPUT_DIR"

#Input previous trim_galore output files and run a loop using BWA to create SAM files that will be converted to BAM files then zip them
for FILE in "$INPUT_DIR"/*_trimmed.fq.gz; do
  [ -e "$FILE" ] || { log "No trimmed FASTQ found in $INPUT_DIR"; break; }
  SAMPLE=$(basename "$FILE" _trimmed.fq.gz)
  log "Processing $SAMPLE"

  SAM="${OUTPUT_DIR}/aln-${SAMPLE}.sam"
  SORTED="${OUTPUT_DIR}/aln-${SAMPLE}.sorted.bam"
  NONHOST="${OUTPUT_DIR}/non_host_reads_${SAMPLE}.bam"
  CLEANED="${OUTPUT_DIR}/cleaned_reads_${SAMPLE}.fastq"

  bwa mem "$REF" "$FILE" > "$SAM"
  samtools view -@ "$THREADS" -Sb "$SAM" | samtools sort -@ "$THREADS" -o "$SORTED"
  samtools index "$SORTED"
  samtools view -@ "$THREADS" -b -f 4 "$SORTED" > "$NONHOST"
  samtools fastq -@ "$THREADS" -0 "$CLEANED" "$NONHOST"

  gzip -f "$SAM" "$SORTED" "$SORTED.bai" "$NONHOST" "$CLEANED"
  log "Finished $SAMPLE"
done

log "BWA complete for all samples."
