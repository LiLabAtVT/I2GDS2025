#!/bin/bash

#SBATCH -t 70:00:00
#SBATCH -p normal_q
#SBATCH -A introtogds
#SBATCH --mail-type=ALL
#SBATCH --mail-user=###vt-pid@vt.edu # Change to whichever email you would like to receive job updates
#SBATCH --cpus-per-task=4
#SBATCH --mem=200GB
#SBATCH --output=spades_%j.out
#SBATCH --error=spades_%j.err

#Path to main folder (where the github folders are located)
cd /projects/intro2gds/I2GDS2025/G4_Viruses/github/

#Set Conda Environment
source ~/.bashrc
conda activate g4_viruses

#Create an input and output directory for SPAdes samples, set the thread count, and create a log
INPUT_DIR="outputs/bwa_outputs"
OUTPUT_DIR="outputs/spades_outputs"
LOG_DIR="logs"
THREADS=16

LOGFILE="$LOG_DIR/spades_${SLURM_JOB_ID}.log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"; }

log "Starting SPAdes assemblies"
mkdir -p "$OUTPUT_DIR"

for FILE in "$INPUT_DIR"/cleaned_reads_sample*_test_data.fastq.gz; do
  [ -e "$FILE" ] || { log "No cleaned reads found"; break; }
  SAMPLE=$(basename "$FILE" .fastq.gz | sed 's/cleaned_reads_//')
  log "Running SPAdes for $SAMPLE"

  spades.py --s1 "$FILE" -t "$THREADS" -o "$OUTPUT_DIR"/"${SAMPLE}" --only-assembler
  log "Finished SPAdes for $SAMPLE"
done

log "SPAdes assemblies complete."
