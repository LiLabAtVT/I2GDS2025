#!/bin/bash

#SBATCH -t 70:00:00
#SBATCH -p normal_q
#SBATCH -A introtogds
#SBATCH --mail-type=ALL
#SBATCH --mail-user=###vt-pid@vt.edu # Change to whichever email you would like to receive job updates
#SBATCH --cpus-per-task=4
#SBATCH --mem=200GB
#SBATCH --output=kraken2_%j.out
#SBATCH --error=kraken2_%j.err

#Set downloaded directory (where the github folders are located)
cd /projects/intro2gds/I2GDS2025/G4_Viruses/github/

#Set conda environment
source ~/.bashrc
conda activate g4_viruses

#Parameters
DB="/projects/intro2gds/I2GDS2025/G4_Viruses/databases/kraken2/k2_db"
SPADES_DIR="outputs/spades_outputs"
OUTPUT_BASE="outputs/kraken2_outputs"
LOG_DIR="logs"
THREADS=16

#Logging setup
LOGFILE="$LOG_DIR/kraken2_${SLURM_JOB_ID:-manual}.log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"; }

log "Starting Kraken2 classification job on $(hostname)"
log "Using database: $DB"
log "Scanning SPAdes assemblies in: $SPADES_DIR"

#Main loop
for CONTIG_PATH in "${SPADES_DIR}"/sample*_test_data/contigs.fasta; do
    #Skip if no files found
    [ -e "$CONTIG_PATH" ] || { log "No contigs.fasta files found in $SPADES_DIR"; break; }

    #Extract sample name (e.g. sample1_test_data)
    SAMPLE_DIR=$(basename "$(dirname "$CONTIG_PATH")")
    SAMPLE="${SAMPLE_DIR%%_test_data}"

    log "Processing sample: $SAMPLE"

    #Define per-sample output directory under kraken2_outputs
    OUT_DIR="${OUTPUT_BASE}/${SAMPLE_DIR}"
    mkdir -p "$OUT_DIR"

    #Define output file paths
    REPORT="${OUT_DIR}/${SAMPLE}_assembly_report_test_data.txt"
    OUTPUT="${OUT_DIR}/${SAMPLE}_assembly_kraken_test_data.out"
    CLASSIFIED="${OUT_DIR}/${SAMPLE}_assembly_classified_test_data.fastq"

    #Run Kraken2 classification
    k2 classify \
        --db "$DB" \
        "$CONTIG_PATH" \
        --threads "$THREADS" \
        --report "$REPORT" \
        --output "$OUTPUT" \
        --classified-out "$CLASSIFIED" \
        2>&1 | tee -a "$LOGFILE"

    #Compress large outputs
    gzip -f "$OUTPUT" "$CLASSIFIED"

    log "Finished processing $SAMPLE"
    log "--------------------------------"
done

log "All samples processed successfully."
