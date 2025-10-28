# Environment Setup
`interact -A introtogds -p normal_q -t 1:00:00`

# Software Installation

# Data Download

# Analysis Pipeline

# Trim Galore

Trim Galore is used to conduct a quality control on the sample files. 

<details>
  <summary>Click to expand code</summary>
  
```
#!/bin/bash

#SBATCH -t 70:00:00
#SBATCH -p normal_q
#SBATCH -A introtogds
#SBATCH --mail-type=ALL
#SBATCH --mail-user=###mitchellgercken@vt.edu
#SBATCH --cpus-per-task=4
#SBATCH --mem=200GB
#SBATCH --output=trim_galore_%j.out
#SBATCH --error=trim_galore_%j.err

#Path to main folder
cd /projects/intro2gds/I2GDS2025/G4_Viruses/github/

#Resolve script directory and repo root

INPUT_DIR="test_data"
OUTPUT_DIR="outputs/trimmed_outputs"
LOG_DIR="logs"
THREADS=8

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

#--- Logging function ---
LOGFILE="$LOG_DIR/trim_galore_${SLURM_JOB_ID}.log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"; }

log "Starting Trim Galore job on $(hostname)"
log "Input directory: $INPUT_DIR"
log "Output directory: $OUTPUT_DIR"

#--- Activate conda environment ---
source ~/.bashrc
conda activate g4_viruses

#--- Main loop ---

FASTQ_FILES=(test_data/sample*_test_data.fastq.gz)
[ ${#FASTQ_FILES[@]} -gt 0 ] || { log "No FASTQ files found in $INPUT_DIR"; exit 1; }

for FILE in "${FASTQ_FILES[@]}"; do
    SAMPLE="${FILE%%_test_data.fastq.gz}"
    log "Processing sample: $SAMPLE"
    
    trim_galore "$FILE" -j "$THREADS" -o "$OUTPUT_DIR"
    
    log "Finished sample: $SAMPLE"
done

log "Trim Galore completed successfully."
```


## BWA- Burrow-Wheeler Aligner that aligns DNA sequences against a large reference genome. This uses three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first is for Illumina reads up to 100 bp, the second and third uses reads up from 70 bp to a few megabases. 

## Diamond- Sequence aligner for protein and translated DNA searches. This uses pairwise alignment of proteins and translated DNA at 100x-10,000x speed of BLAST. 

## SPAdes- This is used for assembling and analyzing sequencing data from Illumina and IonTorrent technologies. This is for smaller genomes (bacteria) and not intended for larger genomes (larger eukaryotes). This supports paired-end reads, mate-pairs (including high quality Nextera Mate Pairs) and unpaired reads.

## Fastq
