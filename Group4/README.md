# Environment Setup
`interact -A introtogds -p normal_q -t 1:00:00`

# Software Installation

# Data Download

# Analysis Pipeline

## Trim Galore

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
</details>

## BWA

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
#SBATCH --output=bwa_%j.out
#SBATCH --error=bwa_%j.err

cd /projects/intro2gds/I2GDS2025/G4_Viruses/github/

#Set Conda Environment
source ~/.bashrc
conda activate g4_viruses

REF="/projects/intro2gds/I2GDS2025/G4_Viruses/databases/bwa/human_ref.fna"
INPUT_DIR="outputs/trimmed_outputs"
OUTPUT_DIR="outputs/bwa_outputs"
THREADS=16

LOGFILE="logs/bwa_filter_${SLURM_JOB_ID}.log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"; }

log "Starting BWA filtering on $(hostname)"
log "Reference: $REF"
mkdir -p "$OUTPUT_DIR"

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

log "All BWA filtering complete."
```

</details>

## Diamond

## Spades

## Fastq
