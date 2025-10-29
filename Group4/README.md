# Project Overview

The goal of this project is to replicate the pipeline as seen in 

Nousias, Orestis, Mark McCauley, Maximilian R. Stammnitz, et al. 2025. “Shotgun Sequencing of Airborne eDNA Achieves Rapid Assessment of Whole Biomes, Population Genetics and Genomic Variation.” Nature Ecology & Evolution 9 (6): 1043–60. https://doi.org/10.1038/s41559-025-02711-w.

<img width="486" height="717" alt="image" src="https://github.com/user-attachments/assets/4bbc007f-2aa5-410c-9c36-12747132898e" />



# Environment Setup
Create an ARC environment to complete jobs with

`interact -A introtogds -p normal_q -t 1:00:00`

# Software Installation
Software was downloaded via Conda. All of the required packages needed to run the pipeline are listed within the environment.yml file in the materials directory. 

Initialize Conda on ARC:
```
module load Miniconda3
```
To create the Conda environment:
```
conda env create -f environment.yml -n g4_viruses
```

# Analysis Pipeline

The following steps outline the main pipeline to analyze the reads from the reference paper. The final results will allow you to see species types in each sample. 

## Trim Galore

Trim Galore is used to clean high-throughput sequencing reads by automatically trimming adapters and low-quality bases. It serves as a wrapper around Cutadapt and FastQC, combining adapter removal with quality control checks in a single step. By removing unwanted sequences and short or poor-quality reads, Trim Galore improves the overall accuracy and reliability of downstream analyses such as read alignment and assembly. 

<details>
  <summary>Click to expand code</summary>
  
```
#!/bin/bash

#SBATCH -t 70:00:00
#SBATCH -p normal_q
#SBATCH -A introtogds
#SBATCH --mail-type=ALL
#SBATCH --mail-user=###vt-pid@vt.edu # Change to whichever email you would like to receive job updates
#SBATCH --cpus-per-task=4
#SBATCH --mem=200GB
#SBATCH --output=trim_galore_%j.out
#SBATCH --error=trim_galore_%j.err

#Path to main folder (where the github folders are located)
cd /projects/intro2gds/I2GDS2025/G4_Viruses/github/

#Set variables for loop

#create an input and output directory for trim_galore samples, set the thread count, and create a log
INPUT_DIR="test_data"
OUTPUT_DIR="outputs/trimmed_outputs"
LOG_DIR="logs"
THREADS=8

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

#Logging function
LOGFILE="$LOG_DIR/trim_galore_${SLURM_JOB_ID}.log"
#have log set exact date and time for each iteration
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"; }

log "Starting Trim Galore job on $(hostname)"
log "Input directory: $INPUT_DIR"
log "Output directory: $OUTPUT_DIR"

#Activate conda environment
source ~/.bashrc
conda activate g4_viruses

#Main loop
#Input test data files and run trim_galore on them, outputting them to a new directory
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

Burrow-Wheeler Aligner for short-read alignment. This maps DNA sequences against a large reference genome, such as the human genome. This uses 3 algorithms- BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence reads up to 100bp, while the other two for longer sequences ranged from 70bp to a few megabases.

<details>
  <summary>Click to expand code</summary>

```
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
THREADS=16

LOGFILE="logs/bwa_filter_${SLURM_JOB_ID}.log"
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

log "All BWA filtering complete."
```

</details>

## SPAdes 
SPAdes is used for analyzing sequencing data from Illumina and IonTorrent technologies. This is for smaller genomes (bacteria) and not intended for larger genomes (larger eukaryotes). This supports paired-end reads, mate-pairs and unpaired reads. 
<details>
  <summary>Click to expand code</summary>

```
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
THREADS=16

LOGFILE="logs/spades_${SLURM_JOB_ID}.log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"; }

log "Starting SPAdes assemblies"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || exit

for FILE in outputs/bwa_outputs/cleaned_reads_sample*_test_data.fastq.gz; do
  [ -e "$FILE" ] || { log "No cleaned reads found"; break; }
  SAMPLE=$(basename "$FILE" .fastq.gz | sed 's/cleaned_reads_//')
  log "Running SPAdes for $SAMPLE"

  spades.py --s1 "$FILE" -t "$THREADS" -o "${SAMPLE}" --only-assembler
  log "Finished SPAdes for $SAMPLE"
done

log "SPAdes assemblies complete."
```

</details>

## DIAMOND 
DIAMOND is a sequence aligner for protein and translated DNA searches. This uses pairwise alignment of proteins and translated DNA at 100x-10,000x speed of BLAST. 
<details>
  <summary>Click to expand code</summary>
  
```
#!/bin/bash

#SBATCH -t 70:00:00
#SBATCH -p normal_q
#SBATCH -A introtogds
#SBATCH --mail-type=ALL
#SBATCH --mail-user=###vt-pid@vt.edu # Change to whichever email you would like to receive job updates
#SBATCH --cpus-per-task=4
#SBATCH --mem=200GB
#SBATCH --output=diamond_%j.out
#SBATCH --error=diamond_%j.err

#Path to main folder (where the github folders are located)
cd /projects/intro2gds/I2GDS2025/G4_Viruses/github/

#Set conda environment
source ~/.bashrc
conda activate g4_viruses

#Create an output directory for diamond samples, set the input database, set the thread count, and create a log
DB="/projects/intro2gds/I2GDS2025/G4_Viruses/databases/diamond/nr"
SPADES_DIR="outputs/spades_outputs"
THREADS=16

LOGFILE="diamond_${SLURM_JOB_ID}.log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"; }

log "Starting DIAMOND BLASTx job on $(hostname)"
log "Database: $DB"

for CONTIG in "$SPADES_DIR"/sample*_test_data/contigs.fasta; do
  [ -e "$CONTIG" ] || { log "No contigs.fasta found in $SPADES_DIR"; break; }

  SAMPLE_DIR=$(basename "$(dirname "$CONTIG")")
  SAMPLE="${SAMPLE_DIR%%_test_data}"
  OUT_FILE="${SAMPLE_DIR}_assembly_test_data.daa"

  log "Running DIAMOND for $SAMPLE"

  diamond blastx \
    -d "$DB" \
    -q "$CONTIG" \
    -o "$OUT_FILE" \
    -p "$THREADS" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames stitle \
    --max-target-seqs 1 \
    --more-sensitive

  log "Finished DIAMOND for $SAMPLE"
done

log "DIAMOND BLASTx complete."
```

</details>

## Kraken2- 
Kraken2 is a very fast way to assign taxonomic labels using k-mers to metagenomic DNA sequences. Kraken2 splits sequences into smaller fragments of DNA as "k-mers". The k-mers are then compared in a hashing table to determine similarity to reference genomes in the database.

<details>
  <summary>Click to expand code</summary>

```
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

#Path to main folder (where the github folders are located)
cd /projects/intro2gds/I2GDS2025/G4_Viruses/github/

#Set conda environment
source ~/.bashrc

conda activate g4_viruses

#Parameters
DB="/projects/intro2gds/I2GDS2025/G4_Viruses/databases/kraken2/k2_db"
SPADES_DIR="outputs/spades_outputs"
THREADS=16

#Log file setup
LOGFILE="logs/kraken2_${SLURM_JOB_ID}.log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"; }

log "Starting Kraken2 classification job on $(hostname)"
log "Using database: $DB"
log "Scanning SPAdes assemblies in: $SPADES_DIR"

#Main loop
for CONTIG_PATH in "${SPADES_DIR}"/sample*_test_data/contigs.fasta; do
    #Skip if no files found
    [ -e "$CONTIG_PATH" ] || { log "No contigs.fasta files found in $SPADES_DIR"; break; }

    #Extract sample name (e.g., sample1_test_data)
    SAMPLE_DIR=$(basename "$(dirname "$CONTIG_PATH")")
    SAMPLE="${SAMPLE_DIR%%_test_data}"

    log "Processing sample: $SAMPLE"

    #Define output directory and files
    OUT_DIR="${SAMPLE_DIR}"
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
        --classified-out "$CLASSIFIED"

    #Compress large outputs
    gzip -f "$OUTPUT" "$CLASSIFIED"

    log "Finished processing $SAMPLE"
    log "--------------------------------"
done

log "All samples processed successfully."
```

</details>

## For Our Reviewers
The scripts needed to run the pipeline are located in the scripts/ directory. The expected results from each step are located in the expected_output_files/ directory, followed by each step. The scripts should be set up so you have the option of how many steps of the pipeline you would like to run. **Each script uses the outputs from the last. If you want to skip a step you can set the input directory variable to expected_outputs/your-step_outputs**. If you would like to get the evaluation done quickly, steps 2 and 3 take around ~30mins to complete. **Please change the code in each SLURM script under the #Path to main folder line to be the directory where your repository folder is downloaded. It will be different for everyone.**
**To download repository**
```
git clone https://github.com/LiLabAtVT/I2GDS2025/Group4/repository
```
The scripts are intended to be run in the logs directory. Each script will output three files: **1)** the default SLURM .out file, **2)** a .log file for the process of each script to be recorded, and **3)** a .err file for script errors to be reported if any occur. These three files will be named according to each step. For example, if you run step 2 (Trim Galore) and step 3 (BWA) the files will be named trim_galore_JOBID.out and bwa_JOBID.out. If anyone has an issue with the BWA, Kraken2, or DIAMOND scripts failing due to not having permission to any of the reference databases, please email mitchellgercken@vt.edu requesting access. 
If you would like to run the entire pipeline from Trim Galore to DIAMOND and/or Kraken2, a master_pipeline.sh script is included for your convenience.
