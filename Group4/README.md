# Project Overview

The goal of this project is to replicate the pipeline as seen in:

Nousias, Orestis, Mark McCauley, Maximilian R. Stammnitz, et al. 2025. “Shotgun Sequencing of Airborne eDNA Achieves Rapid Assessment of Whole Biomes, Population Genetics and Genomic Variation.” Nature Ecology & Evolution 9 (6): 1043–60. https://doi.org/10.1038/s41559-025-02711-w.

We aim to create figures indicating total read count per species of the samples, indicating relative abundance. 

<img width="486" height="717" alt="image" src="https://github.com/user-attachments/assets/4bbc007f-2aa5-410c-9c36-12747132898e" />



# Environment Setup
Create an ARC environment to complete jobs with

`interact -A introtogds -p normal_q -t 1:00:00`

# Software Installation
Software was downloaded via Conda. All of the required packages needed to run the pipeline are listed within the environment.yml file in the repository/ directory. 

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
INPUT_DIR="/projects/intro2gds/I2GDS2025/TestData_LinuxPeerEval/G4_testdata"
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
FASTQ_FILES=("$INPUT_DIR"/sample*_test_data.fastq.gz)
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
LOG_DIR="logs"
THREADS=16

#have log set exact date and time for each iteration
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
LOG_DIR="logs"
THREADS=16

#have log set exact date and time for each iteration
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

#Set downloaded directory (where the github folders are located)
cd /projects/intro2gds/I2GDS2025/G4_Viruses/github/

#Set conda environment
source ~/.bashrc
conda activate g4_viruses

#Create an output directory for diamond samples, set the input database, set the thread count, and create a log
DB="/projects/intro2gds/I2GDS2025/G4_Viruses/databases/diamond/nr"
SPADES_DIR="outputs/spades_outputs"
OUTPUT_BASE="outputs/diamond_outputs"
LOG_DIR="logs"
THREADS=16

#have log set exact date and time for each iteration
LOGFILE="$LOG_DIR/diamond_${SLURM_JOB_ID}.log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"; }

log "Starting DIAMOND BLASTx job on $(hostname)"
log "Database: $DB"

for CONTIG_PATH in "$SPADES_DIR"/sample*_test_data/contigs.fasta; do
  [ -e "$CONTIG_PATH" ] || { log "No contigs.fasta found in $SPADES_DIR"; break; }

    #Extract sample name (e.g. sample1_test_data)
    SAMPLE=$(basename "$(dirname "$CONTIG_PATH")")
    SAMPLE_DIR="${SAMPLE}"
    
    #Define per-sample output directory under diamond_outputs
    OUT_DIR="${OUTPUT_BASE}/${SAMPLE_DIR}"
    mkdir -p "$OUT_DIR"
    OUTPUT="${OUT_DIR}/${SAMPLE}_assembly_diamond_test_data.daa"

  log "Running DIAMOND for $SAMPLE_DIR"

  diamond blastx \
    -d "$DB" \
    -q "$CONTIG" \
    -o "$OUTPUT" \
    -p "$THREADS" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames stitle \
    --max-target-seqs 1 \
    --more-sensitive

  log "Finished DIAMOND for $SAMPLE_DIR"
done

log "DIAMOND BLASTx complete."
```

</details>

## Kraken2
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
#have log set exact date and time for each iteration
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
```

</details>

## For Our Reviewers
The scripts needed to run the pipeline are located in the scripts/ directory. The expected results from each step are located in the expected_output_summary.txt file. The scripts should be set up so you have the option of how many steps of the pipeline you would like to run (except for the DIAMOND script, the script takes around 5 hours to complete). **Each script uses the outputs from the last.** If you would like to get the evaluation done quickly, steps 2 and 3 take around ~30mins to complete. **Please change the code in each SLURM script under the #Path to main folder line to be the directory where your repository folder is downloaded. It will be different for everyone.**

**To download repository:**
```
git clone --no-checkout https://github.com/LiLabAtVT/I2GDS2025/
cd I2GDS2025
git sparse-checkout init --cone
git sparse-checkout set Group4
git checkout
cd Group4/repository
```

The scripts are intended to be run in the logs directory. Each script will output three files: **1)** the default SLURM .out file, **2)** a .log file for the process of each script to be recorded, and **3)** a .err file for script errors to be reported if any occur. These three files will be named according to each step. For example, if you run steps 2 (Trim Galore) and 3 (BWA) the files will be named trim_galore_JOBID.out and bwa_JOBID.out. If anyone has an issue with the BWA or Kraken2 scripts failing due to not having permission to any of the reference databases, please email mitchellgercken@vt.edu requesting access. 

## R Script
The R script below generates a pie chart for both the Kraken2 and DIAMOND output files from the Linux-based scripts above. The R packages tidyverse, ggrepel, ggplot2, and data.table are required for the script to run.
<details>
  <summary>Click to expand code</summary>
  
```
library(tidyverse)
library(data.table)
library(ggrepel)
library(ggplot2)

indir <- "PATH/TO/INPUTDIR" # Set to input directory

# List input files
kraken_files <- list.files(indir, pattern = "_kraken_subset_1m.out$", full.names = TRUE) # read kraken2 input files
taxmap_files <- list.files(indir, pattern = "_tax-map.tsv$", full.names = TRUE) # read taxonomy file for kraken2 files
diamond_files <- list.files(indir, pattern = "_subset_1m_assembly.daa$", full.names = TRUE) # read DIAMOND input files

# Helper function to get sample name
get_sample_name <- function(filepath) {
  gsub("_assembly_kraken_subset_1m.out$", "", basename(filepath)) # extract sample name from input files
}

# Settings
desired_slices <- 10      # total slices including "Other"
other_min_pct <- 10       # minimum percentage for the Other slice

for (kraken_file in kraken_files) {

  sample_name <- get_sample_name(kraken_file)
  message("\nProcessing sample: ", sample_name)

  # Find matching taxmap and diamond files
  taxmap_file <- taxmap_files[grepl(sample_name, taxmap_files)]
  diamond_file <- diamond_files[grepl(sample_name, diamond_files)]

  if (length(taxmap_file) == 0) { message("No taxmap found for sample: ", sample_name); next }
  if (length(diamond_file) == 0) { message("No diamond file found for sample: ", sample_name); next }

  message("Using taxmap: ", basename(taxmap_file))
  message("Using diamond: ", basename(diamond_file))

  #############################################
  #### ---------- KRAKEN PROCESS ---------- ###
  #############################################

  kraken <- read_tsv(
    kraken_file,
    col_names = c("classified", "seq_id", "taxid", "sequence_length", "lca_mapping"),
    show_col_types = FALSE
  )

  kraken_out <- kraken %>%
    filter(classified == "C") %>% # Filter Kraken 2 file to only include classified contigs
    select(seq_id, taxid) %>% # Select only the seq_id and taxid columns
    arrange(taxid) # arrange table by taxid

  kraken_pie <- kraken_out %>%
    count(taxid, name = "group") # count the number that each taxid occurs

  taxmap <- read.csv(
    taxmap_file,
    header = TRUE,
    sep = "\t"
  ) %>%
    select(taxid, Tax_name) # Read taxonomy input file, select only the taxid and Tax_name columns

  # Merge with taxmap
  kraken_pie <- merge(kraken_pie, taxmap, by = "taxid", all.x = TRUE) # merge kraken_pie table with taxonomy map file based on matching taxid

  # Filtering
  kraken_pie <- kraken_pie %>%
    filter(group > 5) %>%
    filter(!taxid %in% c(1, 2, 131567, 9606)) # filter out any taxid that has a group value lower than 5, and further filtering based on occurrences in all kraken files (human DNA, root, etc.)

  # Dynamic threshold for top slices
  kraken_pie <- kraken_pie %>% arrange(desc(group))
  if (nrow(kraken_pie) > desired_slices) {
    threshold <- kraken_pie$group[desired_slices] # decide threshold used for plotting based on the value of "desired_slices"

    kraken_pie <- kraken_pie %>%
      mutate(Tax_name = ifelse(group < threshold, "Other", Tax_name)) %>% # Condense taxid that did not make the threshold cut into a slice labeled "Other"
      group_by(Tax_name) %>%
      summarise(group = sum(group), .groups = "drop") %>%
      arrange(desc(group))
  }

  # Force minimum size for Other
  total_count <- sum(kraken_pie$group)
  other_index <- which(kraken_pie$Tax_name == "Other")
  if (length(other_index) == 1) {
    current_other_pct <- kraken_pie$group[other_index] / total_count * 100
    remaining_pct <- 100 - other_min_pct
    remaining_total <- total_count - kraken_pie$group[other_index]
    adjust_factor <- remaining_pct / (remaining_total / total_count * 100)

    kraken_pie <- kraken_pie %>%
      mutate(group = ifelse(
        Tax_name == "Other",
        total_count * other_min_pct / 100,
        group * adjust_factor
      ))
  }

  # Compute percentages
  kraken_pie <- kraken_pie %>%
    mutate(
      pct = group / sum(group) * 100,
      pct_label = paste0(round(pct, 1), "%")
    )

  # KRAKEN PLOT
  kraken_plot <- ggplot(kraken_pie, aes(x = "", y = group, fill = Tax_name)) +
    geom_bar(stat = "identity", width = 1, color = "black") +
    coord_polar("y", start = 0) +
    theme_void() +
    labs(
      title = paste("Kraken2 Pie Chart for", sample_name), # Set title of plot to either Kraken 2/DIAMOND Pie Chart for sample#
      fill = "Species" # Set Legend title to Species
    ) +
    geom_label_repel( # Add percentage labels to the pie chart overlapping the chart itself
      aes(label = pct_label),
      position = position_stack(vjust = 0.5),
      show.legend = FALSE,
      color = "black",
      box.padding = 0.6,
      point.padding = 0.3,
      segment.color = "black",
      segment.size = 0.5
    )

  # Save Kraken plot
  ggsave(
    filename = file.path(indir, paste0("KrakenPie_", sample_name, ".png")), # Save Kraken 2 plot to a .PNG can be changed to .PDF if needed
    plot = kraken_plot,
    width = 7, height = 7, dpi = 300,
    bg = "white"
  )


  #############################################
  #### ---------- DIAMOND PROCESS ---------- ###
  #############################################
## Comments for Kraken process are the same as DIAMOND process
  diamond <- fread(
    diamond_file,
    col.names = c(
      "seq_id", "sseqid", "pident", "length", "mismatch", "gapopen",
      "qstart", "qend", "sstart", "send", "evalue", "bitscore",
      "taxid", "species", "diamond_title"
    )
  )

  diamond_out <- diamond %>%
    select(seq_id, taxid, species) %>%
    arrange(desc(taxid))

  diamond_pie <- diamond_out %>%
    count(species, name = "group") %>%
    filter(group > 50) %>%
    filter(species != "N/A") %>%
    arrange(desc(group))

  # Dynamic threshold for top slices
  if (nrow(diamond_pie) > desired_slices) {
    threshold <- diamond_pie$group[desired_slices]

    diamond_pie <- diamond_pie %>%
      mutate(species = ifelse(group < threshold, "Other", species)) %>%
      group_by(species) %>%
      summarise(group = sum(group), .groups = "drop") %>%
      arrange(desc(group))
  }

  # Force minimum size for Other
  total_count <- sum(diamond_pie$group)
  other_index <- which(diamond_pie$species == "Other")
  if (length(other_index) == 1) {
    current_other_pct <- diamond_pie$group[other_index] / total_count * 100
    remaining_pct <- 100 - other_min_pct
    remaining_total <- total_count - diamond_pie$group[other_index]
    adjust_factor <- remaining_pct / (remaining_total / total_count * 100)

    diamond_pie <- diamond_pie %>%
      mutate(group = ifelse(
        species == "Other",
        total_count * other_min_pct / 100,
        group * adjust_factor
      ))
  }

  # Compute percentages
  diamond_pie <- diamond_pie %>%
    mutate(
      pct = group / sum(group) * 100,
      pct_label = paste0(round(pct, 1), "%")
    )

  # DIAMOND PLOT
  diamond_plot <- ggplot(diamond_pie, aes(x = "", y = group, fill = species)) +
    geom_bar(stat = "identity", width = 1, color = "black") +
    coord_polar("y", start = 0) +
    theme_void() +
    labs(
      title = paste("DIAMOND Pie Chart for", sample_name),
      fill = "Species"
    ) +
    geom_label_repel(
      aes(label = pct_label),
      position = position_stack(vjust = 0.5),
      show.legend = FALSE,
      color = "black",
      box.padding = 0.6,
      point.padding = 0.3,
      segment.color = "black",
      segment.size = 0.5
    )

  # Save DIAMOND plot
  ggsave(
    filename = file.path(indir, paste0("DiamondPie_", sample_name, ".png")),
    plot = diamond_plot,
    width = 7, height = 7, dpi = 300,
    bg = "white"
  )

  message("Finished sample: ", sample_name)
}
```
</details>
