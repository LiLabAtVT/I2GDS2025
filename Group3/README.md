# Analyzing bisulfite single cell single-cell DNA methylomes

This project aims to replicate the single-cell DNA methylation analysis described in:

Zhang, Q., Ma, S., Liu, Z. et al. (2023). Droplet-based bisulfite sequencing for high-throughput profiling of single-cell DNA methylomes. Nature Communications 14, 4672. https://doi.org/10.1038/s41467-023-40411-w

Reference paper: [download here](https://www.nature.com/articles/s41467-023-40411-w)  
These data can be downloaded [from NCBI.](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204691)

## Pipeline description
To replicate the analysis, we use Virginia Tech’s ARC Computing Cluster.
The pipeline includes the following steps:

 1.	UMI-tools: extraction of unique molecular identifiers (UMIs) and generation of a whitelist of valid cell barcodes
 2.	TrimGalore: adapter removal and quality trimming of raw reads
 3.	idemp: demultiplexing of reads based on cell barcodes
 4.	Bismark: alignment of bisulfite-treated reads to the reference genome
 5.	BamTools: merging and format conversion of multiple BAM files
 6.	Picard: removal of PCR duplicates and generation of final cleaned alignment files
<img width="3076" height="1815" alt="image" src="https://github.com/user-attachments/assets/929b9e6a-11d3-4dd9-bf86-7e735d83b8fc" />

### 00 Environment setup and tool installation
Before running the scripts, follow these steps to create and activate the environment and install the necessary tools.

**Environment creation:** create a new Conda environment and assign it any name. In this example, we name it “bisulfite” and specify Python 3.10 to ensure compatibility with the tools.
```
conda create -n bisulfite -c conda-forge -c bioconda python=3.9
```
**Environment activation:** activate the newly created environment so that all installed packages and tools will be available for use.
```
source activate bisulfite
```
**Tool installation:**: install the software packages needed for the first two steps of the pipeline: UMI-tools (for UMI extraction) and Trim Galore (for adapter and quality trimming). Channel priority settings are adjusted to avoid dependency conflicts during installation.
```
# Install UMI-tools
python -m pip install --no- cache-dir "umi_tools==1.1.3"

# Temporarily set flexible channel priority
conda config --set channel_priority flexible

# Install Trim Galore
conda install -c bioconda trim-galore -y

# Reset channel priority to strict
conda config --set channel_priority strict
```

### 01 Data Download and subset extraction (SRAtools)

SRAtools (Sequence Read Archive Toolkit) is a set of command-line tools provided by NCBI (National Center for Biotechnology Information) for downloading, viewing, and converting sequencing data from the SRA (Sequence Read Archive) database. First, we locate the file we need on the NCBI website, and then use the `prefetch` tool to download it.
For example, if the file we need is SRR19391270, the specific command is:
```
prefetch	SRR19391270
```
After downloading the file, we will obtain a file named SRR19391270.sra. Next, we need to use the `fasterq-dump` program to convert it into the FASTQ format, which can be processed further. The script is as follows:
```
fasterq-dump --split-files -e 8 SRR19391270.sra -O .
```
This command can split the SRA file into two paired-end sequencing files for subsequent processing.

Because the original file is very large, the `Seqtk` tool was used to extract a portion of the data for scripesting. The command for data extraction is as follows.
```
seqtk sample -s77 SRR19391270_1.fastq 0.01 > sub_R1.fastq
seqtk sample -s77 SRR19391270_2.fastq 0.01 > sub_R2.fastq
```
The s77 represents random seed, which is used to randomly extract the reads. 0.01 means 1% of data will be extracted, and the extracted data will be output to sub_R1.fastq and sub_R2.fastq.


### 02 UMI extraction and whitelist creation (UMI-tools)
#### Introduction
`UMI-tools` is a software to process sequencing data containing Unique Molecular Identifiers (UMIs). It helps remove PCR duplicates and accurately count unique molecules, improving the reliability of single-cell and RNA-seq analyses.
#### Explanation
In the experiment, we obtained many cells, but we only wanted to analyze those with the highest quality. Therefore, `UMI-tools` was used to extract the cells with the top 5,000 highest read counts based on the barcode. First, extract the whitelist from the original file. The specific command is:
```
#!/bin/bash
#SBATCH --job-name=extract_trim
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=extract_trim_%j.out
#SBATCH --error=extract_trim_%j.err

sub_R1="sub_R1.fastq"
sub_R2="sub_R2.fastq"

umi_tools whitelist \
  --stdin "$sub_R2" \
  --extract-method=regex \
  --bc-pattern='(?P<discard_1>TAAGTAGAAGATGGTATATGAGAT){s<=4}(?P<cell_1>.{15})(?P<umi_1>.{0})' \
  --subset-reads=10000000 \
  --error-correct-threshold=1 \
  --set-cell-number=100 \
  --log2stderr > whitelist_top100.txt
```
Here, up to 10,000,000 reads are analyzed with the threshold set to 1. The fixed sequence TAAGTAGAAGATGGTATATGAGAT is used to locate the barcode region, and the 15 bases immediately following this sequence are extracted as the barcode. Based on the read counts, the top 5,000 most frequent barcodes are selected, and a whitelist file is generated accordingly. Then, the fastq file will be extracted based on the whitelist. The script is as follows:
```
umi_tools extract \
  --extract-method=regex \
  --bc-pattern='(?P<discard_1>TAAGTAGAAGATGGTATATGAGAT){s<=4}(?P<cell_1>.{15})(?P<umi_1>.{0})' \
  --error-correct-cell \
  --stdin "$sub_R2" --stdout "${sub_R2%.fastq}_extracted.fastq" \
  --read2-in "$sub_R1" --read2-out "${sub_R1%.fastq}_extracted.fastq" \
  --filter-cell-barcode --whitelist=whitelist_top100.txt
  ```
Here, the barcodes are extracted again based on the fixed position and output to the extracted.fastq file. These extracted barcodes are then matched against the whitelist generated in the previous step.

#### Expected Output
This step produces two `.fastq` files and a whitelist text file that includes barcodes. A log file accompanies the text file.

### 03 Adapter trimming (TrimGalore)
#### Introduction
`Trim Galore` is used to remove sequencing adapter sequences, trim low-quality bases, and filter out short reads. Here, we use it to perform trimming and cleaning on the extracted FASTQ files.v
#### Explanation

The script is as follows.
  ```
trim_galore --paired --quality 20 --length 20 "sub_R1"_extracted.fastq "sub_R2"_extracted.fastq
  ```
This script removes bases with a quality score below Q20, discards reads shorter than 20 bp, and automatically detects and synchronizes trimming of paired-end reads.

#### Expected Output

This step should produce two `.fq` files with accompanying report `.txt` files. An `err` and `log` file for `Trim_Galore` should also be present.

### 04 Read demultiplexing (idemp)
The purpose of idemp is to demultiplex reads based on cell barcodes. After extracting UMIs and generating a whitelist with UMI-tools, idemp uses that whitelist to assign each read to its corresponding cell, creating separate FASTQ files for each cell.

**Input:**
 
 •	Trimmed FASTQ files (Trim Galore output)
 
 •	Whitelist file (UMI-tools output)
 
 •	Output directory path

**Script:**
```
OUTDIR="03_idemp"
mkdir -p "$OUTDIR"

# input
WL="extract_trim_out/whitelist_top100.txt"
R1="extract_trim_out/sub_R1_extracted_val_1.fq"
R2="extract_trim_out/sub_R2_extracted_val_2.fq"

# 1. barcode_r.txt
awk '{print $1, NR}' "$WL" > "$OUTDIR/barcode_r.txt"

# 2. index.fastq
awk 'NR%4==1{
  if (match($0, /_([ACGTN]+)_/, m)) bc=m[1]; else bc="N";
  q=""; for(i=1;i<=length(bc);i++) q=q"F";
  print "@idx" ++c; print bc; print "+"; print q;
}' "$R1" > "$OUTDIR/index.fastq"

# 3. run idemp
idemp \
  -b "$OUTDIR/barcode_r.txt" \
  -I1 "$OUTDIR/index.fastq" \
  -R1 "$R1" \
  -R2 "$R2" \
  -m 0 \
  -o "$OUTDIR"
```

**What the script does:**
 1.	Loads the idemp module on ARC or activates the environment containing idemp.
 2.	Runs idemp to assign reads to individual cells according to their barcodes.
 3.	Writes demultiplexed FASTQ files, one per cell barcode, into the output directory.


**Output:**

 •	Demultiplexed paired FASTQ files
 
 •	Log file with summary of read counts per barcode


**Note:**

Whitelist file from UMI-tools needs to be the same format expected by idemp. We need to check the quote before running as idemp can generate many small files for large datasets.


### 05 Bisulfite alignment (Bismark)
The purpose of Bismark is to align bisulfite-treated reads to the reference genome. Bismark performs both alignment and methylation context reporting, ensuring that C-->T (and G-->A) conversions introduced by bisulfite treatment are properly handled during mapping.

**Input:**

 •	Demultiplexed FASTQ files (idemp output)

 •	Reference genome directory

 •	Output directory path

**Script:**
```
GENOME_DIR="/projects/intro2gds/I2GDS2025/G3_SingleCell/Rui/reference/Bisulfite_Genome_human"
INDIR="03_idemp"
OUTDIR="04_bismark_align"
THREADS="${SLURM_CPUS_PER_TASK:-8}"


mkdir -p "$OUTDIR"


# ---- Align R1 reads (PBAT mode) ----
for f in ${INDIR}/*_R1_*.fastq.gz; do
  base=$(basename "$f" .fastq.gz)
  echo "[$(date)] Aligning (PBAT) $base ..."
  bismark --genome "$GENOME_DIR" --pbat -p "$THREADS" \
          -o "$OUTDIR" "$f"
done

# ---- Align R2 reads (default mode) ----
for f in ${INDIR}/*_R2_*.fastq.gz; do
  base=$(basename "$f" .fastq.gz)
  echo "[$(date)] Aligning (default) $base ..."
  bismark --genome "$GENOME_DIR" -p "$THREADS" \
          -o "$OUTDIR" "$f"
done
```

**What the script does:**
 1.	Loads the Bismark and Bowtie2 modules.
 2.	Specifies the reference genome path.
 3.	Loops through all FASTQ files in the input directory and runs Bismark alignment.
 4.	Converts the resulting SAM files to BAM format and sorts them using samtools.
 5.	Outputs BAM files and alignment reports.


**Output:**

 •	Sorted BAM files for each demultiplexed cell or sample
 
 •	Bismark alignment summary reports


**Note:**

The reference genome must be bisulfite-converted and indexed before running alignment. We will need to specify one input file for single-end reads or include -1 and -2 flags with corresponding files for aligning paired-end reads.


### 06 File merging and format conversion (BamTools)

The previous step provide sorted BAM files and that must be merged. **BamTools** is a collection of tools for working with BAM files. More information about BamTools can be found in this [Toolkit Tutorial document](https://raw.githubusercontent.com/wiki/pezmaster31/bamtools/Tutorial_Toolkit_BamTools-1.0.pdf). A [repository for BamTools is also on GitHub](https://github.com/pezmaster31/bamtools/tree/master). The tool used here is _merge_. It is also possible to use additional BamTools tools for filtering and sorting if required.

BamTools is available on the VT ARC listed in the [Table of Software](https://www.docs.arc.vt.edu/software/01table.html). The version loaded on ARC is 2.5.2 (_as of October 2025_).

To complete the merge, use BamTools in a SLURM script. The main commands for BamTools in the script include specifying the directory containing the input files and the directory for the merged output files.

In the script, the input and output directories are specified:
```
INDIR="04_bismark_align"
OUTDIR="05_bamtools"
```
The script also includes an array in order to cycle through each of the cell's files using the `cell_ids.txt` file to provide the numbers.

The specific merge code is:
```
INDIR="/projects/intro2gds/I2GDS2025/TestData_LinuxPeerEval/G3_TestData/04_bismark_align"
OUTDIR="/projects/intro2gds/I2GDS2025/TestData_LinuxPeerEval/G3_TestData/05_bamtools"

mkdir -p "$OUTDIR"

# Merge + Sort
for R1 in ${INDIR}/sub_R1_extracted_val_1.fq_*_bismark_bt2.bam
do
  ID=$(basename "$R1" | sed 's/.*\.fq_//; s/_bismark_bt2.bam//')
  R2="${INDIR}/sub_R2_extracted_val_2.fq_${ID}_bismark_bt2.bam"
  MERGED="${OUTDIR}/${ID}.merged.bam"
  SORTED="${OUTDIR}/${ID}_sorted.bam"
  
  # Merge R1 + R2
  bamtools merge -in "$R1" -in "$R2" -out "$MERGED"

  # Sort merged BAM
  samtools sort "$MERGED" -o "$SORTED"

  # Delete unsorted intermediate file to save space
  rm -f "$MERGED"

  echo "[$(date)] Done ${ID} -> ${SORTED}"
done
```
`merge` is the command in BamTools. The two `-in` statements indicated the files to be combined into the merged file in `-out`.

### 06 Deduplication (Picard)

With single cell genomic analysis, PCR amplification of sequences is necessary to allow sequencing but can also lead to duplication of sequences. Removing PCR duplicates is a cleaning procedue. For this pipeline, Picard is used. Picard is a Java based set of tools for manipulating sequencing data in BAM files as well as other file formats. Information about Picard can be found [at the Broad Institute](https://broadinstitute.github.io/picard/index.html). Resources can also be found [in a GitHub repository](https://github.com/broadinstitute/picard). For this pipeline, we used the `MarkDuplicates` tool and `REMOVE_DUPLICATES=true` in the SLURM script.

The script starts with loading `picard` and `samtools`. SAMTools will be used for sorting and indexing the BAM files. Both `picard` and `samtools` are installed on the VT ARC listed in the [Table of Software](https://www.docs.arc.vt.edu/software/01table.html). Picard is version 3.3.0 and SAMTools is version 1.21 (_as of October 2025_).
```
module load picard
module load samtools
```
A SLURM array is used, as in the previous step, to cycle through all 1000 cells in the single cell dataset, each with its own BAM file.
```
BAM_FILE=$(ls $IN_DIR/*.merged.bam | sed -n "${SLURM_ARRAY_TASK_ID}p")
```
Subsequently, `samtools` is used to sort the the BAM files according to genomic coordinates and outputs a new sorted file.
```
samtools sort -@4 -o $OUT_DIR/${BASE}.sorted.bam "$BAM_FILE"
```
Then `picard` is called to remove the duplicates. `I=` indicates where the inputfiles are located with `O=` providing the directory and new name of the files with duplicates removed. `M=` sets up a metric file.
```
INDIR="/projects/intro2gds/I2GDS2025/TestData_LinuxPeerEval/G3_TestData/05_bamtools"
OUTDIR="/projects/intro2gds/I2GDS2025/TestData_LinuxPeerEval/G3_TestData/06_picard"

mkdir -p "$OUTDIR"

for BAM in ${INDIR}/*_sorted.bam
do
  ID=$(basename "$BAM" _sorted.bam)
  OUT_BAM="${OUTDIR}/${ID}_dedup.bam"
  METRICS="${OUTDIR}/${ID}_metrics.txt"

  echo "[$(date)] Processing ${ID} ..."

  picard MarkDuplicates \
    I="$BAM" \
    O="$OUT_BAM" \
    M="$METRICS" \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=SILENT

  echo "[$(date)] Done ${ID} -> ${OUT_BAM}"
done
```
And finally, `samtools` is used to create an index file. This index file can be used to accelerate finding specific regions of the BAM files allowing faster processing in downstream steps.
```
samtools index $OUT_DIR/${BASE}.dedup.bam
```
