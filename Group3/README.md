# Analyzing bisulfite single cell single-cell DNA methylomes

This project aims to replicate the single-cell DNA methylation analysis described in:

Zhang, Q., Ma, S., Liu, Z. et al. (2023). Droplet-based bisulfite sequencing for high-throughput profiling of single-cell DNA methylomes. Nature Communications 14, 4672. https://doi.org/10.1038/s41467-023-40411-w

Reference paper: [download here](https://www.nature.com/articles/s41467-023-40411-w)  
These data can be downloaded [from NCBI.](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204691)

## Pipeline description
To replicate the analysis, we use Virginia Tech’s ARC Computing Cluster.
The pipeline includes the following steps:

 •	UMI-tools: extraction of unique molecular identifiers (UMIs) and generation of a whitelist of valid cell barcodes

 •	TrimGalore: adapter removal and quality trimming of raw reads

 •	idemp: demultiplexing of reads based on cell barcodes

 •	Bismark: alignment of bisulfite-treated reads to the reference genome

 •	BamTools: merging and format conversion of multiple BAM files

 •	Picard: removal of PCR duplicates and generation of final cleaned alignment files

### 00 Data Download

### 01 UMI extraction and whitelist creation (UMI-tools)

### 02 Quality and adapter trimming (TrimGalore)

### 03 Read demultiplexing (idemp)
The purpose of idemp is to demultiplex reads based on cell barcodes. After extracting UMIs and generating a whitelist with UMI-tools, idemp uses that whitelist to assign each read to its corresponding cell, creating separate FASTQ files for each cell.

**Input:**
 
 •	Trimmed FASTQ files (Trim Galore output)
 
 •	Whitelist file (UMI-tools output)
 
 •	Output directory path


**What the script does:**
 1.	Loads the idemp module on ARC or activates the environment containing idemp.
 2.	Runs idemp to assign reads to individual cells according to their barcodes.
 3.	Writes demultiplexed FASTQ files, one per cell barcode, into the output directory.


**Output:**

 •	Demultiplexed paired FASTQ files
 
 •	Log file with summary of read counts per barcode


**Note:**

Whitelist file from UMI-tools needs to be the same format expected by idemp. We need to check the quote before running as idemp can generate many small files for large datasets.


### 04 Bisulfite alignment (Bismark)
The purpose of Bismark is to align bisulfite-treated reads to the reference genome. Bismark performs both alignment and methylation context reporting, ensuring that C-->T (and G-->A) conversions introduced by bisulfite treatment are properly handled during mapping.

**Input:**

 •	Demultiplexed FASTQ files (idemp output)

 •	Reference genome directory

 •	Output directory path


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


### 05 File merging and format conversion (BamTools)

### 06 Deduplication and QC (Picard)
