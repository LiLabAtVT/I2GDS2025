
# Metagenomic Preprocessing and Taxonomic Classification
Our pipeline aims to replicate the preprocessing and taxonomic classification steps from
"Short- and long-read metagenomics of urban and rural South African gut microbiomes reveal 
a transitional composition and undescribed taxa" by Tamburini et al.

Link to paper: [https://www.nature.com/articles/s41467-021-27917-x](https://www.nature.com/articles/s41467-021-27917-x)

## 00. Download data from SRA


## 01. TrimGalore


## 02. Super Deduper


## 03. BWA

### 03a. Download and index human reference genome

In order to run BWA, we have to have a human reference genome to compare to. 
For our analysis, we downloaded the same human reference genome used in the paper
we are trying to replicate. 

There are a few different steps to download the genome and convert 
it to a usable format for BWA.

1. Download reference genome
```
mkdir -p /projects/intro2gds/I2GDS2025/G5_MG_AMR/03a_Human_ref/ # change to preferred directory path
cd /projects/intro2gds/I2GDS2025/G5_MG_AMR/03a_Human_ref/
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit;
```
2. Convert to fasta format
```
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa;
chmod +x twoBitToFa;
./twoBitToFa hg19.2bit hg19.fa;
```
3. Create bowtie index for BWA to use for comparisons
This final step of obtaining the reference genome was completed with a bash script: `03a_index_human_ref.sh`
This script loads the ARC module used to access BWA commands, set the current directory to the folder
where you downloaded the reference genome, and runs the command used to index the reference genome. 

### Run BWA Analysis

BWA analysis was run using a batch script: `03_BWA.sh`

**Inputs needed:**
* Input directory path: where all the sequence files going into BWA are located
* Output directory path: Where all the sequence files coming out of BWA should be saved
* Path to human reference genome

**What the script does:**
1. Loads the BWA module on ARC
2. Defines the input directory where the samples going into BWA are located (in this case, the output sequences
from super deduper) and the output directory where the files with human mapped sequences removed will be saved.
3. Changes the working directory to the input directory, loops through all sequence files and runs them through 
BWA, and then uses samtools to process the output from BWA and save the unmapped sequences as fastq files.

**Notes:**
The output from BWA for each sample is three files instead of two like the input. This is because sometimes
when mapping sequences to the reference genome only one half of a pair will map to the reference. In this case, only 
the one sequence that mapped is removed and the now unpaired sequence is saved to the singletons file. Thus, each
sample has two paired-read files and a singletons file as output for BWA.


## 04. Final quality check with FastQC


## 05. Taxonomic Classification with Kraken2

### 05a. Build (download) reference databases for Kraken2