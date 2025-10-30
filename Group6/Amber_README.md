# Data download 

## 1. Download NCBI SRA (Sequence Read Archive) Toolkit, version 3.2.1
    https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit

## 2. Copy link to preferred SRA Toolkit\
   Link to AlmaLinux 64 bit architecture-  https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.2.1/sratoolkit.3.2.1-alma_linux64.tar.gz

## 3. Paste link into ARC terminal, choose preferred location
```
   curl -O https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.2.1/sratoolkit.3.2.1-alma_linux64.tar.gz
```
## 4. Unzip toolkit 
```
tar -xvf sratoolkit.3.2.1-alma_linux64.tar.gz
```
Check that download worked 
```
cd sratoolkit.3.2.1.-alma_linux64/
ls
cd bin
ls
```
## 5. Collect data from reference paper using NCBI website 
The National Center for Biotechnology Information (NCBI) [https://www.ncbi.nlm.nih.gov/] advances science and health by providing access to biomedical and genomic information.

Reference Paper: https://doi-org.ezproxy.lib.vt.edu/10.1016/j.watres.2024.121425  

Use the left menu bar and choose "SRA" and type in BioProject #  
For our reference paper, the bioproject # is **PRJNA669820**  
Press "Search"  
Then on the "Send To:" menu, click "Go"  

## 6. Get SRR File #s
There should be a new screen now, and click "Accession List"  
This will download a txt file to your computer.  
Create a .txt file in Linux and paste in all the SRR #s  

***We have chosen a subset of 10 SRR #s for the remaining downstream analysis of this pipeline.***
```
cat SRR_File.txt 
SRR12900993
SRR12901004
SRR12901025
SRR12901034
SRR12901035
SRR12901036
SRR12901037
SRR12901038
SRR12901047
SRR12901048
```
## 7. Create a .sh file to run a batch script with SRA Toolkit and SRR File #s 
Use `prefetch` and `fasterq-dump` to obtain fastq files from the SRA run accession list.  
A more detailed explanation can be foudnd here   
https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump  

Create a file to write in the batch script 
```
vim SRR_slrum  # change to desired name

```
Press "I" to begin writing in a vim file    
Press "esc" button when want to exite and use ":wq" to save and exit  
```
#!/bin/bash

#SBATCH --account=metagen    # change to desired account 

#SBATCH --partition=normal_q    # change to desired partition

#SBATCH -t 12:00:00

#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --mail-user=aaamber@vt.edu    # change to desired email address 


mkdir -p fastq_files # makes a folder for output of fastq files 

# Add your SRA Toolkit bin directory to PATH

export PATH=/projects/intro2gds/I2GDS2025/G6_AMR_ARG/Amber/SRA_Demo/sratoolkit.3.2.1-alma_linux64/bin:$PATH

for acc in $(cat SRR_File.txt); do

    prefetch "$acc"

    fasterq-dump "$acc" --split-files --threads 4 -O fastq_files

done
```
## 8. Submit your job 
```
sbatch SRR_slurm.sh  # use desired name
```
Fastq files should now be available in the folder created from the batch script. 
There will be two files (_1 and _2) per SRR # which represent the forward and reverse reads. 

# Fastqc – assess quality of raw sequencing reads 

### 1. Create Miniconda environment on the ARC and add packages into your environment 
Add Miniconda 
```
module load Miniconda3
```
Add fastqc to the environment 
```
conda create -- name fastqc
```
Install fastqc 
```
conda install bioconda::fastqc
```
To activate
```
source activate fastqc
```
### 2. Create a file to write in the fastqc batch script
```
#!/bin/bash
#SBATCH --account=prudenlab   # change to desired account 
#SBATCH --partition=normal_q   # change to desird partition 
#SBATCH -t 3:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aaamber@vt.edu   # change to desired email address 

mkdir -p /projects/intro2gds/I2GDS2025/G6_AMR_ARG/Amber/fastqc_results    # makes a directory to store the fastqc outputs 

module load Miniconda3

source activate fastqc 


cd /projects/intro2gds/I2GDS2025/G6_AMR_ARG/Amber/RefPaperSubset    # location where fastq files are located

for file in *.fastq.gz; do       # can change to .fastq if files are not zipped 

    fastqc "$file" -o /projects/intro2gds/I2GDS2025/G6_AMR_ARG/Amber/fastqc_results/  # runs fastqc on each file and puts the output in desired path

done

```
### 3. Check quality score with FastQC Report  
Each file will produce a fastqc.html and fastqc.zip  
Download the .html file to your computer and open to evaluate the FastQC Report for EACH sample  
Remember, each sample has a corresponding _1 and _2 file. This means you need to check the report for each _1 and _2 file. 

Here is a table of what the Phred score (y-axis on "Per base sequence quality" graph) represents:  
<img width="525" height="189" alt="image" src="https://github.com/user-attachments/assets/a9e8ea41-48d3-43cf-8375-53d6a9ac17d0" />

### 4. Combine files using MultiQC
An alternative to checking a fastqc report for each _1 and _2 file of each sample is creating a report that contains all files  
This is done using `multiqc` 
Multiqc will combine all fastqc files together and create a single .html file that summarizes all samples  

Install mulitqc  
```
module load Miniconda3
```
Install into environment  
```
conda install bioconda::multiqc 
```
Activate multiqc in the command line  
```
source activate multiqc
```
The input for multiqc is the directory that contains the fastqc files made from running fastqc  
Run multiqc on that folder
```
multiqc fastqc_results/    # change folder to desired name
```
Download multiqc_report.html file and open it in web browser  


