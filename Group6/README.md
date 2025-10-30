# Classifying Antimicrobial Resistance Genes from Short-Read Metagenomic Data
Data and required reference files are located here:
```
/projects/intro2gds/I2GDS2025/G6_AMR_ARG/ForReviewers
```

# 1. Data download 

## 1.1. Download NCBI SRA (Sequence Read Archive) Toolkit, version 3.2.1
    https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit

## 1.2. Copy link to preferred SRA Toolkit\
   Link to AlmaLinux 64 bit architecture-  https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.2.1/sratoolkit.3.2.1-alma_linux64.tar.gz

## 1.3. Paste link into ARC terminal, choose preferred location
```
   curl -O https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.2.1/sratoolkit.3.2.1-alma_linux64.tar.gz
```
## 1.4. Unzip toolkit 
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
## 1.5. Collect data from reference paper using NCBI website 
The National Center for Biotechnology Information (NCBI) [https://www.ncbi.nlm.nih.gov/] advances science and health by providing access to biomedical and genomic information.

Reference Paper: https://doi-org.ezproxy.lib.vt.edu/10.1016/j.watres.2024.121425  

Use the left menu bar and choose "SRA" and type in BioProject #  
For our reference paper, the bioproject # is **PRJNA669820**  
Press "Search"  
Then on the "Send To:" menu, click "Go"  

## 1.6. Get SRR File #s
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
## 1.7. Create a .sh file to run a batch script with SRA Toolkit and SRR File #s 
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
## 1.8. Submit your job 
```
sbatch SRR_slurm.sh  # use desired name
```
Fastq files should now be available in the folder created from the batch script. 
There will be two files (_1 and _2) per SRR # which represent the forward and reverse reads. 

# 2. Fastqc and Multiqc – assess quality of raw sequencing reads 

### 2.1. Create Miniconda environment on the ARC and add packages into your environment 
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
### 2.2. Create a file to write in the fastqc batch script
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
### 2.3. Check quality score with FastQC Report  
Each file will produce a fastqc.html and fastqc.zip  
Download the .html file to your computer and open to evaluate the FastQC Report for EACH sample  
Remember, each sample has a corresponding _1 and _2 file. This means you need to check the report for each _1 and _2 file. 

Here is a table of what the Phred score (y-axis on "Per base sequence quality" graph) represents:  
<img width="525" height="189" alt="image" src="https://github.com/user-attachments/assets/a9e8ea41-48d3-43cf-8375-53d6a9ac17d0" />

### 2.4. Combine files using MultiQC
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


# 3. Trimmomatic (Cleaning Data)
### 3.1 Introduction
Trimmomatic is a tool used to clean data of low quality reads, and is designed specifically for Illuminia paired-end or single-end sequencing data.\
For more info, visit the [Trimmomatic Website](http://www.usadellab.org/cms/?page=trimmomatic)

### 3.2 Creating Environment
First, create a new environment on ARC. Naming it "trimmomatic" matches the script below:
```
conda create -n trimmomatic
```

After creating the environment, activate it:
```
source activate trimmomatic
```

### 3.3 Downloading Trimmomatic
Next, download trimmomatic:
```
conda install bioconda::trimmomatic
```

### 3.4 Inputs Required
- **Input data**: Paired end FASTQ/FASTQ.gz files
- **Illuminia Specific Adapter File**: Adapter file (TruSeq3-PE.fa)

### 3.5 Trimmomatic Script

### 3.5.1 Explanation of Trimmomatic.sh
**Make a list of sample names**.\
This will be used to cycle through your samples using a for loop. \
For example, if we have the file "SRR12900993_1.fastq.gz", then the code below will produce "SRR12900993".
```
sample=$(ls *.fastq* | awk -F/ '{gsub(/_[12].fastq.gz/, "", $NF); print $NF}' | sort | uniq)
```

**Trimmomatic Command** \
This is the paired-end version, hence the "PE". "-Xmx16G" is added to give the command addtional memory. "phred33" is the current way Illumina scores the quality of the reads:
```
trimmomatic PE -Xmx16G -threads 8 -phred33
```

**Adapter Parameters** \
You can customize how stringent you want your trimming of adapters to be. The code uses the default parameters 
- ILLUMINACLIP:fastaWithAdaptersEtc:seed mismatches:palindrome clip threshold:simple clip threshold
	- fastaWithAdaptersEtc: specifies the path to a fasta file containing all the adapters. This is TruSeq3-PE.fa.
	- seedMismatches: specifies the maximum mismatch count which will still allow a full match to be performed
	- palindromeClipThreshold: specifies how accurate the match between the two 'adapter ligated' reads must be for PE 			palindrome read alignment.
	- simpleClipThreshold: specifies how accurate the match between any adapter etc. sequence must be against a read.

```
ILLUMINACLIP:${ADAPTERS}:2:30:10
```

Description from [Trimmomatic Website](http://www.usadellab.org/cms/?page=trimmomatic)

**Trimmomatic parameters**\
You can customize how stringent you want your cutting to be. The code uses the default parameters
- LEADING: Cut bases off the start of a read, if below a threshold quality
- TRAILING: Cut bases off the end of a read, if below a threshold quality
- SLIDINGWINDOW:windowSize:requiredQuality, will cut if below threshold
	- windowSize: specifies the number of bases to average across
	- requiredQuality: specifies the average quality required.
- MINLEN: Drop the read if it is below a specified length

```
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

Descriptions from [Trimmomatic Website](http://www.usadellab.org/cms/?page=trimmomatic)

### 3.5.2 Script
With the 10 input files, this should take ~1 hour to run. Called Trimmomatic_G6.sh in the "code-files" folder.
<details>
<summary> Trimmomatic.sh</summary>

```
#!/bin/bash
#SBATCH --account=YourAccount #Change
#SBATCH --partition=normal_q
#SBATCH --mem=16G
#SBATCH -t 30:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=username@vt.edu #Change

module load Miniconda3

source activate trimmomatic

cd /projects/intro2gds/I2GDS2025/G6_AMR_ARG/SampleSubset #Change to location of input files

OUTPUT_DIR=/projects/intro2gds/I2GDS2025/your/output/folder #Change
ADAPTERS=TruSeq3-PE.fa

sample=$(ls *.fastq* | awk -F/ '{gsub(/_[12].fastq.gz/, "", $NF); print $NF}' | sort | uniq)

for samples in $sample; do
	echo ${samples}
	trimmomatic PE -Xmx16G -threads 8 -phred33 \
		 "${samples}"_1.fastq.gz "${samples}"_2.fastq.gz \
		 "${OUTPUT_DIR}"/"${samples}"_1_paired.fastq.gz "${OUTPUT_DIR}"/"${samples}"_1_unpaired.fastq.gz \
		 "${OUTPUT_DIR}"/"${samples}"_2_paired.fastq.gz "${OUTPUT_DIR}"/"${samples}"_2_unpaired.fastq.gz \
		 ILLUMINACLIP:${ADAPTERS}:2:30:10 \
 		 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

		echo "finished ${samples}"
done

```
</details>

### 3.6 Output Explanation
Each foward and reverse file will give 2 output files:
- Paired: These reads survived trimming and their "partner" read survived trimming
- Unpaired: These reads survived trimming, but their "partner" read did not

Typically, analysis only moves forward with the paired reads. However, unpaired reads can still be used for some future analysis.

# 4. bbduk (Decontaminating Data)
### 4.1 Introduction
bbduk is a multifaceted tool which has the capabilities to do trimming, read quality filtering, and more; but we will only use it for decontamination. 
Decontamination will remove sequecnes that are known to be human, dog, cat, or mouse. 

### 4.2 Creating Environment
First, create a new environment on ARC. Naming it "bbduk" matches the script below:
```
conda create -n bbduk
```

After creating the environment, activate it:
```
source activate bbduk
```

### 4.3 Downloading bbduk
Next, download bbduk, which is part of the bbmap package:
```
conda install bioconda::bbmap
```

### 4.4 Inputs Required
- **Input Data**: This will be the paired sequences from Trimmomatic
- **Reference Database**: A FATSA file containing sequences known to be human, cat, dog, or mouse. This is "merged_ref_5081444234059102403.fa.gz"


### 4.5 bbduk Script

### 4.5.1 Explaination of bbduk.sh

**bbduk command**\
preallocate=t means we will preallocate memory (t = true), and we will preallocate 180G
```
bbduk.sh preallocate=t -Xmx180G
```

If a read is shorter than 51 after decontamination, then we will throw it out
```
minlength=51
```
The reference database (merged_ref_5081444234059102403.fa.gz) needs to be explicitly referenced
```
ref='/chnage/to/where/merged_ref_5081444234059102403.fa.gz/is/located'
```

Pretty simple! Everything else should be self explanatory
### 4.5.2 bbduk script
With 10 samples, should take ~13 hours. Called "bbduk_G6.sh" in the "code-files" folder.
<details>
<summary> bbduk.sh</summary>

```
#!/bin/bash

#SBATCH --account=YourAccount #Change
#SBATCH --partition=normal_q
#SBATCH --mem=200G
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=username@vt.edu #Change

module load Miniconda3

source activate bbduk

cd /projects/intro2gds/I2GDS2025/location/of/trimmomatic/output #Change

sample=$(ls *_[12]_paired.fastq.gz | awk -F/ '{gsub(/_[12]_paired.fastq.gz/, "", $NF); print $NF}' | sort | uniq)

for one_sample in $sample; do
    echo ${one_sample}

    bbduk.sh preallocate=t -Xmx180G \
        minlength=51 \
        in="${one_sample}"_1_paired.fastq.gz \
        in2="${one_sample}"_2_paired.fastq.gz \
        out=/change/to/output/directory/"${one_sample}"_1_decontam.fastq.gz \ #Change
        out2=/change/to/output/directory/"${one_sample}"_2_decontam.fastq.gz \ #Change
        stats=/change/to/output/directory/"${one_sample}"_decontamination_stats.txt \ #Change
        ref='/chnage/to/where/merged_ref_5081444234059102403.fa.gz/is/located' #Change


done

```



  
</details>

### 4.6 Output Explaination
For each input file, there will be one output FASTQ file, and for each pair of FASTQ files there will be one txt file. This txt file shows how many of each sequence in the reference FASTA file are removed (should be ~0.5-1% of total samples).


# vesearch


# DIAMOND


# Pipeline
<img width="593" height="778" alt="image" src="https://github.com/user-attachments/assets/6f9dc70f-19f5-4431-bddc-5a1dc2c58afc" />
