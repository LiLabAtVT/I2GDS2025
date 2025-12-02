# Classifying Antimicrobial Resistance Genes from Short-Read Metagenomic Data
Raw data is located here:
```
/projects/intro2gds/I2GDS2025/G6_AMR_ARG/ForReviewers/Raw Sequences
```
Bash scripts written throughout are also located in our "code-files" GitHub folder for download

### Reviewers: Only do Steps 4-6! Input data for step 4 is the cleaned Trimmomatic data, located:
```
/projects/intro2gds/I2GDS2025/G6_AMR_ARG/ForReviewers/trimmomatic_output_key
```
If you have issues accessing the data, email dceglio@vt.edu
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
Press "esc" button to exit and use ":wq" to save and exit or ":q" to exit and NOT save 
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
- **Illuminia Specific Adapter File**: Adapter file (TruSeq3-PE.fa, located here:/projects/intro2gds/I2GDS2025/G6_AMR_ARG/ForReviewers)

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

If you can't get it to work, the output files are located here:
```
/projects/intro2gds/I2GDS2025/G6_AMR_ARG/ForReviewers/trimmomatic_output_key
```

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
- **Reference Database**: A FATSA file containing sequences known to be human, cat, dog, or mouse. This is "merged_ref_5081444234059102403.fa.gz" (located:/projects/intro2gds/I2GDS2025/G6_AMR_ARG/ForReviewers)


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
ref='/change/to/where/merged_ref_5081444234059102403.fa.gz/is/located'
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

### 4.6 Output Explanation
For each input file, there will be one output FASTQ file, and for each pair of FASTQ files there will be one txt file. This txt file shows how many of each sequence in the reference FASTA file are removed (should be ~0.5-1% of total samples).

If you can't get bbduk to work, the output files are located here:
```
/projects/intro2gds/I2GDS2025/G6_AMR_ARG/ForReviewers/bbduk_output_key
```

# 5. vsearch
## Merging your trimmed and decontaminated reads

### 5.1 Introduction
Use your bbduk output files (_Sample Name_\_[1 or 2]\_decontam.fastq.gz) and merge the forward (\_1\_) and reverse (\_2\_) reads for each sample. Your input should be the bbduk output files you generated from the previous step, or the contingency bbduk output files provided.

### 5.2 Creating your environment for vsearch
Create an environment where you'll install vsearch, activate said environment
```
module load Miniconda3/24.7.1-0
conda create vsearchenv
source activate vsearchenv
```

You should now see (vsearchenv) \[yourPID@clustername directoryname\]$, indicating that you've entered the vsearchenv environment

### 5.3 Installing vsearch
Install vsearch in environment
```
conda install bioconda::vsearch
```
Press "y" when prompted to complete installation

### 5.4 Running vsearch
Use code below to access bbduk output files and merge forward and reverse reads. Be sure to replace the PID and directory path for your input files. Change allocation (--account) if necessary.
<details>
  <summary>runningvsearch.sh</summary>
  
``` 
#!/bin/bash
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yourPID@vt.edu

module load Miniconda3

source activate vsearchenv

cd #insert path to directory where all bbduk output files are located

samples=$(ls *_[12]_decontam.fastq.gz | awk -F/ '{gsub(/_[12]_decontam.fastq.gz/, "", $NF); print $NF}' | sort | uniq)

for sample in $samples; do

 vsearch --fastq_mergepairs ${sample}_1_decontam.fastq.gz --reverse ${sample}_2_decontam.fastq.gz --fastaout ${sample}_merged.fastq --fastaout_notmerged_fwd ${sample}_unmerged_forward.fastq --fastaout_notmerged_rev ${sample}_unmerged_reverse.fastq --threads 128
cat ${sample}_merged.fastq ${sample}_unmerged_forward.fastq ${sample}_unmerged_reverse.fastq > ${sample}_clean_merged.fastq
rm ${sample}_merged.fastq ${sample}_unmerged_forward.fastq ${sample}_unmerged_reverse.fastq

done  
```
</details>

### 5.5 Troubleshooting and Code/Output Explanation
<ins>Troubleshooting</ins>: If you receive an "Out of Memory" error after submitting the slurm job, try decreasing threads to 10 and adding\
#SBATCH --cpus-per-task=16

<ins>Code/Output Explanation</ins>: This for loop goes through each sample's forward and reverse read and merges them. Due to unequal read lengths, some of the forward and reverse reads will not be merged, so in addition to the "merged" output, there is also an "unmerged forward" and "unmerged reverse" output. The merged, unmerged forward, and unmerged reverse reads for a sample are concatenated to output a "clean merged" file. There is no longer a use for the sample's merged, unmerged forward, and unmerged reverse reads, so they are removed from the directory.

If there's still some trouble with vsearch, contingency output files are located here:
```
/projects/intro2gds/I2GDS2025/G6_AMR_ARG/ForReviewers/vsearch_output_key
```


# 6. DIAMOND
## Using DIAMOND to annotate your reads and identify ARGs

### 6.1 Introduction
You will be running each sample read against a database of known ARG sequences to identify the ARGs present in each sample. Your input should be the vsearch output files you generated from the previous step, or the contingency vsearch output files provided.

### 6.2 Creating your environment for DIAMOND
Create an environment where you'll install DIAMOND, activate said environment
```
module load Miniconda3/24.7.1-0
conda create diamondenv
source activate diamondenv
```
You should now see (diamondenv) \[yourPID@clustername directoryname\]$, indicating that you've entered the diamondenv environment

### 6.3 Installing DIAMOND
Install DIAMOND in environment
```
conda install bioconda::diamond
```
Press "y" when prompted to complete installation

### 6.4 Downloading reference database
Download your reference database for ARG identification

For this data, we utilitze the Comprehensive Antibiotic Resistance Database (CARD), which is a bioinformatic database of resistance genes, their products, and associated phenotypes.

Go to this link [https://card.mcmaster.ca/](url), click Download, and scroll to the section that says "Download CARD Data" (NOT "Download CARD-R Resistomes, etc.)
The reference paper that collected this data annotated it using CARD v3.0.8, so we will do the same. In the "Download CARD Data" section, click "More data downloads..." and scroll down until you find version 3.0.8. Download this compressed folder and extract all files. Then, upload the "protein_fasta_protein_homolog_model.fasta" file to your working directory.

This is file your reference database. 

### 6.5 Converting reference file
Convert your reference file into a DMND file

The DIAMOND tool that annotates each sample requires a DMND file as the reference input. Right now your file is a fasta file. Use the following line of code to convert this file type from fasta to DMND (no need to create a Bash script)
  
``` 
diamond makedb --in protein_fasta_protein_homolog_model.fasta --db protein_fasta_protein_homolog_model.dmnd
```

### 6.6 Running DIAMOND
Use code below to access vsearch output files and annotate your reads against the reference file. Be sure to replace the PID and directory path for your input files. Change allocation (--account) if necessary.

<details>
  <summary>runningdiamond.sh</summary>
  
``` 
#!/bin/bash
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yourPID@vt.edu
#SBATCH --cpus-per-task=16

module load Miniconda3

source activate diamondenv

ARGDB='/projects/intro2gds/I2GDS2025/G6_AMR_ARG/McKenna/protein_fasta_protein_homolog_model.dmnd' #change path to where your converted reference DMND file is located

cd # insert path to directory where all vsearch output files are located

samples=$(ls *_clean_merged.fastq | awk '{gsub(/_clean_merged.fastq/,"",$0); print}' | sort | uniq)

for sample in $samples
do
     	diamond blastx -e 1e-10 --id 80 -k 1 --threads 10 -d $ARGDB -q ${sample}_clean_merged.fastq -o ${sample}_arg_full.csv --outfmt 6
done
```
</details>

### 6.7 Code and Output Explanations
<ins>Code Explanation</ins>: This for loop goes through each merged read and annotates it against your reference database, outputing a csv file for each sample. These files will be in your working directory.

<ins>Output Explanation</ins>: Your output will be a csv file with completed annotations of each sample, including the ARGs present. For the full explanation of each column in the csv file, reference this link [https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options](url) and scroll to "Output options".

If there's any trouble with DIAMOND annotation, example output files are located here:
```
/projects/intro2gds/I2GDS2025/G6_AMR_ARG/ForReviewers/diamond_output_key
```


# Pipeline
<img width="593" height="778" alt="image" src="https://github.com/user-attachments/assets/6f9dc70f-19f5-4431-bddc-5a1dc2c58afc" />


# 7. Working with data in R to create a stacked bar plot

### 7.1 Turning DIAMOND outputs into a workable file
DIAMOND outputs are a little tricky, as it does not compile the number of times a gene is seen in a sample. Therefore, we have to do it ourselves.

First, run DIAMOND again, this time annotating against an rpob reference database. The necessary reference file (RpoB.dmnd) is located in code-files.

<details>
  <summary>diamondrpob.sh</summary>
  
``` 
#!/bin/bash
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yourPID@vt.edu
#SBATCH --cpus-per-task=16

module load Miniconda3

source activate diamond_and_vsearch

RPOB_DB='/projects/ciwars/databases/RpoB.dmnd' #replace with location of your downloaded RpoB.dmnd file

cd #insert path to directory where all vsearch output files are located

samples=$(ls *_clean_merged.fastq | awk '{gsub(/_clean_merged.fastq/,"",$0); print}' | sort | uniq)

for sample in $samples
do
	diamond blastx -e 1e-10 --id 40 -k 1 --threads 10 -d $RPOB_DB -q ${sample}_clean_merged.fastq -o ${sample}_rpob.csv --outfmt 6
done
```
</details>

Next, put rp_abun.py in the same folder as your DIAMOND outputs.
<details>
  <summary>rp_abun.py</summary>
  
``` 
import os
import argparse
import pandas as pd

import warnings
warnings.simplefilter(action='ignore', category=Warning)

class annotated(object):
    def __init__(self, param):
        self.file = param['filename']
        self.length_file = param['length_file']
        self.identity = float(param['identity'])
        self.mlen = int(param['mlen'])
        self.evalue = float(param['evalue'])
 

    def data_process(self):
        data = pd.read_csv(self.file, sep = '\t', header = None) 
        data.sort_values([0,11], inplace = True, ascending = False)
        data.drop_duplicates(subset = 0, keep = 'first', inplace = True)
        data.columns = ['id', 'sub_id', 'identity', 'alignLen', 'mismat', 'gapOpens', 'qStart', 'qEnd', 'sStart', 
                'sEnd', 'evalue', 'bit']

        data = data[data.identity >= self.identity]
        data = data[data.alignLen >= self.mlen]
        data = data[data.evalue <= self.evalue]
        return data
    
    def length_process(self, database):
        len_data = pd.read_csv(self.length_file, sep = '\t', header = None) 
        len_data.columns = ['sub_id', 'gene_len']
        len_temp = len_data['sub_id'].str.split('|', expand=True)
        len_temp.columns = database.colum_names
        len_temp = len_temp.drop(columns = database.dropped_column)
        len_temp['gene_len'] = len_data.iloc[:, 1]
        len_temp['gene_len(bp)'] = (len_temp['gene_len']+1)*3
        return len_temp
    
    def combine_length(self, data, len_data, database):
        temp = data['sub_id'].str.split('|', expand=True)
        temp.columns = database.colum_names
        temp = temp.drop(columns = database.dropped_column)
        temp['length'] = data.iloc[:, 3]
        data = temp.groupby(database.keep_column)['length'].agg([('count','count'), ('length','sum')]).reset_index()
        data = pd.merge(data, len_data, how = "left", on = database.keep_column)
        return data
        

class Database(object):
    def __init__(self, name):
        self.colum_names = []
        self.dropped_column = []
        self.keep_column = []
        self.set_properties(name)
        
    def set_properties(self, name):
        if name == 'deeparg':
            self.colum_names = ['protein_accession', 'extra', 'gene', 'drug', 'gene_family']
            self.dropped_column = ['extra']
            self.keep_column = ['gene', 'drug', 'protein_accession', 'gene_family']
        elif name == 'card':
            self.colum_names = ['extra','protein_accession', 'aro_index', 'gene']
            self.dropped_column = ['extra']
            self.keep_column = ['gene', 'protein_accession', 'aro_index']
        elif name == 'rpob':
            self.colum_names = ['extra','accession_num', 'rpob_gene']
            self.dropped_column = ['extra']
            self.keep_column = ['accession_num', 'rpob_gene']
        elif name == 'mobileOG':
            self.colum_names = ["mobileOG Entry Name", "col1", "col2", "mge_class", "mge_subclass", "col5", "col6"]
            self.dropped_column = ["col1", "col2", "col5", "col6"]
            self.keep_column = ['mobileOG Entry Name', 'mge_class', 'mge_subclass']


def cal_abundance(arg_param, rpob_param, database, output_file = ""):
    arg_db = Database(database)
    rpob_db = Database('rpob')
    arg_obj = annotated(arg_param)
    rpob_obj = annotated(rpob_param)
    
    arg_data = arg_obj.data_process()
    arg_len_data = arg_obj.length_process(arg_db)
    arg_data = arg_obj.combine_length(arg_data, arg_len_data, arg_db)
    rpob_data = rpob_obj.data_process()
    rpob_len_data = rpob_obj.length_process(rpob_db)
    rpob_data = arg_obj.combine_length(rpob_data, rpob_len_data, rpob_db)
    arg_db.keep_column.append("count")
    output = arg_data[arg_db.keep_column]

    Nrpob = rpob_data.loc[rpob_data['length'] >= 0]['count'].sum()
    Lrpob = rpob_len_data['gene_len(bp)'].mean()
    

    output['rpob_Normalization'] = (arg_data['count']/arg_data['gene_len(bp)'])/(Nrpob/Lrpob)
    output['rpob_count'] = Nrpob 
    output['sample'] = os.path.splitext(os.path.basename(arg_param['filename']))[0]
    if not output_file:
        output_file = os.path.join(os.path.dirname(arg_param['filename']), os.path.splitext(os.path.basename(arg_param['filename']))[0] + "_abundance.txt")
    
    
    output.to_csv(output_file, sep = "\t", index = False)

    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--arg", type = str, required = True, help = "Path to arg annotation file (required)")
    parser.add_argument("-la", "--len_arg", type = str, required = True, help = "Path to arg database length file (required)") 
    parser.add_argument("-r", "--rpob", type = str, required = True, help = "Path to rpob annotation file (required)")
    parser.add_argument("-lr", "--len_rpob", type = str, help = "Path to rpob database length file")
    parser.add_argument("-o", "--out", type = str, help = "Path to output file")

    parser.add_argument("--db", type = str, default = "deeparg", help = "name of ARG database such as card/deeparg [default deeparg]")
    parser.add_argument("--arg_identity", type = float, default = 80, help = "minimum identity for alignments [default 80]",)
    parser.add_argument("--arg_mlen", type = float, default = 25, help = "diamond minimum length for considering a hit [default 25aa]",)
    parser.add_argument("--arg_evalue", type = float, default = 1e-10, help = "minimum e-value for alignments [default 1e-10]",)

    parser.add_argument("--rpob_identity", type = float, default = 40, help = "minimum identity for alignments [default 40]",)
    parser.add_argument("--rpob_mlen", type = float, default = 25, help = "diamond minimum length for considering a hit [default 25aa]",)
    parser.add_argument("--rpob_evalue", type = float, default = 1e-10, help = "minimum e-value for alignments [default 1e-10]",)

    args = parser.parse_args()
    db_name = args.db
    arg_parameters = dict(
        filename = args.arg,
        length_file = args.len_arg, 
        identity = args.arg_identity,
        mlen = args.arg_mlen,
        evalue = args.arg_evalue
    )

    rpob_parameters = dict(
        filename = args.rpob,
        length_file = args.len_rpob,
        identity = args.rpob_identity,
        mlen = args.rpob_mlen,
        evalue = args.rpob_evalue
    )
    
    
    cal_abundance(arg_param = arg_parameters, rpob_param = rpob_parameters, database = db_name, output_file = args.out)

if __name__ == '__main__':
    main()
```
</details>


Then, run cal_abundances.sh. You'll first need to create an environment named "mypy3" and install python and pandas. To do so, use the following lines of code.
```
conda create -n mypy3 python=3.12 pip 
source activate mypy3
conda install ipykernel
pip install plotly kaleido
conda install pandas
```
Run cal_abundances.sh to calculate the abundance of each ARG and RPOB gene in every sample.
<details>
<summary> cal_abundances.sh</summary>

```
#!/usr/bin/sh

#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --mem=20G
#SBATCH -t 1:00:00
#SBATCH -n 50
#SBATCH -N 1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yourPID@vt.edu

module load Miniconda3

eval "$(conda shell.bash hook)"

conda activate mypy3

#REF_arg='/projects/ciwars/ChezLiz/CARD4.0.1.dmnd'
REF_arg='/projects/ciwars/ChezLiz/protein_fasta_protein_homolog_model.len.txt' ##comment either this out or the deeparg length file depending on what you used. 
#REFA='/projects/ciwars/databases/bacmet_len.txt'

#REF_16s='/projects/ciwars/databases/gg_13_5.len'
REF_rpob='/projects/ciwars/databases/RpoB.ref.len'

#samples=$(ls $src/projects/ciwars/avdarling/mergedreadsfinal/fastp/S2* | awk '{split($_,x,"_fastp.fastq.gz"); print x[1]}') #change path to where your files are located and awk line depending on how you named your files
#samples=$(ls /projects/ciwars/thomasbyrne/*_clean_merged.fastq | awk '{split($_,x,"_clean_merged.fastq"); print x[1]}') #for my subset
cd /projects/intro2gds/I2GDS2025/G6_AMR_ARG/Daniel/Subset/diamond_output
samples=$(ls S*_arg_full.csv | awk -F/ '{gsub(/_arg_full.csv/, "", $NF); print $NF}' | sort | uniq)

for sample in $samples; do
    #sample=$(basename "$sample")
    #echo $sample
    printf "%s \n" ${sample}
    file=$(basename -- ${sample})
    printf "%s \n" ${file}
    python rp_abun.py -a ${sample}_arg_full.csv -r ${sample}_rpob.csv -la $REF_arg -lr $REF_rpob --rpob_identity 40 --db card -o /projects/intro2gds/I2GDS2025/G6_AMR_ARG/Daniel/Subset/diamond_output/${sample}_abundance_new.csv
done
```
</details>

Finally, run merging_my_noramlized.sh. This script merges the previous abundance_new.csv files for each sample into a single, final output: I2GDS_G6_AMR_Diamond.txt. This output will be one of the inputs for the R script.
<details>
<summary> merging_my_noramlized.sh</summary>

```
#!/bin/bash

#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --mem=1G
#SBATCH -t 1:00:00
#SBATCH -n 50
#SBATCH -N 1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yourPID@vt.edu
# Output file name

output_file="I2GDS_G6_AMR_Diamond.txt"

cd /projects/intro2gds/I2GDS2025/G6_AMR_ARG/Daniel/Subset/diamond_output

# Find and merge files

for file in *_abundance_new.csv; do
    # Skip the output file to avoid self-merging
    if [ "$file" != "$output_file" ]; then

        # Print filename as a comment in the merged file. nope i made it a column so that i can just treat it as a csv directly
        #echo "# $file" >> "$output_file"

        # Append contents to the merged file, skipping the header line
        #tail -n +2 "$file" >> "$output_file"

        #so i think i just make it 1 instead of 2 but we'll see ig
        tail -n +1 "$file" >> "$output_file"
    fi
done


echo "Merge complete. Merged file: $output_file"
```

</details>

### 7.2 Working in R
Inputs: `CARD4.0.1_aro_cat.csv`, `I2GDSMetadata.csv`, `Antibiotic to Drug.csv`, and `I2GDS_G6_finalDiamond.csv`
All inputs are provided in the code-files folder. I2GDS_G6_finalDiamond.csv is the exact same as the txt file output from the previous merging step, just converted to a csv file using Excel and renamed.

Libraries needed: `tidyverse`, `RColorBrewer`
To install these libraries, open RStudio and click the "Tools" drop-down menu at the top of the window. Then select "Install Packages".
In the "Packages" box, type tidyverse and click install. Do the same for RColorBrewer.

In R, you will take the output from the previous merging step (converted from txt to csv file) and use additional metadata, ARG annotation, and Antibiotic to Drug categorization files to created a stacked bar chart. This stacked bar chart shows which ARGs were found in different sample types, classifying the ARGs by drug resistance.

To produce this stacked bar graph, run the R code below. Comments within the code explain the general workflow.

<details>
<summary> I2GDS_G6.R</summary>
	
```
#Loading packages
library(tidyverse)
library(RColorBrewer)

#Set working directory - change yourself!
setwd("C:/Users/Jvlan/Downloads")

#Reading in all data needed
#The drugs object comes from the deeparg database, where they assigned each specific antibiotic to a broader class
#The CARD_aro object is from the CARD database, and shows which gene confers resistance to each specific antibiotic
CARD_aro <- read.csv("CARD4.0.1_aro_cat.csv")
Args <- read.csv("I2GDS_G6_finalDiamond.csv")
metadata <- read.csv("I2GDSMetadata.csv")
drugs <- read.csv("Antibiotic to Drug.csv")

#To make the graph look better, we assign the genes in the CARD database to a broader drug class using the deeparg conversion
for (i in 1:nrow(CARD_aro)) {
  for (j in 1:nrow(drugs)) {
    if(grepl(drugs[j,1], CARD_aro[i, 4])) {CARD_aro[i, "Drug"] <- drugs[j,2]}
  }
}

for (i in 1:nrow(CARD_aro)) {
  if(grepl(";", CARD_aro[i, 4])) {CARD_aro[i, "Drug"] <- "Multidrug"}
}

CARD_aro$Drug[is.na(CARD_aro$Drug)] <- "Other"

#Merging the DIAMOND data and the CARD metadata about each gene
Args <- merge(Args, CARD_aro, by.x = "protein_accession", by.y = "Protein.Accession")

#Renaming the samples to only include the same name
Args$sample = substr(Args$sample, 1, nchar(Args$sample)-9)

#Merging the DIAMOND data and sample metadata
Args <- merge(Args, metadata, by.x = "sample", by.y = "Run")

Args$rpob_Normalization[Args$rpob_Normalization==""] <- 0

#Assigning data its appropriate type
Args$Drug.Class <- as.factor(Args$Drug.Class)
Args$isolation_source <- as.factor(Args$isolation_source)
Args$count <- as.numeric(Args$count)
Args$rpob_Normalization <- as.numeric(Args$rpob_Normalization)
Args$Drug <- as.factor(Args$Drug)

#Grouping the DIAMOND data to find how often each type of resistance was in each sample
df <- Args %>%
  group_by(sample, isolation_source, Drug) %>%
  summarize(sum_count = sum(count))

#Checking most abundant drug class
a <- df %>%
  group_by(Drug) %>%
  summarize(drug_count = sum(sum_count))

#Reassigning the orders to make the graph look better. Drug was done in order of highest to lowest counts, isolation source to match the paper
df$Drug <- factor(df$Drug, levels = c("Multidrug", "Other", "tetracycline", "macrolide-lincosamide-streptogramin", "aminoglycoside", "aminocoumarin",
                                      "peptide", "sulfonamide", "quinolone", "beta_lactam",
                                      "phenicol", "mupirocin", "pleuromutilin"))

df$isolation_source <- factor(df$isolation_source, levels = c("Influent", "Primary", "Act. Sludge",
                                                              "Secondary", "Floc Sed", "Ozonation",
                                                              "BAC/GAC", "Chlorination", "3 days", "Background"))


#Creating an object so we don't have to use R's ugly default colors
n <- 13
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#Actually creating the stacked bar chart
ggplot(df, aes(x = isolation_source, y = sum_count, fill = Drug)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values=col_vector) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Total ARG copies") +
  xlab("")
```

</details>

The output should look like:
<img width="1891" height="1087" alt="image" src="https://github.com/user-attachments/assets/a9bf7d2b-53ac-4e44-a142-add51023255c" />

The paper's figure looks like:
<img width="720" height="548" alt="image" src="https://github.com/user-attachments/assets/86f76cf9-ff1b-4c56-8485-d64415595501" />

We did not normalize to 16S gene copies, as that data was not easily accessable. We also only used 10 samples from 1 sampling event, while the paper had multiple. These differences could explain why the data looks slightly different.
