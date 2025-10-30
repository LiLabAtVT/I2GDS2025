## Flowchart
<img width="593" height="778" alt="image" src="https://github.com/user-attachments/assets/6f9dc70f-19f5-4431-bddc-5a1dc2c58afc" />



## 3. Trimmomatic (Cleaning Data)
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
You can customize how stringent you want your trimming of adapters to be. The code uses the default parameters \
- ILLUMINACLIP:fastaWithAdaptersEtc:seed mismatches:palindrome clip threshold:simple clip threshold
	- fastaWithAdaptersEtc: specifies the path to a fasta file containing all the adapters. This is TruSeq3-PE.fa.
	- seedMismatches: specifies the maximum mismatch count which will still allow a full match to be performed
	- palindromeClipThreshold: specifies how accurate the match between the two 'adapter ligated' reads must be for PE 			palindrome read alignment.
	- simpleClipThreshold: specifies how accurate the match between any adapter etc. sequence must be against a read.

```
ILLUMINACLIP:${ADAPTERS}:2:30:10
```

Description from [Trimmomatic Website](http://www.usadellab.org/cms/?page=trimmomatic)

**Trimmomatic parameters** \
You can customize how stringent you want your cutting to be. The code uses the default parameters \
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
With the 10 input files, this should take ~1 hour to run
<details>
<summary> Trimmomatic.sh</summary>

```
#!/bin/bash
#SBATCH --account=YourAccount #Change
#SBATCH --partition=normal_q
#SBATCH --mem=16G
#SBATCH -t 30:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=START,END,FAIL
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

## 4. bbduk (Decontaminating data)

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

### 4.5.2 bbduk script
<details>
<summary> bbduk.sh</summary>

```
#!/bin/bash

#SBATCH --account=YourAccount #Change
#SBATCH --partition=normal_q
#SBATCH --mem=200G
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=START,END,FAIL
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
For each input file, there will be one output FASTQ file, and for each pair of FASTQ files there will be one txt file. This file shows how many of each sequence in the reference FASTA file are removed (should be ~0.5-1%)
