## Flowchart
<img width="568" height="725" alt="image" src="https://github.com/user-attachments/assets/8700e9e5-7c02-42b6-a23a-876d598eb109" />


## Trimmomatic (Cleaning Data)
### X.1 Introduction
Trimmomatic is a tool used to clean data of low quality reads, and is designed specifically for Illuminia paired-end or single-end sequencing data.\
For more info, visit the [Trimmomatic Website](http://www.usadellab.org/cms/?page=trimmomatic)

### X.2 Creating Environment
First, create a new environment on ARC. Naming it "trimmomatic" matches the script below:
```
conda create -n trimmomatic
```

After creating the environment, activate it:
```
source activate trimmomatic
```

### X.3 Downloading Trimmomatic
Next, download trimmomatic:
```
conda install bioconda::trimmomatic
```

### X.4 Inputs Required
- **Input data**: Paired end FASTQ/FASTQ.gz files
- **Illuminia Specific Adapter File**: Adapter file (TruSeq3-PE.fa)

### X.5 Trimmomatic Script

#### X.5.1 Explanation of Trimmomatic.sh
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

#### X.5.2 Script
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

### X.6 Output Explanation
Each foward and reverse file will give 2 output files:
- paired: These reads survived trimming and their "partner" read survived trimming
- unpaired: These reads survived trimming, but their "partner" read did not

Typically, analyze only moves forward with the paired reads. However, unpaired reads can still be used for some future analysis.

## bbduk (Decontaminating data)

<details>
<summary> script</summary>

```
#!/bin/bash

#SBATCH --account=prudenlab
#SBATCH --partition=normal_q
#SBATCH --mem=200G
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=dceglio@vt.edu

module load Miniconda3

source activate bbduk

cd /projects/intro2gds/I2GDS2025/G6_AMR_ARG/Daniel/Subset/trimmomatic_output

sample=$(ls *_[12]_paired.fastq.gz | awk -F/ '{gsub(/_[12]_paired.fastq.gz/, "", $NF); print $NF}' | sort | uniq)

for one_sample in $sample; do
    echo ${one_sample}

    bbduk.sh preallocate=t -Xmx180G \
        minlength=51 \
        in="${one_sample}"_1_paired.fastq.gz \
        in2="${one_sample}"_2_paired.fastq.gz \
        out=/projects/intro2gds/I2GDS2025/G6_AMR_ARG/Daniel/Subset/bbduk_output/"${one_sample}"_1_decontam.fastq.gz \
        out2=/projects/intro2gds/I2GDS2025/G6_AMR_ARG/Daniel/Subset/bbduk_output/"${one_sample}"_2_decontam.fastq.gz \
        stats=/projects/intro2gds/I2GDS2025/G6_AMR_ARG/Daniel/Subset/bbduk_output/"${one_sample}"_decontamination_stats.txt \
        ref='/projects/apruden_lab/thomasbyrne/metagenomics_workshop/merged_ref_5081444234059102403.fa.gz'


done

```



  
</details>
