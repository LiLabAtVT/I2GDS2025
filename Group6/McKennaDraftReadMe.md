## vsearch
### Merging your trimmed and decontaminated reads

**Overview**\
Use your bbduk output files (_Sample Name_\_[1 or 2]\_decontam.fastq.gz) and merge the forward (\_1\_) and reverse (\_2\_) reads for each sample. Your input should be the bbduk output files you generated from the previous step, or the contingency bbduk output files provided.

**Step 1**\
Create an environment where you'll install vsearch, activate said environment
```
module load Miniconda3/24.7.1-0
conda create vsearchenv
source activate vsearchenv
```

You should now see (vsearchenv) \[yourPID@clustername directoryname\]$, indicating that you've entered the vsearchenv environment

**Step 2**\
Install vsearch in environment
```
conda install bioconda::vsearch
```
Press "y" when prompted to complete installation

**Step 3**\
Use code below to access bbduk output files and merge forward and reverse reads. Be sure to replace the PID and directory path for your input files.
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

<ins>Troubleshooting</ins>: If you receive an "Out of Memory" error after submitting the slurm job, try decreasing threads to 10 and adding\
#SBATCH --cpus-per-task=16\
<ins>Code Explanation</ins>: This for loop goes through each sample's forward and reverse read and merges them. Due to unequal read lengths, some of the forward and reverse reads will not be merged, so in addition to the "merged" output, there is also an "unmerged forward" and "unmerged reverse" output. The merged, unmerged forward, and unmerged reverse reads for a sample are concatenated to output a "clean merged" file. There is no longer a use for the sample's merged, unmerged forward, and unmerged reverse reads, so they are removed from the directory.

## DIAMOND
### Using DIAMOND to annotate your reads and identify ARGs

**Overview**\
You will be running each sample read against a database of known ARG sequences to identify the ARGs present in each sample. Your input should be the vsearch output files you generated from the previous step, or the contingency vsearch output files provided.

**Step 1**\
Create an environment where you'll install DIAMOND, activate said environment
```
module load Miniconda3/24.7.1-0
conda create diamondenv
source activate diamondenv
```
You should now see (diamondenv) \[yourPID@clustername directoryname\]$, indicating that you've entered the diamondenv environment

**Step 2**\
Install DIAMOND in environment
```
conda install bioconda::diamond
```
Press "y" when prompted to complete installation

**Step 3**\
Download your reference database for ARG identification

For this data, we utilitze the Comprehensive Antibiotic Resistance Database (CARD), which is a bioinformatic database of resistance genes, their products, and associated phenotypes.

Go to this link [https://card.mcmaster.ca/](url), click Download, and scroll to the section that says "Download CARD Data" (NOT "Download CARD-R Resistomes, etc.)
The reference paper that collected this data annotated it using CARD v3.0.8, so we will do the same. In the "Download CARD Data" section, click "More data downloads..." and scroll down until you find version 3.0.8. Download this compressed folder and extract all files. Then, upload (mention protein homolog file)
