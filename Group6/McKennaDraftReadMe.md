## vsearch
### Merging your trimmed and decontaminated reads

**Overview**\
Use your bbduk output files (_Sample Name_\_[1 or 2]\_decontam.fastq.gz) and merge the forward (\_1\_) and reverse (\_2\_) reads for each sample. Your input should be the bbduk output files you generated from the previous step, or the contingency bbduk output files provided.

**Step 1**\
Create an environment to download vsearch, activate said environment
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

<ins>Code Explanation</ins>: This for loop goes through each sample's forward and reverse read and merges them. Due to unequal read lengths, some of the forward and reverse reads will not be merged, so in addition to the "merged" output, there is also an "unmerged forward" and "unmerged reverse" output. The merged, unmerged forward, and unmerged reverse reads for a sample are concatenated to output a "clean merged" file. There is no longer a use for the sample's merged, unmerged forward, and unmerged reverse reads, so they are removed from the directory.
