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

You should now see (vsearchenv) [yourPID@ ____ to indicate that you've entered the environment

**Step 2**\
Install vsearch in environment
```
conda install bioconda::vsearch
```
Press "y" when prompted to complete installation
