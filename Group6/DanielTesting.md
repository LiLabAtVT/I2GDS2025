## Flowchart
<img width="320" height="831" alt="image" src="https://github.com/user-attachments/assets/e529bdb8-67f6-4f32-94da-f3e258bca886" />

## Trimmomatic
<details>
<summary> script</summary>


```
#!/bin/bash

#SBATCH --account=prudenlab
#SBATCH --partition=normal_q
#SBATCH --mem=16G
#SBATCH -t 30:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=dceglio@vt.edu

module load Miniconda3

source activate trimmomatic

cd /projects/intro2gds/I2GDS2025/G6_AMR_ARG/SampleSubset

OUTPUT_DIR=/projects/intro2gds/I2GDS2025/G6_AMR_ARG/Daniel/Subset/trimmomatic_output
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


## bbduk

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
