## Introduction
This project aims to replicate the pipline for analyzing a set of bacterial genomes (Xylella fastidiosa (Xf)) from pure cultures and herbarium specimen

Reference: Century-old herbarium specimen provides insights into Pierce’s disease of grapevines emergence in the Americas
https://doi.org/10.1016/j.cub.2024.11.029

## Overview of workflow

1. **Data download**
   - Metagenomic reads from herbarium specimen(s)
   - Raw reads of 44 modern Xf strains
   - Reference genome(s) of *X. fastidiosa*

2. **Mapping Xf in metagenomes**
   - Use **Bowtie2** to map herbarium metagenomic reads to the Xf reference genome

3. **Assembly**
   - **SPAdes**
     - 2.1 Assemble mapped reads from herbarium specimen
     - 2.2 Assemble raw reads of modern strains

4. **Quality control**
   - **CheckM** for QC of assembled genomes/scaffolds 

5. **Annotation**
   - **Prokka** for gene annotation

## Environment setup

## Data download - SRAtools
1.1 Metagenome dataset (PRJNA1114123) and reference sequence genome (GCA_000007245.1) was directly downloaded from NCBI


1.2 The accession list of 44 modern strains was download from NCBI and then using SRAtools for downloading. 


1.3 As downloaded files were in ".sra" format, command "fasterq-dump" was used for transfering SRA to FASTQ 



## Mapping Xf genome in herbarium sample - Bowtie2
As Xf is a bacterial pathogen colonizing in xylem vessel, it is possible to sequece and assemble its historical genome from historical material. Reads were mapped against the reference geonme sequence using Bowtie2.




## SPAdes


## CheckM


## Prokka
