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

## Data download





## Mapping Xf genome in herbarium via Bowtie2



## SPAdes


## CheckM


## Prokka
