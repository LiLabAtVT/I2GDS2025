## Introduction
This project aims to replicate the pipline for analyzing a set of bacterial genomes (Xylella fastidiosa (Xf)) from pure cultures and herbarium specimen

Reference: Century-old herbarium specimen provides insights into Pierce’s disease of grapevines emergence in the Americas
https://doi.org/10.1016/j.cub.2024.11.029

## Overview of workflow

00. data download
    01. Metagenomes of herbarium specimen
    02. Raw reads of 44 modern strains of Xf
    03. Reference genome of Xf

1. Bowtie2 map Xf genome in herbarium specimen sequence
2. SPAdes
   2.1 Assemble mapped reads from herbarium specimen
   2.2 Assemble raw reads of modern strains       
3. CheckM for quality control of assembled sequences
4. Prokka for gene annotation


## Download raw sequence and reference genome





## Mapping Xf genome in herbarium via Bowtie2



## SPAdes


## CheckM


## Prokka
