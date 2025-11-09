# 🧬 Single-Nucleus RNA-Seq Analysis of Arabidopsis thaliana Roots Infected with *Phytophthora capsici*#

This repository contains teaching materials and analysis scripts for processing and analyzing single-nucleus RNA-seq (snRNA-seq) data from Arabidopsis thaliana roots at 24 hours post-infection with the oomycete pathogen Phytophthora capsici.
The workflow demonstrates how to perform sample integration, clustering, marker identification, differential expression, and visualization using Seurat.

## 🌱 Project Overview##

Goal: Explore how Arabidopsis root cell types respond to pathogen infection.

Approach: snRNA-seq profiling and Seurat-based data integration.

Outcome: Identification of 12 transcriptionally distinct cell clusters and cell-type-specific immune responses.

This repository is structured as a teaching example, highlighting key steps needed to analyze snRNA-seq data and visualize biological insights.

🧩 Requirements

Before running the workflow, ensure R and the required packages are installed.

# CRAN packages
install.packages(c("ggplot2", "dplyr", "cowplot"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Seurat", "clusterProfiler", "org.At.tair.db", "enrichplot"))

📁 Repository Structure
├── 1_Integration_UMAP/
│   ├── Integration.R         # Integration, PCA, UMAP, clustering
│   ├── UMAP_ClusterID.R                 # Assign biological cluster identities
│
├── 2_CellMarkers_DEG/
│   ├── CellMarkers_DEG_Analysis.R       # Marker genes and differential expression
│
├── 3_Visualization/
│   ├── DotPlot_FunctionalGenes.R        # Hormone/defense gene dot plots
│   ├── FeaturePlot_CellMarkers.R        # Marker gene UMAP visualization
│   └── FeaturePlot_SelectedMarkers.R    # Additional gene-level feature plots

## 🧬 Pipeline Steps (Teaching Version)##
# Step 1 — Data Integration & UMAP Clustering #

Script: 1_Integration_UMAP/Integration_Seurat_SCT.R

Integrates infected and control datasets using SCTransform.

Performs PCA, chooses dimensions, and clusters nuclei.

Generates UMAP visualization of root cell populations.

#Outputs:#

Integrated_Seurat_Object.rds

UMAP_Clusters.png

# Step 2 — Marker Gene Identification & DE Analysis #

Script: 2_CellMarkers_DEG/CellMarkers_DEG_Analysis.R

Finds cluster-specific marker genes.

Detects differentially expressed genes (Infected vs. Control).

Outputs:

All_markers.csv

All_markers.csv

# Step 3 — Visualization of Biological Responses #

Scripts:
3_Visualization/DotPlot_FunctionalGenes.R
3_Visualization/FeaturePlot_CellMarkers.R

Dot plots for top marker genes.

UMAPs of selected marker genes and cell identity signatures

Outputs:

Script	Visualization
DotPlot_FunctionalGenes.R	Defense / hormone expression patterns
FeaturePlot_CellMarkers.R	Key marker expression on UMAP

# 🧭 Summary of Biological Insight #
Analysis Step	What We Learn
Integration & Clustering	Defines root cell populations and structure
Marker Analysis	Determines cell-type identity and infection response
Visualization	Reveals where defense pathways are activated in root tissue
