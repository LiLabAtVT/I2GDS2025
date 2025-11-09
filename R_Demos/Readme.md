# 🧬 Single-Nucleus RNA-Seq Analysis of *Arabidopsis thaliana* Roots Infected with *Phytophthora capsici* #

This repository contains teaching materials and analysis scripts for processing and analyzing single-nucleus RNA-seq (snRNA-seq) data from *Arabidopsis thaliana* roots at 24 hours post-infection with the oomycete pathogen Phytophthora capsici.
The workflow demonstrates how to perform sample integration, clustering, marker identification, differential expression, and visualization using Seurat.

## 🌱 Project Overview ##
Component	Description
Goal	Understand how Arabidopsis root cell types respond to pathogen infection.
Approach	snRNA-seq profiling and Seurat-based integration pipelines.
Outcome	Identification of 12 transcriptionally distinct cell clusters and cell-type-specific immune responses.

This repository is designed as a teaching example, focusing on the essential analysis and visualization steps.

# 🧩 Requirements #

Before running the workflow, install the following R packages:
```bash
# CRAN packages
install.packages(c("ggplot2", "dplyr", "cowplot"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Seurat", "clusterProfiler", "org.At.tair.db", "enrichplot"))

📁 Repository Structure
├── Integration.R                      # Data integration, PCA, clustering, UMAP
├── UMAP_Cluster_ID_Script.R           # Assign biological cluster identities

├── CellMarkers_DEG_Analysis.R         # Cluster markers + Differential expression

├── DotPlot_FunctionalGenes.R          # Visualization of defense/hormone pathways
├── Feature_plot_Celltype_Markers.R    # UMAP plots of selected cell type markers
```

## 🧬 Pipeline Steps (Teaching Version) ##

# Step 1 — Data Integration & UMAP Clustering #

Script: Integration.R

Performs SCTransform-based sample integration.

Runs PCA, clustering, and UMAP dimensional reduction.

Outputs:

Integrated_Seurat_Object.rds

UMAP_Clusters.pdf (or .png depending on export)

# Step 2 — Cluster Marker Genes & Differential Expression #

Script: CellMarkers_DEG_Analysis.R

Identifies cluster-specific marker genes.

Performs Infected vs. Non-infected differential expression analysis.

Outputs:

All_Markers.csv

DEG_Infected_vs_Control.csv

# Step 3 — Visualization of Biological Responses #

Scripts:

Script	Visualization Generated
DotPlot_FunctionalGenes.R	Hormone / defense signaling gene expression patterns
Feature_plot_Celltype_Markers.R	UMAP plots showing key marker genes by cluster

🧭 Summary of Biological Insight


| Analysis Step                | Biological Interpretation                                               |
| ---------------------------- | ----------------------------------------------------------------------- |
| **Integration & Clustering** | Defines transcriptionally distinct root cell populations.               |
| **Marker Analysis**          | Identifies genes defining each cell identity and infection response.    |
| **Visualization**            | Reveals where plant defense pathways are activated at the tissue level. |
