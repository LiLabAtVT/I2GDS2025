install.packages("ape")
install.packages("phangorn")
install.packages("ggtree")
install.packages("BiocManager")
install.packages("dplyr")
install.packages("ggplot2")
BiocManager::install("ggtree")

library(ape)
library(phangorn)
library(ggtree)
library(dplyr)
library(ggplot2)

#set working directory
setwd("~/Desktop/Scripts/Rscript")

# 1. read core gene alignment file
alignment <- read.FASTA("core_gene_alignment_filtered.aln")

# remove "_prokka" in rownames
names(alignment) <- gsub("_prokka$", "", names(alignment))

# 2. construct NJ tree
dist_matrix <- dist.dna(alignment, model = "raw")
nj_tree <- nj(dist_matrix)

# 3. convert to phyDat format (for Maximum likelihood tree)
dna <- phyDat(alignment, type = "DNA")
fit <- pml(nj_tree, dna)

# 4. optimize model
fit_GTR <- update(fit, k = 4, inv = 0.2)
fit_opt <- optim.pml(fit_GTR, optNni = TRUE, optBf = TRUE, optQ = TRUE,
                     optGamma = TRUE, optInv = TRUE)


# 5. ML tree
ml_tree <- fit_opt$tree

# 6. bootstrap
bs_trees <- bootstrap.pml(fit_opt, bs = 100, optNni = TRUE, multicore = FALSE)
bs_values <- prop.clades(ml_tree, bs_trees) * 100
ml_tree$node.label <- round(bs_values, 0)


# 7. read metadata
metadata <- read.csv("SraRunTable.csv", stringsAsFactors = FALSE)
metadata <- metadata[,c(1,21,35)]
colnames(metadata) <- c("tip", "isolate", "Subspecies")

metadata$Subspecies[metadata$Subspecies == "Xylella fastidiosa subsp. multiplex"] <- "multiplex"
metadata$Subspecies[metadata$Subspecies == "Xylella fastidiosa subsp. fastidiosa"] <- "fastidiosa"


# 8. plot tree（ggtree）
tree_plot <- ggtree(ml_tree,branch.length = "scaled") %<+% metadata +
  geom_tiplab(aes(label = isolate, color = Subspecies), size = 3) +
  scale_x_continuous(limits = c(0, 0.028)) +
  theme_tree2() +
  scale_color_manual(values = c("multiplex" = "blue",
                                "fastidiosa" = "darkred")) +
  theme(legend.position = c(0.85, 0.2)) +
  labs(color = "Subspecies")


tree_plot

