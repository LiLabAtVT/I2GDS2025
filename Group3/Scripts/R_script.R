library(data.table)
library(Matrix)
library(irlba)
library(FNN)
library(umap)
library(ggplot2)

# read files
bin_dir <- "/projects/lu_lab/Rui/qiang_NC_data/11.05/09_allc/bin1Mb/" 

files <- list.files(bin_dir, pattern = "*.tsv", full.names = TRUE)
cell_ids <- gsub(".bin1Mb.tsv", "", basename(files))

message("Reading bin files ... ")

bin_list <- list()

for (i in seq_along(files)) {
  dt <- fread(files[i])
  setnames(dt, c("bin_size","bin_id","mC","C","rate"))
  bin_list[[cell_ids[i]]] <- dt[, .(bin_id, rate)]
}
# feature matrix（cells × bins）
all_bins <- unique(unlist(lapply(bin_list, function(x) x$bin_id)))
mat <- Matrix(0,
              nrow = length(bin_list),
              ncol = length(all_bins),
              sparse = TRUE)
rownames(mat) <- names(bin_list)
colnames(mat) <- all_bins
message("Building feature matrix ...")
for (i in seq_along(bin_list)) {
  df <- bin_list[[i]]
  idx <- match(df$bin_id, all_bins)
  mat[i, idx] <- df$rate
}
# Global methylation normalization
message("Normalizing ...")
global_rate <- rowMeans(mat, na.rm = TRUE)
mat_norm <- mat / global_rate
# PCA
message("Running PCA ...")
set.seed(123)
pca_res <- irlba(mat_norm, nv = 50)
pc_scores <- pca_res$u %*% diag(pca_res$d)
rownames(pc_scores) <- rownames(mat_norm)

# UMAP 

message("Running UMAP ...")
um_res <- umap(pc_scores, n_neighbors = 30, min_dist = 0.3)

umap_df <- data.frame(
  UMAP1 = um_res$layout[,1],
  UMAP2 = um_res$layout[,2],
  cell = rownames(pc_scores)
)
# K-means 

message("Running k-means clustering (k=3) ...")
set.seed(123)

km3 <- kmeans(umap_df[, c("UMAP1","UMAP2")], centers = 3, nstart = 50)
umap_df$cluster3 <- factor(km3$cluster)
# plot

p <- ggplot(umap_df, aes(UMAP1, UMAP2, color = cluster3)) +
      geom_point(size = 1.5, alpha = 0.9) +
      theme_classic() +
      ggtitle("Drop-BS Single Cell Clustering (K-means = 3)") +
      scale_color_brewer(palette = "Set1")

print(p)

ggsave("dropbs_kmeans3_umap.png", p, width = 7, height = 6, dpi = 300)

message("Done! Output saved: dropbs_kmeans3_umap.png")
