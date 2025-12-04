#Bacteria Figure Script
install.packages(phyloseq)
install.packages(ggplot2)

library(phyloseq)
library(ggplot2)


setwd("~/woRk/") #change to your directory name
bacteria_phylum <- read.csv("~/woRk/apr24_pp/level2-table.tsv", sep = "\t") #change filepath

#change "phylum" to "bacteria_phylum"
phylum$Taxonomy <- sub("^.*p_", "", phylum$X.OTU.ID)
phylum$Taxonomy[startsWith(phylum$Taxonomy, "k")] <- "Unassigned"
phylum$Taxonomy <- gsub("_","",phylum$Taxonomy)

#making it relative
#change "phylum" to "bacteria_phylum"
for (i in 2:66) {phylum[,i] <- phylum[,i]/sum(phylum[,i])
}

#pivot to long data frame
#change "phylum" to "bacteria_phylum"
library(tidyr)
long_phylum <- pivot_longer(phylum,cols = starts_with("Sample"), names_to = "samples", values_to = "abundance")

#plot with ggplot2: stacked & percent
library(ggplot2)
ggplot(data = long_phylum, aes(x = samples, y = abundance, fill = Taxonomy)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() + 
  theme(legend.position="bottom", legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size  +
guides(fill=guide_legend(nrow=5, byrow=TRUE))

# save plot: run the line with your name :)
ggsave(file="G1-bacteria-phylum.png") #for Peter


#Bacteria Figure Script
install.packages("tidyverse")
install.packages("ggplot2")

library(tidyverse)
library(ggplot2)


# folder you will be working out of
setwd("~/2025/I2GDS2025.git/Group1/Needed Materials")

# change filepath
fungi_phylum <- read.csv("~/2025/I2GDS2025.git/Group1/Needed Materials/fungilevel2-table.tsv", sep = "\t") 

# renames all cloumn ncol(fungi_phylum)-1 -> num of column -1 (because first column is taxonomy) "1:" -> creates a column from 1: to the number of column(ncol(fungi_phylum)-1)
colnames(fungi_phylum) <- c("OTU_ID", paste0("fsample", 1:(ncol(fungi_phylum)-1)))


fungi_phylum$Taxonomy <- sapply(fungi_phylum$OTU_ID, function(x) {
  if (grepl("p__", x)) {
    gsub("_", "", sub("^.*p__", "", x))
  } else {
    "Unassigned"
  }
})

sample_cols <- colnames(fungi_phylum)[2:(ncol(fungi_phylum)-1)] 

# convert to numeric
fungi_phylum[sample_cols] <- lapply(fungi_phylum[sample_cols], function(x) as.numeric(as.character(x)))

sum(is.na(fungi_phylum[sample_cols]))

fungi_phylum[sample_cols] <- sweep(fungi_phylum[sample_cols], 2, colSums(fungi_phylum[sample_cols], na.rm = TRUE), FUN = "/")

long_fungi <- pivot_longer(fungi_phylum,
                          cols = all_of(sample_cols),
                          names_to = "Sample",
                          values_to = "Abundance")

head(long_fungi)

plot_stacked <- function(df, output_file, title="Phylum Abundance") {
  ggplot(df, aes(x = Sample, y = Abundance, fill = Taxonomy)) +
    geom_bar(stat = "identity", position = "stack", na.rm = TRUE) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(1, "cm"),
      legend.key.height = unit(1, "cm"),
      legend.key.width = unit(1, "cm"),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    guides(fill = guide_legend(nrow = 5, byrow = TRUE)) +
    labs(title = title, y = "Relative Abundance")
  
  ggsave(output_file, width = 12, height = 8)
}

plot_stacked(long_fungi, "G1-fungi-phylum.png", title="Fungi Phylum Relative Abundance")


