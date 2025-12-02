#Loading packages
library(tidyverse)
library(RColorBrewer)

#Set working directory - change yourself!
setwd("C:/Users/Jvlan/Downloads")

#Reading in all data needed
#The drugs object comes from the deeparg database, where they assigned each specific antibiotic to a broader class
#The CARD_aro object is from the CARD database, and shows which gene confers resistance to each specific antibiotic
CARD_aro <- read.csv("CARD4.0.1_aro_cat.csv")
Args <- read.csv("I2GDS_G6_finalDiamond.csv")
metadata <- read.csv("I2GDSMetadata.csv")
drugs <- read.csv("Antibiotic to Drug.csv")

#To make the graph look better, we assign the genes in the CARD database to a broader drug class using the deeparg conversion
for (i in 1:nrow(CARD_aro)) {
  for (j in 1:nrow(drugs)) {
    if(grepl(drugs[j,1], CARD_aro[i, 4])) {CARD_aro[i, "Drug"] <- drugs[j,2]}
  }
}

for (i in 1:nrow(CARD_aro)) {
  if(grepl(";", CARD_aro[i, 4])) {CARD_aro[i, "Drug"] <- "Multidrug"}
}

CARD_aro$Drug[is.na(CARD_aro$Drug)] <- "Other"

#Merging the DIAMOND data and the CARD metadata about each gene
Args <- merge(Args, CARD_aro, by.x = "protein_accession", by.y = "Protein.Accession")

#Renaming the samples to only include the same name
Args$sample = substr(Args$sample, 1, nchar(Args$sample)-9)

#Merging the DIAMOND data and sample metadata
Args <- merge(Args, metadata, by.x = "sample", by.y = "Run")

Args$rpob_Normalization[Args$rpob_Normalization==""] <- 0

#Assigning data its appropriate type
Args$Drug.Class <- as.factor(Args$Drug.Class)
Args$isolation_source <- as.factor(Args$isolation_source)
Args$count <- as.numeric(Args$count)
Args$rpob_Normalization <- as.numeric(Args$rpob_Normalization)
Args$Drug <- as.factor(Args$Drug)

#Grouping the DIAMOND data to find how often each type of resistance was in each sample
df <- Args %>%
  group_by(sample, isolation_source, Drug) %>%
  summarize(sum_count = sum(count))

#Checking most abundant drug class
a <- df %>%
  group_by(Drug) %>%
  summarize(drug_count = sum(sum_count))

#Reassigning the orders to make the graph look better. Drug was done in order of highest to lowest counts, isolation source to match the paper
df$Drug <- factor(df$Drug, levels = c("Multidrug", "Other", "tetracycline", "macrolide-lincosamide-streptogramin", "aminoglycoside", "aminocoumarin",
                                      "peptide", "sulfonamide", "quinolone", "beta_lactam",
                                      "phenicol", "mupirocin", "pleuromutilin"))

df$isolation_source <- factor(df$isolation_source, levels = c("Influent", "Primary", "Act. Sludge",
                                                              "Secondary", "Floc Sed", "Ozonation",
                                                              "BAC/GAC", "Chlorination", "3 days", "Background"))


#Creating an object so we don't have to use R's ugly default colors
n <- 13
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#Actually creating the stacked bar chart
ggplot(df, aes(x = isolation_source, y = sum_count, fill = Drug)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values=col_vector) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Total ARG copies") +
  xlab("")
