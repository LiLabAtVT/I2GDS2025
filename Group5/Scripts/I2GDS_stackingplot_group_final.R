# Load packages
library(tidyverse)
library(RColorBrewer)

# Print all packages and versions used for the script
print(sessionInfo())

# Set working directory
setwd("/Users/lihong/Desktop/bracken_data")

# Function to read bracken files
read_bracken <- function(filename) {
  # Get sample id by removing file extension from file name
  sample_id <- str_remove(filename, "\\.bracken$")
  # Read the Bracken file
  read_tsv(filename, show_col_types = FALSE) %>%
  # Filter to only genus level taxonomy 
  #   - Don't technically need this for our data, but here for if input changes
  filter(taxonomy_lvl == "G") %>%
  # Select only the columns with genus name and abundance value
  select(Genus = name, Abundance = fraction_total_reads) %>%
  # Add the sample id as a column
  mutate(Sample = sample_id)
}

### Read and merge all bracken files ###

# Get list of file names
bracken_files <- list.files(pattern = "\\.bracken$")
# Print how many files were found
cat("Found", length(bracken_files), "bracken files\n")

# Read and filter all files using function and merge them into one long dataframe
all_data <- map_df(bracken_files, read_bracken)
# Write merged data to a file
write_tsv(all_data, "all_genus_abundance.txt")+
cat("Data saved to: all_genus_abundance.txt\n\n")

# Function to process data: keep top N genera, group others
process_and_order <- function(data, top_n) {
  # Get top genera by mean abundance
  top_genera <- data %>%
    # Group data by Genus and calculate mean abundance for each
    group_by(Genus) %>% summarise(mean_abund = mean(Abundance)) %>%
    # Filter to the top 20
    slice_max(mean_abund, n = top_n) %>%
    # Pull just the species names
    pull(Genus)
  
  # Group data and order samples by Prevotella abundance
  data_processed <- data %>%
    # Make a new column of the original data that keeps the genus name if it is 
    # in the top 20 or labels it as "Other" if it isn't
    mutate(Genus_plot = if_else(Genus %in% top_genera, Genus, "Other")) %>%
    # Sum abundance values for each new genus label within each sample to group all "Other"
    group_by(Sample, Genus_plot) %>% summarise(Abundance = sum(Abundance), .groups = "drop")
  
  # Order samples
  sample_order <- data_processed %>%
    # Filter data based on Prevotella abundance
    filter(Genus_plot == "g__Prevotella") %>%
    # Order samples based on descending Prevotella abundance
    arrange(desc(Abundance)) %>%
    # Save only the sample id column
    pull(Sample)
  
  # Use list of sample ids to reorder the full dataframe
  data_processed %>%
    mutate(Sample = factor(Sample, levels = sample_order)) %>%
    arrange(Sample)
}

# Use function to process the data
data_final <- process_and_order(all_data, top_n = 20)

# Set colors
# Count the number of genera
n_genera <- length(unique(data_final$Genus_plot))
# Get the same number of colors as we have genera by interpolating the Paired 
# color palette that has 12 colors in it
colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_genera)
# Assign each genus to a color
names(colors) <- sort(unique(data_final$Genus_plot))
# Change the color of the Other category to gray
if ("Other" %in% names(colors)) colors["Other"] <- "#CCCCCC"

# Create plot
# Using final filtered and ordered data, set x-axis to be the samples, set the
# y-axis to be the abundance value times 100 to make it a percent, and set the fill
# (in this case the bars) to be the genera
p <- ggplot(data_final, aes(x = Sample, y = Abundance * 100, fill = Genus_plot)) +
  # Make a bar plot using the y value specified from our dataframe and set the 
  # width of the bars equal to one
  geom_bar(stat = "identity", width = 1) +
  # Set the colors based on the list of colors created above
  scale_fill_manual(values = colors, name = "Genus") +
  # Set the min and mix for the y-axis to be 0 and 100 respectively
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  # Set axis labels and title
  labs(x = "Sample (Participant)", y = "Relative Abundance (%)",
       title = "Gut Microbiome Composition at Genus Level") +
  # Use simple theme
  theme_minimal() +
  theme(
    # Set the x-axis labels to angle 45 degrees
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    # Put the legend on the right hand side of the plot
    legend.position = "right",
    # Set the box size for the legend to be half a cm
    legend.key.size = unit(0.5, "cm"),
    # Adjust sizing of, shift placement of, and bold title
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    # Bold and adjust sizing of axis labels
    axis.title = element_text(face = "bold", size = 12)
  )

# Print for viewing in RStudio
print(p)
# Save the plot to a pdf file with set height and width
ggsave("genus_abundance_plot.pdf", p, width = 12, height = 6)
# Print out where the plot has been saved
cat("Plot saved to: genus_abundance_plot.pdf\n")

