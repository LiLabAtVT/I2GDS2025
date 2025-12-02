library(tidyverse)
library(data.table)
library(ggrepel)

indir <- "PATH/TO/INPUTDIR"

# List input files
kraken_files <- list.files(indir, pattern = "_kraken_subset_1m.out$", full.names = TRUE)
taxmap_files <- list.files(indir, pattern = "_tax-map.tsv$", full.names = TRUE)
diamond_files <- list.files(indir, pattern = "_subset_1m_assembly.daa$", full.names = TRUE)

# Helper function to get sample name
get_sample_name <- function(filepath) {
  gsub("_assembly_kraken_subset_1m.out$", "", basename(filepath))
}

# Settings
desired_slices <- 10      # total slices including "Other"
other_min_pct <- 10       # minimum percentage for the Other slice

for (kraken_file in kraken_files) {

  sample_name <- get_sample_name(kraken_file)
  message("\nProcessing sample: ", sample_name)

  # Find matching taxmap and diamond files
  taxmap_file <- taxmap_files[grepl(sample_name, taxmap_files)]
  diamond_file <- diamond_files[grepl(sample_name, diamond_files)]

  if (length(taxmap_file) == 0) { message("No taxmap found for sample: ", sample_name); next }
  if (length(diamond_file) == 0) { message("No diamond file found for sample: ", sample_name); next }

  message("Using taxmap: ", basename(taxmap_file))
  message("Using diamond: ", basename(diamond_file))

  #############################################
  #### ---------- KRAKEN PROCESS ---------- ###
  #############################################

  kraken <- read_tsv(
    kraken_file,
    col_names = c("classified", "seq_id", "taxid", "sequence_length", "lca_mapping"),
    show_col_types = FALSE
  )

  kraken_out <- kraken %>%
    filter(classified == "C") %>%
    select(seq_id, taxid) %>%
    arrange(taxid)

  kraken_pie <- kraken_out %>%
    count(taxid, name = "group")

  taxmap <- read.csv(
    taxmap_file,
    header = TRUE,
    sep = "\t"
  ) %>%
    select(taxid, Tax_name)

  # Merge with taxmap
  kraken_pie <- merge(kraken_pie, taxmap, by = "taxid", all.x = TRUE)

  # Filtering
  kraken_pie <- kraken_pie %>%
    filter(group > 5) %>%
    filter(!taxid %in% c(1, 2, 131567, 9606))

  # Dynamic threshold for top slices
  kraken_pie <- kraken_pie %>% arrange(desc(group))
  if (nrow(kraken_pie) > desired_slices) {
    threshold <- kraken_pie$group[desired_slices]

    kraken_pie <- kraken_pie %>%
      mutate(Tax_name = ifelse(group < threshold, "Other", Tax_name)) %>%
      group_by(Tax_name) %>%
      summarise(group = sum(group), .groups = "drop") %>%
      arrange(desc(group))
  }

  # Force minimum size for Other
  total_count <- sum(kraken_pie$group)
  other_index <- which(kraken_pie$Tax_name == "Other")
  if (length(other_index) == 1) {
    current_other_pct <- kraken_pie$group[other_index] / total_count * 100
    remaining_pct <- 100 - other_min_pct
    remaining_total <- total_count - kraken_pie$group[other_index]
    adjust_factor <- remaining_pct / (remaining_total / total_count * 100)

    kraken_pie <- kraken_pie %>%
      mutate(group = ifelse(
        Tax_name == "Other",
        total_count * other_min_pct / 100,
        group * adjust_factor
      ))
  }

  # Compute percentages
  kraken_pie <- kraken_pie %>%
    mutate(
      pct = group / sum(group) * 100,
      pct_label = paste0(round(pct, 1), "%")
    )

  # KRAKEN PLOT
  kraken_plot <- ggplot(kraken_pie, aes(x = "", y = group, fill = Tax_name)) +
    geom_bar(stat = "identity", width = 1, color = "black") +
    coord_polar("y", start = 0) +
    theme_void() +
    labs(
      title = paste("Kraken2 Pie Chart for", sample_name),
      fill = "Species"
    ) +
    geom_label_repel(
      aes(label = pct_label),
      position = position_stack(vjust = 0.5),
      show.legend = FALSE,
      color = "black",
      box.padding = 0.6,
      point.padding = 0.3,
      segment.color = "black",
      segment.size = 0.5
    )

  # Save Kraken plot
  ggsave(
    filename = file.path(indir, paste0("KrakenPie_", sample_name, ".png")),
    plot = kraken_plot,
    width = 7, height = 7, dpi = 300,
    bg = "white"
  )


  #############################################
  #### ---------- DIAMOND PROCESS ---------- ###
  #############################################

  diamond <- fread(
    diamond_file,
    col.names = c(
      "seq_id", "sseqid", "pident", "length", "mismatch", "gapopen",
      "qstart", "qend", "sstart", "send", "evalue", "bitscore",
      "taxid", "species", "diamond_title"
    )
  )

  diamond_out <- diamond %>%
    select(seq_id, taxid, species) %>%
    arrange(desc(taxid))

  diamond_pie <- diamond_out %>%
    count(species, name = "group") %>%
    filter(group > 50) %>%
    filter(species != "N/A") %>%
    arrange(desc(group))

  # Dynamic threshold for top slices
  if (nrow(diamond_pie) > desired_slices) {
    threshold <- diamond_pie$group[desired_slices]

    diamond_pie <- diamond_pie %>%
      mutate(species = ifelse(group < threshold, "Other", species)) %>%
      group_by(species) %>%
      summarise(group = sum(group), .groups = "drop") %>%
      arrange(desc(group))
  }

  # Force minimum size for Other
  total_count <- sum(diamond_pie$group)
  other_index <- which(diamond_pie$species == "Other")
  if (length(other_index) == 1) {
    current_other_pct <- diamond_pie$group[other_index] / total_count * 100
    remaining_pct <- 100 - other_min_pct
    remaining_total <- total_count - diamond_pie$group[other_index]
    adjust_factor <- remaining_pct / (remaining_total / total_count * 100)

    diamond_pie <- diamond_pie %>%
      mutate(group = ifelse(
        species == "Other",
        total_count * other_min_pct / 100,
        group * adjust_factor
      ))
  }

  # Compute percentages
  diamond_pie <- diamond_pie %>%
    mutate(
      pct = group / sum(group) * 100,
      pct_label = paste0(round(pct, 1), "%")
    )

  # DIAMOND PLOT
  diamond_plot <- ggplot(diamond_pie, aes(x = "", y = group, fill = species)) +
    geom_bar(stat = "identity", width = 1, color = "black") +
    coord_polar("y", start = 0) +
    theme_void() +
    labs(
      title = paste("DIAMOND Pie Chart for", sample_name),
      fill = "Species"
    ) +
    geom_label_repel(
      aes(label = pct_label),
      position = position_stack(vjust = 0.5),
      show.legend = FALSE,
      color = "black",
      box.padding = 0.6,
      point.padding = 0.3,
      segment.color = "black",
      segment.size = 0.5
    )

  # Save DIAMOND plot
  ggsave(
    filename = file.path(indir, paste0("DiamondPie_", sample_name, ".png")),
    plot = diamond_plot,
    width = 7, height = 7, dpi = 300,
    bg = "white"
  )

  message("Finished sample: ", sample_name)
}