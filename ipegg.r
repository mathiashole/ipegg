#!/usr/bin/env Rscript

# Library
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(tidyverse)
  library(gggenes)
})

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Initialize variables
input_file <- NULL
output_file <- "output.pdf"
ordenar <- FALSE
tree_file <- NULL
normalize_string <- NULL
remove_arg <- NULL
replace_names <- list()

# Parse arguments
for (i in seq_along(args)) {
  if (args[i] == "--input" || args[i] == "-i") {
    input_file <- args[i + 1]
  } else if (args[i] == "--output" || args[i] == "-o") {
    output_file <- args[i + 1]
  } else if (args[i] == "--ordenar") {
    ordenar <- TRUE
  } else if (args[i] == "--tree" || args[i] == "-t") {
    tree_file <- args[i + 1]
  } else if (args[i] == "--normalizeString" || args[i] == "-ns") {
    normalize_string <- args[i + 1]
  } else if (args[i] == "--remove" || args[i] == "-rm") {
    remove_arg <- args[i + 1]
    # remove_words <- unlist(strsplit(args[i + 1], ","))
  } else if (args[i] == "--replace" || args[i] == "-rp") {
    replacement_args <- args[i + 1]
  }
}

# validation
if (is.null(input_file)) {
  cat("Usage: plot_domains.R --input <file.tsv> [--output <file.pdf>] [--sort --contree <.contree>] [--normalize_string word1,word2] [--remove word3,word4] [--replace old=new,...]\n")
  quit(status = 1)
}

# read tsv file
data <- read.delim(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

df_sorted <- data %>%
  group_by(V1) %>%           # Group by column V1
  arrange(V7, .by_group = TRUE) %>% # Sort within each group according to column V7
  ungroup()                  # Remove clustering to avoid unwanted effects later


if (!is.null(normalize_string)) {
  rules <- unlist(strsplit(normalize_string, ","))

  for (rule in rules) {
    parts <- unlist(strsplit(rule, ":", fixed = TRUE))
    if (length(parts) == 2) {
      pattern <- parts[1]
      replacement <- parts[2]
      df_sorted$V5 <- ifelse(
        grepl(pattern, df_sorted$V5, ignore.case = TRUE),
        replacement,
        df_sorted$V5
      )
    }
  }
}

replace_names <- c()
if (!is.null(replacement_args)) {
  replacement_pairs <- unlist(strsplit(replacement_args, ","))
  replace_names <- sapply(replacement_pairs, function(x) {
    kv <- unlist(strsplit(x, "=", fixed = TRUE))
    if (length(kv) == 2) {
      setNames(kv[2], kv[1])
    } else {
      stop(paste("Invalid replace pair:", x))
    }
  }, simplify = FALSE)
  replace_names <- unlist(replace_names)
} # replace (e.g. "PS00022=EGF-like-1,SM00710=PbH1,SSF51126=Pectin lyase-like,SSF57184=Growth factor receptor domain")

#=====================================================================
# Remove words

words_to_remove <- c()
if (!is.null(remove_arg)) {
  words_to_remove <- unlist(strsplit(remove_arg, ","))
}

df_filtered <- df_sorted
if (length(words_to_remove) > 0) {
  df_filtered <- df_sorted %>%
    filter(!grepl(paste(words_to_remove, collapse = "|"), V5, ignore.case = TRUE))
} # Remove (e.g. "SIGNAL_PEPTIDE_H_REGION,SIGNAL_PEPTIDE_C_REGION,SIGNAL_PEPTIDE_N_REGION,TRANSMEMBRANE,NON_CYTOPLASMIC_DOMAIN,CYTOPLASMIC_DOMAIN,PTHR24044,mobidb-lite,PF11024,PF11038,PF11040,PF22274,PF22279,G3DSA:2.160.20.10"

#=====================================================================
# Order tree

if(!is.null(tree_file)) {
  tree_data <- read.tree(tree_file)
  tree_plot <- ggtree(tree_data) +
    geom_tiplab() +
    geom_treescale()

    tree_order <- tree_plot$data  %>%
      filter(isTip)  %>%
      arrange(y)  %>%
      pull(label)
}

#=====================================================================
# Rename and plot

# Renombrar columnas
df_rename <- df_filtered %>%
  rename(
    block_id = V1, 
    domain = V5, 
    start = V7, 
    end = V8, 
    label = V6
  )

nbh <- df_rename %>%
  mutate(strand = "+") %>%  # Agregar columna strand
  select(label, nucleotide = domain, block_id, start, end, strand) %>%
  distinct()

myset <- df_rename %>%
  mutate(from = start, to = end, strand = "+") %>%  # Agregar columna strand
  select(label, block_id, domain, start, end, from, to, strand) %>%
  distinct()

seq_limits <- nbh %>%
  mutate(block_id = factor(block_id)) %>%  # Convert to factor
  group_by(block_id) %>%
  summarize(
    seq_start = min(start),
    seq_end = max(end)
  )

nbh <- nbh %>% mutate(block_id = factor(block_id))
myset <- myset %>% mutate(block_id = factor(block_id))

nbh_plot <- ggplot() +
  # Central line per sequence
  geom_segment(
    data = seq_limits,
    aes(x = seq_start, xend = seq_end,
        y = as.numeric(block_id),  # Same numeric approach
        yend = as.numeric(block_id)),
    color = "gray50",
    linewidth = 0.8,
    alpha = 0.8
  ) +
  # main gene
  geom_rect(
    data = (nbh %>% distinct()),
    aes(xmin = start, xmax = end,
        ymin = as.numeric(block_id) - 0.2,
        ymax = as.numeric(block_id) + 0.2),
    fill = "grey80",
    alpha = 0.7
  ) +
  # Domain
  geom_rect(
    data = myset,
    aes(xmin = from, xmax = to,
        ymin = as.numeric(block_id) - 0.25,
        ymax = as.numeric(block_id) + 0.25,
        fill = domain),
    color = NA
  ) +
  # Scales and themes
  scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous(
    breaks = NULL,
    expand = expansion(add = 0.5)  # Better space
  ) +
  labs(x = "Posición", y = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    legend.position = "right"
  )

output_file_png <- "interproScan_plot.png"
output_file <- "interproScan_plot.pdf"
ggsave(output_file_png, plot = nbh_plot, width = 18, height = 8, dpi = 600)
ggsave(output_file, plot = nbh_plot, width = 18, height = 10)