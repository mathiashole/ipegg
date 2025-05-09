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


# If there are multiple --normalize, combine them
if (!is.null(normalize_string)) {
  normalization_rules <- strsplit(normalize_string, ",")[[1]]
  normalization_rules <- lapply(normalize_string, function(x) {
    parts <- strsplit(x, ":", fixed = TRUE)[[1]]
    list(pattern = parts[1], replacement = parts[2])
  })

  for (rule in normalization_rules) {
    df_sorted$V5 <- ifelse(
      grepl(rule$pattern, df_sorted$V5, ignore.case = TRUE),
      rule$replacement,
      df_sorted$V5
    )
  }
} # Normalization string (e.g. "SIGNAL_PEPTIDE:Signal_peptide,SignalP-TM:Signal_peptide")

if (length(normalization_rules) > 0) {
  for (rule in normalization_rules) {
    df_sorted$V5 <- ifelse(
      grepl(rule$pattern, df_sorted$V5, ignore.case = TRUE),
      rule$replacement,
      df_sorted$V5
    )
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

# Generar el gráfico
nbh_plot <- ggplot(
  (nbh %>% distinct()),
  aes(xmin = start, xmax = end, y = block_id)
) +
  geom_gene_arrow() +
  geom_subgene_arrow(
    data = myset,
    aes(xmin = start, xmax = end, xsubmin = from, xsubmax = to, fill = domain, y = block_id)
  ) +
  geom_subgene_label(
    data = myset,
    aes(xmin = start, xmax = end, xsubmin = from, xsubmax = to, label = domain, y = block_id),
    min.size = 0
  ) +
  scale_fill_brewer(palette = "Dark2") +
  theme_genes() %+replace% 
  theme(
    panel.grid.major.y = element_line(colour = NULL),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  )

output_file_png <- "interproScan_plot.png"
output_file <- "interproScan_plot.pdf"
ggsave(output_file_png, plot = nbh_plot, width = 18, height = 8, dpi = 600)
ggsave(output_file, plot = nbh_plot, width = 18, height = 10)