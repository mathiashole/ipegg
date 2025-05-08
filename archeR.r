#!/usr/bin/env Rscript

# Library
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
})

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Initialize variables
input_file <- NULL
output_file <- "output.pdf"
ordenar <- FALSE
contree_file <- NULL
normalize_string <- NULL
remove_args <- NULL
replace_names <- list()

# Parse arguments
for (i in seq_along(args)) {
  if (args[i] == "--input" || args[i] == "-i") {
    input_file <- args[i + 1]
  } else if (args[i] == "--output" || args[i] == "-o") {
    output_file <- args[i + 1]
  } else if (args[i] == "--ordenar") {
    ordenar <- TRUE
  } else if (args[i] == "--contree" || args[i] == "-t") {
    contree_file <- args[i + 1]
  } else if (args[i] == "--normalizeString" || args[i] == "-ns") {
    normalize_string <- args[i + 1]
  } else if (args[i] == "--remove" || args[i] == "-rm") {
    remove_args <- args[i + 1]
    # remove_words <- unlist(strsplit(args[i + 1], ","))
  } else if (args[i] == "--replace" || args[i] == "-rp") {
    replacement_args <- args[i + 1]
    # replacements <- unlist(strsplit(args[i + 1], ","))
    # replace_names <- sapply(replacements, function(x) {
    #   kv <- unlist(strsplit(x, "="))
    #   setNames(kv[2], kv[1])
    # }, simplify = FALSE)
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