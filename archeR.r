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
remove_words <- NULL
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
  } else if (args[i] == "--remove" || args[i] == "-r") {
    remove_words <- unlist(strsplit(args[i + 1], ","))
  } else if (args[i] == "--replace" || args[i] == "-rp") {
    replacements <- unlist(strsplit(args[i + 1], ","))
    replace_names <- sapply(replacements, function(x) {
      kv <- unlist(strsplit(x, "="))
      setNames(kv[2], kv[1])
    }, simplify = FALSE)
  }
}

# validation
if (is.null(input_file)) {
  cat("Usage: plot_domains.R --input <file.tsv> [--output <file.pdf>] [--sort --contree <.contree>] [--normalize_string word1,word2] [--remove word3,word4] [--replace old=new,...]\n")
  quit(status = 1)
}

# read tsv file
dominios <- read.delim(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# If there are multiple --normalize, combine them
if (length(normalize_patterns) > 0) {
  normalization_rules <- lapply(normalize_patterns, function(x) {
    parts <- strsplit(x, ":", fixed = TRUE)[[1]]
    list(pattern = parts[1], replacement = parts[2])
  })
} else {
  normalization_rules <- list()
}