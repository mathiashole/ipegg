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
keywords <- NULL
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
  } else if (args[i] == "--keywords" || args[i] == "-k") {
    keywords <- unlist(strsplit(args[i + 1], ","))
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
  cat("Uso: plot_dominios.R --input <archivo.tsv> [--output <archivo.pdf>] [--ordenar --contree <.contree>] [--keywords palabra1,palabra2] [--remove palabra3,palabra4] [--replace vieja=nueva,...]\n")
  quit(status = 1)
}

# read tsv file
dominios <- read.delim(input_file, header = TRUE, stringsAsFactors = FALSE)