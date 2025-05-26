#!/usr/bin/env Rscript

# Library
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(tidyverse)
})

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Initialize variables
input_file <- NULL
output_file <- NULL
ordenar <- FALSE
tree_file <- NULL
normalize_string <- NULL
remove_arg <- NULL
backgroud_arg <- NULL
itol_data <- FALSE

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
  } else if (args[i] == "--replace" || args[i] == "-rp") {
    normalize_string <- args[i + 1]
  } else if (args[i] == "--remove" || args[i] == "-rm") {
    remove_arg <- args[i + 1]
    # remove_words <- unlist(strsplit(args[i + 1], ","))
  } else if (args[i] == "--background" || args[i] == "-bg") {
    background_arg <- args[i + 1]
  } else if (args[i] == "--itol" || args[i] == "-it") {
    itol_data <- TRUE
  } # else if (args[i] == "--replace" || args[i] == "-rp") {
  #   replacement_args <- args[i + 1]
  # }
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
} ## Normalized (e.g. --replace "SIGNAL_PEPTIDE:Signal_peptide,SignalP-TM:Signal_peptide,SignalP-noTM:Signal_peptide,PS00022:EGF-like-1,SM00710:PbH1,SSF51126:Pectin-lyase-like,SSF57184:GF receptor")
# replace (e.g. "PS00022=EGF-like-1,SM00710=PbH1,SSF51126=Pectin lyase-like,SSF57184=Growth factor receptor domain")

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
    linewidth = 0.5,
    alpha = 0.8
  )

# if(!is.null(shadow_domain)) {
#   # Shadow structure
#   nbh_plot <- nbh_plot +
#     geom_rect(
#       data = (nbh %>% distinct()),
#       aes(xmin = start, xmax = end,
#           ymin = as.numeric(block_id) - 0.2,
#           ymax = as.numeric(block_id) + 0.2),
#       fill = "grey80",
#       alpha = 0.7
#     )
# }

nbh_plot <- nbh_plot +
  geom_rect(
    data = myset,
    aes(xmin = from, xmax = to,
        ymin = as.numeric(block_id) - 0.35,
        ymax = as.numeric(block_id) + 0.35,
        fill = domain),
    color = NA
  ) +
  scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous(
    breaks = NULL,
    expand = expansion(add = 0.5)
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

if (is.null(output_file)) {
  base_name <- sub("\\.[^.]*$", "", basename(input_file))
  base_name <- paste0(base_name, "_interproScan_plot")
} else {
  base_name <- sub("\\.[^.]*$", "",output_file)
}

output_file_png <- paste0(base_name, ".png")
output_file_pdf <- paste0(base_name, ".pdf")
output_file_itol <- paste0(base_name, ".itol")

cat("Input:", input_file, "\n")
cat("Output:", output_file_png, "\n")
cat("Output:", output_file_pdf, "\n")

ggsave(output_file_png, plot = nbh_plot, width = 18, height = 8, dpi = 600)
ggsave(output_file_pdf, plot = nbh_plot, width = 18, height = 10)

#======================================================================
# itol


convert_to_itol <- function(data, output_file = "itol_domains.txt", dataset_label = "Protein Domains") {
  # Prepare the data: sort and calculate protein lengths
  data_prepared <- data %>%
    arrange(block_id, start) %>%
    group_by(block_id) %>%
    mutate(protein_length = max(end)) %>%
    ungroup()
  
  # Assign unique colors to each domain (consistent with "Dark2" Brewer palette)
  domain_levels <- levels(factor(data_prepared$domain))
  domain_colors <- setNames(
    brewer.pal(max(3, length(domain_levels)), "Dark2")[seq_along(domain_levels)],
    domain_levels
  )
  
  # Create the file header
  header <- c(
    "DATASET_DOMAINS",
    "#Protein domain datasets are visualized as schematic representations of proteins",
    "",
    "#=================================================================#",
    "#                    MANDATORY SETTINGS                           #",
    "#=================================================================#",
    "SEPARATOR COMMA",
    paste0("DATASET_LABEL,", dataset_label),
    "COLOR,#ff0000",  # Default color (can be changed later in iTOL)
    "",
    "#=================================================================#",
    "#                    OPTIONAL SETTINGS                            #",
    "#=================================================================#",
    "BACKBONE_COLOR,#dddddd",
    "BACKBONE_HEIGHT,8",
    "SHOW_DOMAIN_LABELS,1",
    "LABEL_SIZE_FACTOR,0.8",
    "LABEL_AUTO_COLOR,1",
    "BORDER_WIDTH,0.5",
    "BORDER_COLOR,#000000",
    "",
    "#=================================================================#",
    "#       Actual data follows after the \"DATA\" keyword              #",
    "#=================================================================#",
    "DATA",
    ""
  )
  
  # Generate iTOL domain lines
  data_lines <- data_prepared %>%
    group_by(block_id, protein_length) %>%
    summarize(
      domains = paste(
        sprintf("RE|%d|%d|%s|%s", start, end, domain_colors[domain], domain),
        collapse = ","
      ),
      .groups = "drop"
    ) %>%
    mutate(itol_line = sprintf("%s,%d,%s", block_id, protein_length, domains)) %>%
    pull(itol_line)
  
  # Write the output file
  writeLines(c(header, data_lines), output_file)
  
  message("iTOL file successfully generated at: ", output_file)
  return(invisible(domain_colors))
}

if(itol_data) {
  # Initial message
  message("\n🔵 [1/3] Preparing data for iTOL...")

  convert_to_itol(myset, output_file = "my_domains.itol", dataset_label = "My Protein Domains")
}