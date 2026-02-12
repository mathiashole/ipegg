#!/usr/bin/env Rscript

# Library
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(tidyverse)
  library(yaml)
})

# Read configuration from YAML file
#--------------------------------------------------------------------------

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2 || args[1] != "--config") {
  stop("Usage: plot_domains_yaml.R --config config.yaml")
}

config_file <- args[2]

if (!file.exists(config_file)) {
  stop("Error: Configuration file not found: ", config_file)
}

config <- yaml.load_file(config_file)

# Assign config values to variables
#-------------------------------------------------------------------------

input_file <- config$input$file
output_pref <- config$output$prefix

tree_file <- config$tree$file # Optional: Path to the .contree file for ordering sequences

ordenar        <- isTRUE(config$options$ordenar)
itol_data      <- isTRUE(config$options$itol)
generate_stats <- isTRUE(config$options$statistics)

replace_rules  <- config$normalize$replace
words_to_remove <- config$remove

itol_label <- config$itol$dataset_label %||% "Protein Domains"

if (is.null(input_file)) {
  stop("Input file not defined in config.yaml")
}

# read tsv file
data <- read.delim(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

df_sorted <- data %>%
  group_by(V1) %>% # Group by column V1
  arrange(V7, .by_group = TRUE) %>% # Sort within each group according to column V7
  ungroup() # Remove clustering to avoid unwanted effects later

# Normalize domain names based on replace rules
#---------------------------------------------------------------------------

if (!is.null(replace_rules)) {
  message("Normalizing domain names...")
  for (pattern in names(replace_rules)) {
    replacement <- replace_rules[[pattern]] # Get the replacement value for the current pattern
    df_sorted$V5 <- ifelse( # Use ifelse to replace values in V5 based on the pattern
      grepl(pattern, df_sorted$V5, ignore.case = TRUE),
      replacement,
      df_sorted$V5
    )
  }
}

# Remove specified words from domain names
#---------------------------------------------------------------------------

df_filtered <- df_sorted

if (!is.null(words_to_remove) && length(words_to_remove) > 0) {
  message("Removing unwanted domains...")
  df_filtered <- df_sorted %>%
    filter(!grepl(paste(words_to_remove, collapse = "|"),
                  V5, ignore.case = TRUE))
}

# Tree ordering (optional)
#---------------------------------------------------------------------------

if (!is.null(tree_file)) {
  message("Reading tree file...")
  tree_data <- read.tree(tree_file)
  tree_plot <- ggtree(tree_data)
  
  tree_order <- tree_plot$data %>%
    filter(isTip) %>%
    arrange(y) %>%
    pull(label)
}

# Rename columns for clarity
#---------------------------------------------------------------------------

df_rename <- df_filtered %>%
  rename(
    block_id = V1, 
    domain = V5, 
    start = V7, 
    end = V8, 
    label = V6
  )

# Prepare data for plotting
#---------------------------------------------------------------------------

nbh <- df_rename %>%
  mutate(strand = "+") %>%  # Add strand column
  select(label, nucleotide = domain, block_id, start, end, strand) %>%
  distinct() # Ensure distinct rows

myset <- df_rename %>%
  mutate(from = start, to = end, strand = "+") %>%  # Add strand column
  select(label, block_id, domain, start, end, from, to, strand) %>%
  distinct() # Ensure distinct rows

# Background domain plotting
#---------------------------------------------------------------------------

seq_limits <- nbh %>%
  mutate(block_id = factor(block_id)) %>%  # Convert to factor
  group_by(block_id) %>%
  summarize(seq_start = min(start), seq_end = max(end))

nbh <- nbh %>% mutate(block_id = factor(block_id))
myset <- myset %>% mutate(block_id = factor(block_id))

x_max_real <- max(seq_limits$seq_end, na.rm = TRUE)
x_max_round <- round(x_max_real / 100) * 100
my_breaks <- seq(0, x_max_round, length.out = 10)

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
  scale_x_continuous(
      breaks = my_breaks,
      labels = round(my_breaks),
      limits = c(0, x_max_real),
      expand = expansion(mult = c(0, 0.05)) 
  ) +
  labs(x = "Position", y = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    legend.position = "right"
  )  

# Output file generation
#---------------------------------------------------------------------------

base_name <- ifelse(is.null(output_pref)) {
  output_pref
} else {
  sub("\\.[^.]*$", "", basename(input_file))
}

output_png <- paste0(base_name, ".png")
output_pdf <- paste0(base_name, ".pdf")
output_itol <- paste0(base_name, ".itol")

ggsave(output_png, plot = nbh_plot, width = 18, height = 8, dpi = 600)
ggsave(output_pdf, plot = nbh_plot, width = 18, height = 10)

message("Plot saved: ", output_png)
message("Plot saved: ", output_pdf)

# Descriptive statistics generation (optional)
#---------------------------------------------------------------------------

if (generate_stats) {
  message("Generating domain statistics...")

  df_stats <- df_rename %>%
    mutate(domain_length = end - start + 1)
  # Table 1: Statistics by domain type
  domain_stats <- df_stats %>%
    group_by(domain) %>%
    summarize(
      count = n(),
      mean_length = round(mean(domain_length), 2),
      sd_length = round(sd(domain_length), 2)
    ) %>%
    arrange(desc(count))

  write.table(domain_stats, paste0(base_name, "_domain_stats.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  # Table 2: Domain count per sequence
  domain_per_seq <- df_stats %>%
    group_by(block_id) %>%
    summarize(num_domains = n())

  write.table(domain_per_seq, paste0(base_name, "_domain_count_per_sequence.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    message("📄 Domain count per sequence written to: ", paste0(base_name, "_domain_count_per_sequence.tsv"))

}

# ItoL file generation (optional)
#---------------------------------------------------------------------------

convert_to_itol <- function(data, output_file, dataset_label) {
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

# # Initialize variables
# input_file <- NULL
# output_file <- NULL
# ordenar <- FALSE
# tree_file <- NULL
# normalize_string <- NULL
# remove_arg <- NULL
# backgroud_arg <- NULL
# itol_data <- FALSE
# generate_stats <- FALSE

# # Parse arguments
# for (i in seq_along(args)) {
#   if (args[i] == "--input" || args[i] == "-i") {
#     input_file <- args[i + 1]
#   } else if (args[i] == "--output" || args[i] == "-o") {
#     output_file <- args[i + 1]
#   } else if (args[i] == "--ordenar") {
#     ordenar <- TRUE
#   } else if (args[i] == "--tree" || args[i] == "-t") {
#     tree_file <- args[i + 1]
#   } else if (args[i] == "--replace" || args[i] == "-rp") {
#     normalize_string <- args[i + 1]
#   } else if (args[i] == "--remove" || args[i] == "-rm") {
#     remove_arg <- args[i + 1]
#     # remove_words <- unlist(strsplit(args[i + 1], ","))
#   } else if (args[i] == "--background" || args[i] == "-bg") {
#     background_arg <- args[i + 1]
#   } else if (args[i] == "--itol" || args[i] == "-it") {
#     itol_data <- TRUE
#   } else if (args[i] == "--statistics" || args[i] == "-stat") {
#     generate_stats <- TRUE
#   }# else if (args[i] == "--replace" || args[i] == "-rp") {
#   #   replacement_args <- args[i + 1]
#   # }
# }

# # validation
# if (is.null(input_file)) {
#   cat("Usage: plot_domains.R --input <file.tsv> [--output <file.pdf>] [--sort --contree <.contree>] [--normalize_string word1,word2] [--remove word3,word4] [--replace old=new,...]\n")
#   quit(status = 1)
# }

# # read tsv file
# data <- read.delim(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# df_sorted <- data %>%
#   group_by(V1) %>%           # Group by column V1
#   arrange(V7, .by_group = TRUE) %>% # Sort within each group according to column V7
#   ungroup()                  # Remove clustering to avoid unwanted effects later


# if (!is.null(normalize_string)) {
#   rules <- unlist(strsplit(normalize_string, ","))

#   for (rule in rules) {
#     parts <- unlist(strsplit(rule, ":", fixed = TRUE))
#     if (length(parts) == 2) {
#       pattern <- parts[1]
#       replacement <- parts[2]
#       df_sorted$V5 <- ifelse(
#         grepl(pattern, df_sorted$V5, ignore.case = TRUE),
#         replacement,
#         df_sorted$V5
#       )
#     }
#   }
# } ## Normalized (e.g. --replace "SIGNAL_PEPTIDE:Signal_peptide,SignalP-TM:Signal_peptide,SignalP-noTM:Signal_peptide,PS00022:EGF-like-1,SM00710:PbH1,SSF51126:Pectin-lyase-like,SSF57184:GF receptor")
# # replace (e.g. "PS00022=EGF-like-1,SM00710=PbH1,SSF51126=Pectin lyase-like,SSF57184=Growth factor receptor domain")

# #=====================================================================
# # Remove words

# words_to_remove <- c()
# if (!is.null(remove_arg)) {
#   words_to_remove <- unlist(strsplit(remove_arg, ","))
# }

# df_filtered <- df_sorted
# if (length(words_to_remove) > 0) {
#   df_filtered <- df_sorted %>%
#     filter(!grepl(paste(words_to_remove, collapse = "|"), V5, ignore.case = TRUE))
# } # Remove (e.g. "SIGNAL_PEPTIDE_H_REGION,SIGNAL_PEPTIDE_C_REGION,SIGNAL_PEPTIDE_N_REGION,TRANSMEMBRANE,NON_CYTOPLASMIC_DOMAIN,CYTOPLASMIC_DOMAIN,PTHR24044,mobidb-lite,PF11024,PF11038,PF11040,PF22274,PF22279,G3DSA:2.160.20.10"

# #=====================================================================
# # Order tree

# if(!is.null(tree_file)) {
#   tree_data <- read.tree(tree_file)
#   tree_plot <- ggtree(tree_data) +
#     geom_tiplab() +
#     geom_treescale()

#     tree_order <- tree_plot$data  %>%
#       filter(isTip)  %>%
#       arrange(y)  %>%
#       pull(label)
# }

# #=====================================================================
# # Rename and plot

# # Rename colums
# df_rename <- df_filtered %>%
#   rename(
#     block_id = V1, 
#     domain = V5, 
#     start = V7, 
#     end = V8, 
#     label = V6
#   )

# nbh <- df_rename %>%
#   mutate(strand = "+") %>%  # Agregar columna strand
#   select(label, nucleotide = domain, block_id, start, end, strand) %>%
#   distinct()

# myset <- df_rename %>%
#   mutate(from = start, to = end, strand = "+") %>%  # Agregar columna strand
#   select(label, block_id, domain, start, end, from, to, strand) %>%
#   distinct()

# #=====================================================================
# # Background domain 

# seq_limits <- nbh %>%
#   mutate(block_id = factor(block_id)) %>%  # Convert to factor
#   group_by(block_id) %>%
#   summarize(
#     seq_start = min(start),
#     seq_end = max(end)
#   )

# nbh <- nbh %>% mutate(block_id = factor(block_id))
# myset <- myset %>% mutate(block_id = factor(block_id))

# nbh_plot <- ggplot() +
#   # Central line per sequence
#   geom_segment(
#     data = seq_limits,
#     aes(x = seq_start, xend = seq_end,
#         y = as.numeric(block_id),  # Same numeric approach
#         yend = as.numeric(block_id)),
#     color = "gray50",
#     linewidth = 0.5,
#     alpha = 0.8
#   )

# nbh_plot <- nbh_plot +
#   geom_rect(
#     data = myset,
#     aes(xmin = from, xmax = to,
#         ymin = as.numeric(block_id) - 0.35,
#         ymax = as.numeric(block_id) + 0.35,
#         fill = domain),
#     color = NA
#   ) +
#   scale_fill_brewer(palette = "Dark2") +
#   scale_y_continuous(
#     breaks = NULL,
#     expand = expansion(add = 0.5)
#   ) +
#   labs(x = "Posición", y = NULL) +
#   theme_minimal() +
#   theme(
#     panel.grid = element_blank(),
#     axis.text.y = element_blank(),
#     axis.ticks = element_blank(),
#     panel.background = element_blank(),
#     legend.position = "right"
#   )

# if (is.null(output_file)) {
#   base_name <- sub("\\.[^.]*$", "", basename(input_file))
#   base_name <- paste0(base_name, "_interproScan_plot")
# } else {
#   base_name <- sub("\\.[^.]*$", "",output_file)
# }

# output_file_png <- paste0(base_name, ".png")
# output_file_pdf <- paste0(base_name, ".pdf")
# output_file_itol <- paste0(base_name, ".itol")

# cat("Input:", input_file, "\n")
# cat("Output:", output_file_png, "\n")
# cat("Output:", output_file_pdf, "\n")

# ggsave(output_file_png, plot = nbh_plot, width = 18, height = 8, dpi = 600)
# ggsave(output_file_pdf, plot = nbh_plot, width = 18, height = 10)

# violin_png <- paste0(base_name, "_violin_plot.png")
# violin_pdf <- paste0(base_name, "_violin_plot.pdf")
# stats_file <- paste0(base_name, "_domain_stats.tsv")
# count_file <- paste0(base_name, "_domain_count_per_sequence.tsv")

# #=====================================================================
# # Descriptive statistics (--statistics)

# if (generate_stats) {
#   message("📊 Generating domain statistics...")

#   df_stats <- df_rename %>%
#     mutate(domain_length = end - start + 1)

#   # Table 1: Statistics by domain type
#   domain_stats <- df_stats %>%
#     group_by(domain) %>%
#     summarize(
#       count = n(),
#       mean_length = round(mean(domain_length), 2),
#       sd_length = round(sd(domain_length), 2)
#     ) %>%
#     arrange(desc(count))

#   write.table(domain_stats, stats_file, sep = "\t", row.names = FALSE, quote = FALSE)
#   message("📄 Domain stats written to: ", stats_file)

#   # Table 2: Domain count per sequence
#   domain_per_seq <- df_stats %>%
#     group_by(block_id) %>%
#     summarize(num_domains = n())

#   write.table(domain_per_seq, count_file, sep = "\t", row.names = FALSE, quote = FALSE)
#   message("📄 Domain count per sequence written to: ", count_file)

#   # Violin plot
#   violin_plot <- ggplot(df_stats, aes(x = reorder(domain, domain_length, FUN = median), y = domain_length)) +
#     geom_violin(fill = "#69b3a2", alpha = 0.7, scale = "width") +
#     geom_jitter(height = 0, width = 0.2, alpha = 0.3, size = 0.3) +
#     labs(
#       x = "Domain",
#       y = "Domain length (bp/aa)",
#       title = "Domain length distribution"
#     ) +
#     theme_minimal(base_size = 12) +
#     theme(
#       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#       plot.title = element_text(face = "bold", size = 14)
#     )

#   ggsave(violin_png, violin_plot, width = 12, height = 6, dpi = 300)
#   ggsave(violin_pdf, violin_plot, width = 12, height = 6)
#   message("🎻 Violin plot saved to: ", violin_png, " and ", violin_pdf)
# }

# #======================================================================
# # itol


# convert_to_itol <- function(data, output_file = "itol_domains.txt", dataset_label = "Protein Domains") {
#   # Prepare the data: sort and calculate protein lengths
#   data_prepared <- data %>%
#     arrange(block_id, start) %>%
#     group_by(block_id) %>%
#     mutate(protein_length = max(end)) %>%
#     ungroup()
  
#   # Assign unique colors to each domain (consistent with "Dark2" Brewer palette)
#   domain_levels <- levels(factor(data_prepared$domain))
#   domain_colors <- setNames(
#     brewer.pal(max(3, length(domain_levels)), "Dark2")[seq_along(domain_levels)],
#     domain_levels
#   )
  
#   # Create the file header
#   header <- c(
#     "DATASET_DOMAINS",
#     "#Protein domain datasets are visualized as schematic representations of proteins",
#     "",
#     "#=================================================================#",
#     "#                    MANDATORY SETTINGS                           #",
#     "#=================================================================#",
#     "SEPARATOR COMMA",
#     paste0("DATASET_LABEL,", dataset_label),
#     "COLOR,#ff0000",  # Default color (can be changed later in iTOL)
#     "",
#     "#=================================================================#",
#     "#                    OPTIONAL SETTINGS                            #",
#     "#=================================================================#",
#     "BACKBONE_COLOR,#dddddd",
#     "BACKBONE_HEIGHT,8",
#     "SHOW_DOMAIN_LABELS,1",
#     "LABEL_SIZE_FACTOR,0.8",
#     "LABEL_AUTO_COLOR,1",
#     "BORDER_WIDTH,0.5",
#     "BORDER_COLOR,#000000",
#     "",
#     "#=================================================================#",
#     "#       Actual data follows after the \"DATA\" keyword              #",
#     "#=================================================================#",
#     "DATA",
#     ""
#   )
  
#   # Generate iTOL domain lines
#   data_lines <- data_prepared %>%
#     group_by(block_id, protein_length) %>%
#     summarize(
#       domains = paste(
#         sprintf("RE|%d|%d|%s|%s", start, end, domain_colors[domain], domain),
#         collapse = ","
#       ),
#       .groups = "drop"
#     ) %>%
#     mutate(itol_line = sprintf("%s,%d,%s", block_id, protein_length, domains)) %>%
#     pull(itol_line)
  
#   # Write the output file
#   writeLines(c(header, data_lines), output_file)
  
#   message("iTOL file successfully generated at: ", output_file)
#   return(invisible(domain_colors))
# }

# if(itol_data) {
#   # Initial message
#   message("\n🔵 [1/3] Preparing data for iTOL...")
#   # Converting to iTOL format
#   convert_to_itol(myset, output_file = output_file_itol, dataset_label = "My Protein Domains")

#     # Verification and success message
#   if(file.exists(output_file_itol)) {
#     message(paste0(
#       "✅ [2/3] iTOL file generated successfully:\n",
#       "   - Path: ", normalizePath(output_file_itol), "\n",
#       "   - Dataset: 'My Protein Domains'"
#     ))
#     message("🔵 [3/3] You can upload this file directly to iTOL Web:")
#     message("   🌐 https://itol.embl.de/upload.cgi")
#   } else {
#     warning("❌ Error: The iTOL file could not be generated. Please check the input data.")
#   }
# }