############################################################
# Mutation Analysis Pipeline (AML Context)
# Author: Hager Salah Abouelnaga
#
# Description:
# This script performs mutation analysis in a translational
# oncology context (AML), including:
# - Variant impact classification
# - Functional category analysis
# - Gene–pathway network construction
#
# This analysis complements RNA-seq findings by linking
# mutation profiles to potential resistance mechanisms
# (e.g., Venetoclax resistance, apoptotic signaling).
#
# Note:
# Data is not included. Replace file paths with your own.
############################################################

# =========================
# Load libraries
# =========================
library(dplyr)
library(ggplot2)
library(tidyverse)

# Additional libraries used in downstream analysis
library(igraph)
library(ggraph)

# =========================
# Config / Paths (EDIT)
# =========================
mutation_file <- "path/to/mutation_data.csv"
output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE)

# =========================
# Read mutation data
# =========================
aml_mutation_data <- read.csv(mutation_file, stringsAsFactors = FALSE)

# Expected columns:
# SYMBOL or Gene
# Impact (HIGH / MODERATE / LOW)
# Category (e.g., Apoptosis, MAPK, PI3K-AKT)

# =========================
# Mutation impact distribution
# =========================
# Assess distribution of mutation severity in AML samples,
# focusing on variants potentially contributing to resistance

impact_summary <- aml_mutation_data %>%
  filter(!is.na(Impact)) %>%
  count(Impact)

ggplot(impact_summary, aes(x = Impact, y = n, fill = Impact)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Variant Impact Distribution (AML)",
       x = "Impact Level",
       y = "Number of Variants")

ggsave(file.path(output_dir, "mutation_impact.png"),
       width = 6, height = 5)

# =========================
# Functional category analysis
# =========================
# Identify biological pathways enriched in mutated genes,
# particularly those related to apoptosis and survival signaling

func_summary <- aml_mutation_data %>%
  filter(!is.na(Category)) %>%
  count(Category, sort = TRUE)

top_func <- func_summary %>% head(15)

ggplot(top_func, aes(x = reorder(Category, n), y = n)) +
  geom_point(size = 4) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top Functional Categories (AML Mutations)",
       x = "Biological Pathway",
       y = "Gene Count")

ggsave(file.path(output_dir, "functional_categories.png"),
       width = 7, height = 5)

# =========================
# Gene–Pathway Network
# =========================
# Construct network linking mutated genes to functional pathways
# This helps visualize pathway-level impact of genomic alterations

edges <- aml_mutation_data %>%
  select(SYMBOL, Category) %>%   # <- keeping SYMBOL (closer to your real script)
  distinct() %>%
  na.omit()

g <- graph_from_data_frame(edges)

ggraph(g, layout = "fr") +
  geom_edge_link(alpha = 0.3) +
  geom_node_point(size = 3) +
  theme_void() +
  ggtitle("Gene–Function Network (AML)")

ggsave(file.path(output_dir, "network_plot.png"),
       width = 7, height = 6)

# =========================
# Optional: Oncoplot
# =========================
# Oncoplot can be generated using maftools if MAF data is available
# This step was used in exploratory analysis but is not shown here
# to maintain clarity and avoid overcrowded visualization

# Example:
# library(maftools)
# maf <- read.maf("path/to/file.maf")
# oncoplot(maf, top = 10)

# =========================
# Final interpretation
# =========================
message("Mutation analysis completed.")

message("High-impact mutations may contribute to altered apoptotic signaling and potential drug resistance in AML.")

message("This analysis integrates with transcriptomic findings to support translational insights.")