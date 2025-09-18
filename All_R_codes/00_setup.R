################################################################################
# Script: 00_setup.R
# Description: Environment setup and library loading
################################################################################

# Clean environment
rm(list = ls())
gc()

# Set options
options(stringsAsFactors = FALSE)
options(scipen = 999)

# Install required packages if not already installed
required_packages <- c(
  "tidyverse", "survival", "survminer", "survivalROC", "ConsensusClusterPlus",
  "WGCNA", "limma", "pheatmap", "ggplot2", "ggpubr", "openxlsx", "rms",
  "Mime1", "IOBR", "GSVA", "clusterProfiler", "org.Hs.eg.db", "reshape2",
  "cowplot", "ggsci", "RColorBrewer", "ComplexHeatmap", "Seurat", "SCP",
  "Scissor", "monocle3", "forestmodel", "compareC", "regplot", "ggDCA",
  "rmda", "pec", "introdataviz", "CytoTRACE", "infercnv"
)

# Function to check and install packages
check_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      if (!pkg %in% rownames(installed.packages())) {
        if (pkg %in% c("ConsensusClusterPlus", "WGCNA", "limma", "ComplexHeatmap", 
                       "clusterProfiler", "org.Hs.eg.db", "GSVA", "infercnv")) {
          BiocManager::install(pkg)
        } else {
          install.packages(pkg)
        }
      }
      library(pkg, character.only = TRUE)
    }
  }
}

# Load all required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
  library(survivalROC)
  library(ConsensusClusterPlus)
  library(WGCNA)
  library(limma)
  library(pheatmap)
  library(ggplot2)
  library(ggpubr)
  library(openxlsx)
  library(rms)
  library(reshape2)
  library(cowplot)
  library(ggsci)
  library(RColorBrewer)
})

# Set working directory
setwd("Rworkspace/")

# Create output directories
dir_create <- function(paths) {
  for (path in paths) {
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
      cat("Created directory:", path, "\n")
    }
  }
}

output_dirs <- c(
  "results/",
  "results/consensus/",
  "results/wgcna/",
  "results/ml_models/",
  "results/survival/",
  "results/immune/",
  "results/single_cell/",
  "results/enrichment/",
  "figures/"
)

dir_create(output_dirs)

# Set random seed for reproducibility
set.seed(123)

# Custom color palettes
custom_colors <- list(
  risk = c("High" = "#982b2b", "Low" = "#0074b3"),
  cluster = c("C1" = "#E57164", "C2" = "#A184BC"),
  cell_type = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3'),
  heatmap = colorRampPalette(c("#3288bd", "#66c2a5", "#ffffbf", "#f46d43", "#9e0142"))(100)
)

# Print session info
cat("\n========== Session Information ==========\n")
sessionInfo()
cat("=========================================\n\n")
