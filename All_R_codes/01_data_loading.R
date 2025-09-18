################################################################################
# Script: 01_data_loading.R
# Description: Functions for loading various datasets
################################################################################

#' Load TCGA GBM/LGG data
#' @param data_path Path to RData file
#' @return List containing expression matrix and clinical data
load_tcga_data <- function(data_path = "data/gbm_lgg_tpm_cl.Rdata") {
  cat("Loading TCGA data from:", data_path, "\n")
  load(data_path)
  
  # Check data integrity
  if (!identical(colnames(tpm), rownames(cl))) {
    stop("Sample IDs do not match between expression and clinical data")
  }
  
  # Process clinical data
  cl$ID <- rownames(cl)
  cl$OS.time <- as.numeric(cl$OS.time)
  cl$OS <- as.numeric(cl$OS)
  
  # Process age and gender
  cl$Age <- ifelse(cl$Age > 65, ">65", "<=65")
  cl$Gender <- ifelse(cl$Gender == "female", "Female", "Male")
  
  cat("Loaded TCGA data: ", nrow(tpm), " genes, ", ncol(tpm), " samples\n")
  
  return(list(
    expression = tpm,
    clinical = cl,
    sample_ids = colnames(tpm)
  ))
}

#' Load CGGA array data (301 samples)
#' @param data_path Path to RData file
#' @return Data frame with combined expression and clinical data
load_cgga_array <- function(data_path = "data/CGGA_mrna_array_301.Rdata") {
  cat("Loading CGGA array data from:", data_path, "\n")
  load(data_path)
  
  # Process clinical data
  surv2 <- surv2[, c("CGGA_ID", "OS", "Censor..alive.0..dead.1.")]
  colnames(surv2) <- c("ID", "OS.time", "OS")
  rownames(surv2) <- surv2$ID
  
  # Check sample matching
  if (!identical(rownames(surv2), colnames(exp2))) {
    exp2 <- exp2[, match(rownames(surv2), colnames(exp2))]
  }
  
  # Combine data
  cgga_array_dat <- cbind(surv2, t(exp2))
  cgga_array_dat <- na.omit(cgga_array_dat)
  
  cat("Loaded CGGA array data: ", nrow(cgga_array_dat), " samples\n")
  
  return(cgga_array_dat)
}

#' Load CGGA RNA-seq 325 data
#' @param data_path Path to RData file
#' @return Data frame with combined expression and clinical data
load_cgga_rna325 <- function(data_path = "data/CGGA_mrna_rna325.Rdata") {
  cat("Loading CGGA RNA-seq 325 data from:", data_path, "\n")
  load(data_path)
  
  # Process clinical data
  surv3 <- surv3[, c("CGGA_ID", "OS", "Censor..alive.0..dead.1.")]
  colnames(surv3) <- c("ID", "OS.time", "OS")
  rownames(surv3) <- surv3$ID
  
  # Check sample matching
  if (!identical(rownames(surv3), colnames(exp3))) {
    exp3 <- exp3[, match(rownames(surv3), colnames(exp3))]
  }
  
  # Combine data
  cgga_rna325_dat <- cbind(surv3, t(exp3))
  cgga_rna325_dat <- na.omit(cgga_rna325_dat)
  
  cat("Loaded CGGA RNA-seq 325 data: ", nrow(cgga_rna325_dat), " samples\n")
  
  return(cgga_rna325_dat)
}

#' Load CGGA RNA-seq 693 data
#' @param expr_path Path to expression file
#' @param clin_path Path to clinical file
#' @return Data frame with combined expression and clinical data
load_cgga_rna693 <- function(expr_path = "data/CGGA.mRNAseq_693.RSEM-genes.20200506.txt",
                              clin_path = "data/CGGA.mRNAseq_693_clinical.20200506.txt") {
  
  cat("Loading CGGA RNA-seq 693 data\n")
  
  # Load expression data
  cgga_rna693_exp <- read.table(expr_path, header = TRUE, row.names = 1)
  cgga_rna693_exp <- log2(cgga_rna693_exp + 1)
  
  # Load clinical data
  cgga_rna693_cl <- read.table(clin_path, sep = "\t", header = TRUE, check.names = FALSE)
  cgga_rna693_cl <- cgga_rna693_cl[, c("CGGA_ID", "OS", "Censor (alive=0; dead=1)")]
  colnames(cgga_rna693_cl) <- c("ID", "OS.time", "OS")
  rownames(cgga_rna693_cl) <- cgga_rna693_cl$ID
  
  # Filter and match samples
  cgga_rna693_cl <- cgga_rna693_cl[!is.na(cgga_rna693_cl$OS), ]
  cgga_rna693_cl <- cgga_rna693_cl[!is.na(cgga_rna693_cl$OS.time), ]
  cgga_rna693_exp <- cgga_rna693_exp[, match(rownames(cgga_rna693_cl), colnames(cgga_rna693_exp))]
  
  # Combine data
  cgga_rna693_dat <- cbind(cgga_rna693_cl, t(cgga_rna693_exp))
  cgga_rna693_dat <- na.omit(cgga_rna693_dat)
  
  cat("Loaded CGGA RNA-seq 693 data: ", nrow(cgga_rna693_dat), " samples\n")
  
  return(cgga_rna693_dat)
}

#' Load GSE16011 data
#' @param expr_path Path to expression file
#' @param clin_path Path to clinical file
#' @return Data frame with combined expression and clinical data
load_gse16011 <- function(expr_path = "data/GSE16011_exp.txt",
                          clin_path = "data/GSE16011_cl.txt") {
  
  cat("Loading GSE16011 data\n")
  
  # Load expression data
  gse16011_exp <- read.table(expr_path, header = TRUE, row.names = 1)
  gse16011_exp <- gse16011_exp[, -1]
  gse16011_exp <- log2(gse16011_exp + 1)
  
  # Load clinical data
  gse16011_cl <- fread(clin_path)
  gse16011_cl <- as.data.frame(gse16011_cl)
  gse16011_cl <- gse16011_cl[, -2]
  rownames(gse16011_cl) <- gse16011_cl$`H:symbol`
  gse16011_cl <- gse16011_cl[, -1]
  gse16011_cl <- t(gse16011_cl)
  gse16011_cl <- as.data.frame(gse16011_cl)
  
  # Process clinical data
  gse16011_cl$ID <- rownames(gse16011_cl)
  gse16011_cl <- gse16011_cl[, c("ID", "overall survival", "vital status")]
  colnames(gse16011_cl) <- c("ID", "OS.time", "OS")
  
  # Filter valid samples
  gse16011_cl <- gse16011_cl[gse16011_cl$OS != "nd", ]
  gse16011_cl <- gse16011_cl[gse16011_cl$OS != "lost_to_follow_up", ]
  gse16011_cl$OS <- ifelse(gse16011_cl$OS == 'dead', 1, 0)
  gse16011_cl <- gse16011_cl[gse16011_cl$OS.time != "nd", ]
  gse16011_cl$OS.time <- as.numeric(gse16011_cl$OS.time) * 365
  
  # Match samples
  gse16011_exp <- gse16011_exp[, match(gse16011_cl$ID, colnames(gse16011_exp))]
  
  # Combine data
  gse16011_dat <- cbind(gse16011_cl, t(gse16011_exp))
  
  cat("Loaded GSE16011 data: ", nrow(gse16011_dat), " samples\n")
  
  return(gse16011_dat)
}

#' Load senescence-related genes
#' @param file_path Path to Excel file containing gene list
#' @return Character vector of gene symbols
load_srg_genes <- function(file_path = "data/SRG_CellAge.xlsx") {
  cat("Loading senescence-related genes from:", file_path, "\n")
  
  cell_age <- read.xlsx(file_path)
  srg_genes <- unique(cell_age$Gene.symbol)
  
  cat("Loaded", length(srg_genes), "senescence-related genes\n")
  
  return(srg_genes)
}

#' Load all datasets for meta-analysis
#' @return List of all datasets
load_all_datasets <- function() {
  cat("\n========== Loading All Datasets ==========\n")
  
  datasets <- list()
  
  # TCGA
  tcga_data <- load_tcga_data()
  datasets$TCGA <- cbind(
    tcga_data$clinical[, c("ID", "OS.time", "OS")],
    t(tcga_data$expression)
  )
  
  # CGGA datasets
  datasets$CGGA_array <- load_cgga_array()
  datasets$CGGA_rna325 <- load_cgga_rna325()
  datasets$CGGA_rna693 <- load_cgga_rna693()
  
  # GEO datasets
  datasets$GSE16011 <- load_gse16011()
  # Add other GEO datasets as needed
  # datasets$Rembrandt <- load_rembrandt()
  # datasets$GSE4412 <- load_gse4412()
  # datasets$E_MTAB_3892 <- load_emtab3892()
  
  cat("\nLoaded", length(datasets), "datasets\n")
  cat("==========================================\n\n")
  
  return(datasets)
}
