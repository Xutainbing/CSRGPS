################################################################################
# Script: 08_single_cell_analysis.R
# Description: Single-cell RNA-seq analysis including Scissor
################################################################################

library(Seurat)
library(SCP)
library(Scissor)
library(monocle3)

#' Load and preprocess single-cell data
#' @param sc_path Path to single-cell data
#' @param min_cells Minimum cells per gene
#' @param min_features Minimum features per cell
#' @return Seurat object
load_single_cell_data <- function(sc_path,
                                 min_cells = 3,
                                 min_features = 200) {
  
  cat("Loading single-cell data from:", sc_path, "\n")
  
  # Load data
  sc_data <- readRDS(sc_path)
  
  # Basic QC if not already done
  if (!"percent.mt" %in% colnames(sc_data@meta.data)) {
    sc_data[["percent.mt"]] <- PercentageFeatureSet(sc_data, pattern = "^MT-")
  }
  
  cat("Single-cell data:", ncol(sc_data), "cells,", nrow(sc_data), "genes\n")
  cat("Samples:", length(unique(sc_data$orig.ident)), "\n")
  
  return(sc_data)
}

#' Integrate single-cell data using SCP
#' @param sc_data Seurat object
#' @param batch_column Batch column name
#' @param integration_method Integration method
#' @param clustering_resolutions Clustering resolutions to test
#' @return Integrated Seurat object
integrate_single_cell <- function(sc_data,
                                 batch_column = "orig.ident",
                                 integration_method = "fastMNN",
                                 clustering_resolutions = seq(0.1, 1, 0.1)) {
  
  cat("\n========== Single-Cell Integration ==========\n")
  cat("Integration method:", integration_method, "\n")
  cat("Batch column:", batch_column, "\n")
  
  sc_integrated <- Integration_SCP(
    srtMerge = sc_data,
    batch = batch_column,
    integration_method = integration_method,
    linear_reduction_dims = 100,
    nHVF = 2000,
    linear_reduction_dims_use = 1:50,
    cluster_resolution = clustering_resolutions,
    nonlinear_reduction = c("umap", "tsne")
  )
  
  cat("Integration complete\n")
  cat("============================================\n\n")
  
  return(sc_integrated)
}

#' Perform Scissor analysis
#' @param bulk_expr Bulk expression matrix
#' @param sc_data Seurat object
#' @param phenotype Binary phenotype vector
#' @param alpha_range Alpha values to test
#' @param cutoff Cutoff for cell selection
#' @param family Regression family
#' @return Scissor results
perform_scissor <- function(bulk_expr,
                          sc_data,
                          phenotype,
                          alpha_range = seq(0.001, 0.01, 0.002),
                          cutoff = 0.2,
                          family = "binomial") {
  
  cat("\n========== Scissor Analysis ==========\n")
  cat("Bulk samples:", ncol(bulk_expr), "\n")
  cat("Single cells:", ncol(sc_data), "\n")
  cat("Phenotype distribution:", table(phenotype), "\n")
  
  # Convert phenotype to binary if needed
  if (is.character(phenotype) || is.factor(phenotype)) {
    phenotype_binary <- as.numeric(factor(phenotype)) - 1
    tag <- levels(factor(phenotype))
  } else {
    phenotype_binary <- phenotype
    tag <- c("Low", "High")
  }
  
  # Run Scissor
  scissor_results <- Scissor(
    bulk_matrix = bulk_expr,
    sc_dataset = sc_data,
    phenotype = phenotype_binary,
    tag = tag,
    alpha = alpha_range,
    cutoff = cutoff,
    family = family
  )
  
  # Add Scissor results to metadata
  scissor_select <- rep("Background", ncol(sc_data))
  names(scissor_select) <- colnames(sc_data)
  scissor_select[scissor_results$Scissor_pos] <- paste0(tag[2], "_associated")
  scissor_select[scissor_results$Scissor_neg] <- paste0(tag[1], "_associated")
  
  sc_data <- AddMetaData(
    sc_data,
    metadata = scissor_select,
    col.name = "scissor"
  )
  
  cat("\nScissor+ cells:", length(scissor_results$Scissor_pos), "\n")
  cat("Scissor- cells:", length(scissor_results$Scissor_neg), "\n")
  cat("========================================\n\n")
  
  return(list(
    sc_data = sc_data,
    scissor_results = scissor_results,
    selected_cells = list(
      positive = scissor_results$Scissor_pos,
      negative = scissor_results$Scissor_neg
    )
  ))
}

#' Visualize Scissor results
#' @param sc_data Seurat object with Scissor results
#' @param reduction Dimension reduction to use
#' @param save_path Path to save plot
#' @return Plot object
plot_scissor_results <- function(sc_data,
                                reduction = "umap",
                                save_path = NULL) {
  
  # UMAP colored by Scissor results
  p1 <- DimPlot(
    sc_data,
    reduction = reduction,
    group.by = "scissor",
    cols = c("grey90", "indianred1", "royalblue"),
    pt.size = 0.5,
    order = c("High_associated", "Low_associated", "Background")
  ) +
    theme(
      legend.position = c(0.8, 0.2),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    ggtitle("Scissor Cell Selection")
  
  # Cell type composition
  if ("celltype" %in% colnames(sc_data@meta.data)) {
    composition_data <- table(sc_data$celltype, sc_data$scissor)
    composition_df <- as.data.frame(composition_data)
    colnames(composition_df) <- c("CellType", "Scissor", "Count")
    
    p2 <- ggplot(composition_df, aes(x = CellType, y = Count, fill = Scissor)) +
      geom_bar(stat = "identity", position = "fill") +
      scale_fill_manual(values = c("grey90", "indianred1", "royalblue")) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
      ) +
      labs(
        x = "",
        y = "Proportion",
        title = "Scissor Cells by Cell Type"
      )
    
    # Combine plots
    p_combined <- cowplot::plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 1.2))
  } else {
    p_combined <- p1
  }
  
  if (!is.null(save_path)) {
    ggsave(save_path, p_combined, width = 14, height = 6, dpi = 300)
    cat("Plot saved to:", save_path, "\n")
  }
  
  return(p_combined)
}

#' Perform trajectory analysis using Monocle3
#' @param sc_data Seurat object
#' @param use_reduction Whether to use existing reduction
#' @return Monocle3 CDS object
trajectory_analysis <- function(sc_data,
                               use_reduction = TRUE) {
  
  cat("\n========== Trajectory Analysis ==========\n")
  
  # Convert to Monocle3 CDS
  data <- GetAssayData(sc_data, assay = 'RNA', slot = 'counts')
  cell_metadata <- sc_data@meta.data
  gene_annotation <- data.frame(
    gene_short_name = rownames(data),
    row.names = rownames(data)
  )
  
  cds <- new_cell_data_set(
    data,
    cell_metadata = cell_metadata,
    gene_metadata = gene_annotation
  )
  
  # Preprocess
  cat("Preprocessing CDS...\n")
  cds <- preprocess_cds(cds, num_dim = 50)
  
  # Use existing UMAP if available
  if (use_reduction && "umap" %in% names(sc_data@reductions)) {
    cat("Using existing UMAP coordinates\n")
    cds_embed <- cds@int_colData$reducedDims$UMAP
    int_embed <- Embeddings(sc_data, reduction = "umap")
    int_embed <- int_embed[rownames(cds_embed), ]
    cds@int_colData$reducedDims$UMAP <- int_embed
  } else {
    # Align if batch effect exists
    if ("orig.ident" %in% colnames(cell_metadata)) {
      cds <- align_cds(cds, alignment_group = "orig.ident")
    }
    
    # Reduce dimensions
    cat("Reducing dimensions...\n")
    cds <- reduce_dimension(cds, reduction_method = "UMAP", cores = 8)
  }
  
  # Cluster cells
  cat("Clustering cells...\n")
  cds <- cluster_cells(cds)
  
  # Learn trajectory
  cat("Learning trajectory...\n")
  cds <- learn_graph(cds)
  
  # Order cells in pseudotime
  cat("Ordering cells in pseudotime...\n")
  cds <- order_cells(cds)
  
  cat("==========================================\n\n")
  
  return(cds)
}

#' Find trajectory-associated genes
#' @param cds Monocle3 CDS object
#' @param q_cutoff Q-value cutoff
#' @param top_n Number of top genes to return
#' @return Trajectory genes
find_trajectory_genes <- function(cds,
                                 q_cutoff = 0.01,
                                 top_n = 100) {
  
  cat("Finding trajectory-associated genes...\n")
  
  # Graph test
  trajectory_genes <- graph_test(
    cds,
    neighbor_graph = "principal_graph",
    cores = 10
  )
  
  # Filter significant genes
  trajectory_genes_sig <- trajectory_genes %>%
    filter(q_value < q_cutoff) %>%
    arrange(desc(morans_I)) %>%
    head(top_n)
  
  cat("Found", nrow(trajectory_genes_sig), "significant genes\n")
  
  return(trajectory_genes_sig)
}

#' Perform CytoTRACE analysis
#' @param sc_data Seurat object
#' @param ncores Number of cores to use
#' @return CytoTRACE results
perform_cytotrace <- function(sc_data, ncores = 10) {
  
  cat("\n========== CytoTRACE Analysis ==========\n")
  
  library(CytoTRACE)
  
  # Extract expression matrix
  mat <- as.matrix(sc_data@assays$RNA@counts)
  
  # Get phenotype if available
  if ("celltype" %in% colnames(sc_data@meta.data)) {
    phenotype <- as.character(sc_data$celltype)
    names(phenotype) <- rownames(sc_data@meta.data)
  } else {
    phenotype <- NULL
  }
  
  # Run CytoTRACE
  cat("Running CytoTRACE...\n")
  results <- CytoTRACE(
    mat = mat,
    ncores = ncores
  )
  
  # Add to metadata
  sc_data$cytotrace_score <- results$CytoTRACE[colnames(sc_data)]
  
  cat("=========================================\n\n")
  
  return(list(
    sc_data = sc_data,
    cytotrace_results = results
  ))
}

#' Perform InferCNV analysis
#' @param sc_data Seurat object
#' @param reference_cells Reference cell barcodes
#' @param gene_order_file Gene position file
#' @param output_dir Output directory
#' @return InferCNV object
perform_infercnv <- function(sc_data,
                            reference_cells,
                            gene_order_file,
                            output_dir = "results/infercnv/") {
  
  cat("\n========== InferCNV Analysis ==========\n")
  
  library(infercnv)
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Extract counts
  counts_matrix <- GetAssayData(sc_data, slot = "counts")
  
  # Create annotation file
  cell_annotations <- data.frame(
    cell = colnames(counts_matrix),
    celltype = sc_data$celltype,
    stringsAsFactors = FALSE
  )
  
  # Identify reference groups
  ref_groups <- unique(cell_annotations$celltype[cell_annotations$cell %in% reference_cells])
  
  # Write files
  write.table(
    counts_matrix,
    file.path(output_dir, "counts_matrix.txt"),
    sep = "\t",
    quote = FALSE
  )
  
  write.table(
    cell_annotations,
    file.path(output_dir, "annotations.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  # Create InferCNV object
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = file.path(output_dir, "counts_matrix.txt"),
    annotations_file = file.path(output_dir, "annotations.txt"),
    delim = "\t",
    gene_order_file = gene_order_file,
    ref_group_names = ref_groups
  )
  
  # Run InferCNV
  cat("Running InferCNV...\n")
  infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff = 0.1,
    out_dir = output_dir,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = TRUE,
    output_format = "pdf",
    num_threads = 8,
    write_expr_matrix = TRUE
  )
  
  cat("=========================================\n\n")
  
  return(infercnv_obj)
}

#' Calculate CNV scores from InferCNV results
#' @param infercnv_results Path to InferCNV observation matrix
#' @return CNV scores
calculate_cnv_scores <- function(infercnv_results) {
  
  cat("Calculating CNV scores...\n")
  
  # Read observation matrix
  obs_matrix <- fread(infercnv_results)
  obs_matrix <- as.data.frame(obs_matrix)
  rownames(obs_matrix) <- obs_matrix$V1
  obs_matrix <- obs_matrix[, -1]
  
  # Scale and calculate scores
  obs_scaled <- scale(t(obs_matrix))
  obs_scaled <- scales::rescale(obs_scaled, to = c(-1, 1))
  
  # Calculate CNV score
  cnv_scores <- colSums(obs_scaled * obs_scaled)
  
  return(data.frame(
    Cell = names(cnv_scores),
    CNV_Score = cnv_scores,
    stringsAsFactors = FALSE
  ))
}
