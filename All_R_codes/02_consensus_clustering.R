################################################################################
# Script: 02_consensus_clustering.R
# Description: Consensus clustering analysis functions
################################################################################

#' Normalize expression matrix for clustering
#' @param expr_matrix Expression matrix (genes in rows, samples in columns)
#' @return Normalized expression matrix
normalize_for_clustering <- function(expr_matrix) {
  cat("Normalizing expression matrix for clustering\n")
  
  # Median centering
  expr_normalized <- sweep(expr_matrix, 1, apply(expr_matrix, 1, median, na.rm = TRUE))
  
  return(expr_normalized)
}

#' Perform consensus clustering
#' @param expr_matrix Expression matrix (genes in rows, samples in columns)
#' @param max_k Maximum number of clusters
#' @param output_dir Output directory for plots
#' @return List containing clustering results
perform_consensus_clustering <- function(expr_matrix, 
                                        max_k = 9,
                                        output_dir = "results/consensus/") {
  
  cat("\n========== Consensus Clustering Analysis ==========\n")
  cat("Matrix dimensions:", nrow(expr_matrix), "genes,", ncol(expr_matrix), "samples\n")
  cat("Testing k from 2 to", max_k, "\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Run ConsensusClusterPlus
  results <- ConsensusClusterPlus(
    expr_matrix,
    maxK = max_k,
    reps = 1000,
    pItem = 0.8,
    tmyPal = c("white", "#C75D30"),
    pFeature = 1,
    clusterAlg = "km",
    distance = "euclidean",
    seed = 102,
    title = output_dir,
    innerLinkage = 'average',
    plot = 'pdf'
  )
  
  # Calculate ICL
  icl <- calcICL(results, title = output_dir, plot = "pdf")
  
  # Calculate PAC scores
  pac_scores <- calculate_pac_scores(results, max_k)
  
  # Find optimal k
  optimal_k <- which.min(pac_scores$PAC) + 1
  cat("Optimal number of clusters:", optimal_k, "\n")
  
  # Extract cluster assignments for optimal k
  cluster_assignments <- results[[optimal_k]][["consensusClass"]]
  
  # Create cluster data frame
  cluster_df <- data.frame(
    Sample = names(cluster_assignments),
    Cluster = paste0('C', cluster_assignments),
    stringsAsFactors = FALSE
  )
  
  cat("Cluster distribution:\n")
  print(table(cluster_df$Cluster))
  
  cat("====================================================\n\n")
  
  return(list(
    results = results,
    optimal_k = optimal_k,
    pac_scores = pac_scores,
    cluster_assignments = cluster_df,
    icl = icl,
    consensus_matrix = results[[optimal_k]][["consensusMatrix"]]
  ))
}

#' Calculate PAC (Proportion of Ambiguous Clustering) scores
#' @param results ConsensusClusterPlus results
#' @param max_k Maximum number of clusters
#' @return Data frame with PAC scores
calculate_pac_scores <- function(results, max_k) {
  cat("Calculating PAC scores\n")
  
  x1 <- 0.1
  x2 <- 0.9
  Kvec <- 2:max_k
  PAC <- rep(NA, length(Kvec))
  names(PAC) <- paste("K=", Kvec, sep = "")
  
  for (i in Kvec) {
    M <- results[[i]]$consensusMatrix
    Fn <- ecdf(M[lower.tri(M)])
    PAC[i - 1] <- Fn(x2) - Fn(x1)
  }
  
  pac_df <- data.frame(
    K = Kvec,
    PAC = PAC,
    stringsAsFactors = FALSE
  )
  
  return(pac_df)
}

#' Plot PAC scores
#' @param pac_scores Data frame with PAC scores
#' @param save_path Path to save plot
#' @return ggplot object
plot_pac_scores <- function(pac_scores, save_path = NULL) {
  p <- ggplot(pac_scores, aes(factor(K), PAC, group = 1)) +
    geom_line(size = 1.2) +
    geom_point(size = 4, shape = 21, color = 'darkred', fill = 'orange') +
    theme_bw(base_rect_size = 1.5) +
    ggtitle('Proportion of Ambiguous Clustering (PAC)') +
    xlab('Cluster Number K') +
    ylab('PAC Score') +
    theme(
      axis.text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 13)
    ) +
    scale_y_continuous(limits = c(0, max(pac_scores$PAC) * 1.1))
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 8, height = 6, dpi = 300)
    cat("PAC plot saved to:", save_path, "\n")
  }
  
  return(p)
}

#' Plot consensus heatmap
#' @param consensus_matrix Consensus matrix
#' @param cluster_assignments Cluster assignments
#' @param save_path Path to save plot
plot_consensus_heatmap <- function(consensus_matrix, 
                                  cluster_assignments,
                                  save_path = NULL) {
  
  # Order samples by cluster
  cluster_order <- order(cluster_assignments$Cluster)
  ordered_matrix <- consensus_matrix[cluster_order, cluster_order]
  
  # Create annotation
  annotation_df <- data.frame(
    Cluster = cluster_assignments$Cluster[cluster_order],
    row.names = cluster_assignments$Sample[cluster_order]
  )
  
  # Create heatmap
  p <- pheatmap(
    ordered_matrix,
    show_colnames = FALSE,
    show_rownames = FALSE,
    treeheight_row = 20,
    treeheight_col = 20,
    clustering_method = 'complete',
    color = colorRampPalette(c("white", "#C75D30"))(100),
    annotation_row = annotation_df,
    annotation_col = annotation_df,
    annotation_colors = list(
      Cluster = c('C1' = '#4E8279', 'C2' = '#B5739D')
    ),
    main = "Consensus Clustering Heatmap"
  )
  
  if (!is.null(save_path)) {
    pdf(save_path, width = 10, height = 10)
    print(p)
    dev.off()
    cat("Consensus heatmap saved to:", save_path, "\n")
  }
  
  return(p)
}

#' Perform PCA on clustered data
#' @param expr_matrix Expression matrix
#' @param cluster_assignments Cluster assignments
#' @param save_path Path to save plot
#' @return PCA plot
perform_cluster_pca <- function(expr_matrix, 
                               cluster_assignments,
                               save_path = NULL) {
  
  cat("Performing PCA analysis\n")
  
  # Match sample order
  expr_matrix <- expr_matrix[, match(cluster_assignments$Sample, colnames(expr_matrix))]
  
  # Perform PCA
  pca_result <- prcomp(t(expr_matrix), scale = TRUE)
  
  # Extract PC scores
  pc_scores <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Cluster = cluster_assignments$Cluster
  )
  
  # Calculate variance explained
  var_explained <- summary(pca_result)$importance[2, ] * 100
  
  # Create PCA plot
  p <- ggplot(pc_scores, aes(PC1, PC2, color = Cluster, fill = Cluster)) +
    geom_point(size = 3, alpha = 0.7, shape = 21) +
    scale_color_manual(values = c("C1" = "#db6968", "C2" = "#2E9FDF")) +
    scale_fill_manual(values = c("C1" = "#db6968", "C2" = "#2E9FDF")) +
    theme_bw() +
    labs(
      x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 2), "%)"),
      title = "PCA of Consensus Clusters"
    ) +
    theme(
      legend.position = c(0.85, 0.85),
      legend.background = element_rect(fill = 'transparent'),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12)
    )
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 8, height = 6, dpi = 300)
    cat("PCA plot saved to:", save_path, "\n")
  }
  
  return(p)
}
