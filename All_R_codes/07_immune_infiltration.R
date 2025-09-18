################################################################################
# Script: 07_immune_infiltration.R
# Description: Immune infiltration and microenvironment analysis
################################################################################

library(IOBR)
library(GSVA)
library(estimate)

#' Perform comprehensive immune deconvolution
#' @param expr_matrix Expression matrix (genes in rows, samples in columns)
#' @param methods Deconvolution methods to use
#' @return List of deconvolution results
immune_deconvolution <- function(expr_matrix,
                                 methods = c("cibersort", "epic", "mcpcounter", 
                                           "xcell", "quantiseq", "estimate", "timer")) {
  
  cat("\n========== Immune Deconvolution Analysis ==========\n")
  
  results <- list()
  
  # CIBERSORT
  if ("cibersort" %in% methods) {
    cat("Running CIBERSORT...\n")
    results$cibersort <- deconvo_tme(
      eset = expr_matrix,
      method = "cibersort",
      arrays = TRUE,
      perm = 1000
    )
  }
  
  # EPIC
  if ("epic" %in% methods) {
    cat("Running EPIC...\n")
    results$epic <- deconvo_tme(
      eset = expr_matrix,
      method = "epic",
      arrays = TRUE
    )
  }
  
  # MCP-counter
  if ("mcpcounter" %in% methods) {
    cat("Running MCP-counter...\n")
    results$mcp <- deconvo_tme(
      eset = expr_matrix,
      method = "mcpcounter"
    )
  }
  
  # xCell
  if ("xcell" %in% methods) {
    cat("Running xCell...\n")
    results$xcell <- deconvo_tme(
      eset = expr_matrix,
      method = "xcell",
      arrays = TRUE
    )
  }
  
  # QUANTISEQ
  if ("quantiseq" %in% methods) {
    cat("Running QUANTISEQ...\n")
    results$quantiseq <- deconvo_tme(
      eset = expr_matrix,
      method = "quantiseq",
      tumor = TRUE,
      arrays = TRUE,
      scale_mrna = TRUE
    )
  }
  
  # ESTIMATE
  if ("estimate" %in% methods) {
    cat("Running ESTIMATE...\n")
    results$estimate <- deconvo_tme(
      eset = expr_matrix,
      method = "estimate"
    )
  }
  
  # TIMER
  if ("timer" %in% methods) {
    cat("Running TIMER...\n")
    results$timer <- deconvo_tme(
      eset = expr_matrix,
      method = "timer",
      group_list = rep("gbm", ncol(expr_matrix))
    )
  }
  
  # Combine all results
  cat("Combining deconvolution results...\n")
  tme_combine <- results[[1]]
  for (i in 2:length(results)) {
    tme_combine <- inner_join(tme_combine, results[[i]], by = "ID")
  }
  
  cat("====================================================\n\n")
  
  return(list(
    individual = results,
    combined = tme_combine
  ))
}

#' Compare immune scores between risk groups
#' @param immune_scores Immune infiltration scores
#' @param risk_groups Risk group assignments
#' @param method Statistical test method
#' @return Comparison results
compare_immune_by_risk <- function(immune_scores, 
                                  risk_groups,
                                  method = "wilcox.test") {
  
  cat("Comparing immune scores between risk groups\n")
  
  # Merge data
  merged_data <- merge(
    immune_scores,
    risk_groups,
    by.x = "ID",
    by.y = "row.names"
  )
  
  # Get immune cell columns
  immune_cols <- colnames(immune_scores)[!colnames(immune_scores) %in% c("ID")]
  
  # Perform comparisons
  comparison_results <- data.frame()
  
  for (cell_type in immune_cols) {
    # Skip if all NA
    if (all(is.na(merged_data[[cell_type]]))) {
      next
    }
    
    # Statistical test
    if (method == "wilcox.test") {
      test_result <- wilcox.test(
        as.formula(paste(cell_type, "~ Risk")),
        data = merged_data
      )
    } else if (method == "t.test") {
      test_result <- t.test(
        as.formula(paste(cell_type, "~ Risk")),
        data = merged_data
      )
    }
    
    # Calculate means
    mean_high <- mean(merged_data[merged_data$Risk == "High", cell_type], na.rm = TRUE)
    mean_low <- mean(merged_data[merged_data$Risk == "Low", cell_type], na.rm = TRUE)
    
    # Store results
    result <- data.frame(
      Cell_Type = cell_type,
      Mean_High = mean_high,
      Mean_Low = mean_low,
      Fold_Change = mean_high / mean_low,
      P_Value = test_result$p.value,
      Significant = ifelse(test_result$p.value < 0.001, "***",
                          ifelse(test_result$p.value < 0.01, "**",
                                ifelse(test_result$p.value < 0.05, "*", "ns"))),
      stringsAsFactors = FALSE
    )
    
    comparison_results <- rbind(comparison_results, result)
  }
  
  # Order by p-value
  comparison_results <- comparison_results[order(comparison_results$P_Value), ]
  
  return(comparison_results)
}

#' Create immune infiltration heatmap
#' @param immune_scores Immune infiltration scores
#' @param clinical_data Clinical data with annotations
#' @param scale_data Whether to scale data
#' @param save_path Path to save plot
#' @return Heatmap object
plot_immune_heatmap <- function(immune_scores,
                               clinical_data,
                               scale_data = TRUE,
                               save_path = NULL) {
  
  cat("Creating immune infiltration heatmap\n")
  
  # Prepare data
  rownames(immune_scores) <- immune_scores$ID
  immune_matrix <- as.matrix(immune_scores[, -1])
  
  # Scale if requested
  if (scale_data) {
    immune_matrix <- t(scale(t(immune_matrix)))
    immune_matrix[immune_matrix > 2] <- 2
    immune_matrix[immune_matrix < -2] <- -2
  }
  
  # Match sample order with clinical data
  immune_matrix <- immune_matrix[, match(rownames(clinical_data), colnames(immune_matrix))]
  
  # Create annotation
  annotation_col <- clinical_data[, c("Risk", "Type", "Age"), drop = FALSE]
  
  # Define colors
  annotation_colors <- list(
    Risk = c("High" = "#982b2b", "Low" = "#0074b3"),
    Type = c("GBM" = "#E57164", "LGG" = "#A184BC"),
    Age = c(">65" = "#FF9999", "<=65" = "#9999FF")
  )
  
  # Create heatmap
  p <- pheatmap(
    immune_matrix,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    color = colorRampPalette(c("#3288bd", "white", "#d53e4f"))(100),
    gaps_col = which(diff(as.numeric(factor(clinical_data$Risk))) != 0),
    fontsize_row = 6,
    main = "Immune Infiltration Landscape"
  )
  
  if (!is.null(save_path)) {
    pdf(save_path, width = 12, height = 10)
    print(p)
    dev.off()
    cat("Heatmap saved to:", save_path, "\n")
  }
  
  return(p)
}

#' Calculate immune checkpoint expression
#' @param expr_matrix Expression matrix
#' @param checkpoint_genes Immune checkpoint genes
#' @param risk_groups Risk group assignments
#' @return Checkpoint expression analysis results
analyze_immune_checkpoints <- function(expr_matrix,
                                      checkpoint_genes = c("PDCD1", "CD274", "PDCD1LG2", 
                                                          "CTLA4", "HAVCR2", "LAG3", 
                                                          "TIGIT", "IDO1", "CD276"),
                                      risk_groups) {
  
  cat("Analyzing immune checkpoint expression\n")
  
  # Filter for available genes
  available_checkpoints <- checkpoint_genes[checkpoint_genes %in% rownames(expr_matrix)]
  cat("Available checkpoint genes:", length(available_checkpoints), "/", length(checkpoint_genes), "\n")
  
  # Extract expression
  checkpoint_expr <- expr_matrix[available_checkpoints, , drop = FALSE]
  
  # Merge with risk groups
  checkpoint_data <- data.frame(
    t(checkpoint_expr),
    Risk = risk_groups[match(colnames(checkpoint_expr), rownames(risk_groups)), "Risk"]
  )
  
  # Compare expression between groups
  comparison_results <- data.frame()
  
  for (gene in available_checkpoints) {
    # Wilcoxon test
    test_result <- wilcox.test(
      as.formula(paste(gene, "~ Risk")),
      data = checkpoint_data
    )
    
    # Calculate medians
    median_high <- median(checkpoint_data[checkpoint_data$Risk == "High", gene], na.rm = TRUE)
    median_low <- median(checkpoint_data[checkpoint_data$Risk == "Low", gene], na.rm = TRUE)
    
    result <- data.frame(
      Gene = gene,
      Median_High = median_high,
      Median_Low = median_low,
      Fold_Change = median_high / median_low,
      P_Value = test_result$p.value,
      stringsAsFactors = FALSE
    )
    
    comparison_results <- rbind(comparison_results, result)
  }
  
  return(list(
    expression_data = checkpoint_data,
    comparison = comparison_results
  ))
}

#' Calculate immunophenoscore (IPS)
#' @param expr_matrix Expression matrix
#' @param plot Whether to generate plots
#' @return IPS scores
calculate_ips <- function(expr_matrix, plot = FALSE) {
  
  cat("Calculating immunophenoscore (IPS)\n")
  
  ips_results <- deconvo_tme(
    eset = expr_matrix,
    method = "ips",
    plot = plot
  )
  
  return(ips_results)
}

#' Perform ssGSEA for immune signatures
#' @param expr_matrix Expression matrix
#' @param immune_signatures List of immune gene signatures
#' @param method ssGSEA method
#' @return ssGSEA scores
immune_ssgsea <- function(expr_matrix,
                         immune_signatures = NULL,
                         method = "ssgsea") {
  
  cat("Performing ssGSEA for immune signatures\n")
  
  # Use default signatures if not provided
  if (is.null(immune_signatures)) {
    # Load default immune signatures
    immune_signatures <- list(
      CD8_T_cells = c("CD8A", "CD8B", "CD3D", "CD3E", "CD3G"),
      CD4_T_cells = c("CD4", "CD3D", "CD3E", "CD3G"),
      B_cells = c("CD19", "CD79A", "MS4A1", "CD20"),
      NK_cells = c("KLRB1", "KLRD1", "KLRK1", "KLRC1"),
      Macrophages = c("CD68", "CD163", "MARCO"),
      DC = c("ITGAX", "CD1C", "BATF3"),
      Tregs = c("FOXP3", "IL2RA", "IKZF2"),
      MDSC = c("CD14", "CD33", "IL4R", "ITGAM")
    )
  }
  
  # Run ssGSEA
  ssgsea_scores <- gsva(
    expr_matrix,
    immune_signatures,
    method = method,
    kcdf = "Gaussian",
    mx.diff = FALSE,
    verbose = TRUE
  )
  
  # Normalize scores (0-1)
  ssgsea_normalized <- ssgsea_scores
  for (i in 1:nrow(ssgsea_normalized)) {
    ssgsea_normalized[i, ] <- (ssgsea_normalized[i, ] - min(ssgsea_normalized[i, ])) / 
                              (max(ssgsea_normalized[i, ]) - min(ssgsea_normalized[i, ]))
  }
  
  return(list(
    raw_scores = ssgsea_scores,
    normalized_scores = ssgsea_normalized
  ))
}

#' Visualize immune scores by risk group
#' @param immune_scores Immune infiltration scores
#' @param risk_groups Risk group assignments
#' @param score_types Types of scores to plot
#' @param save_path Path to save plot
#' @return Boxplot object
plot_immune_boxplot <- function(immune_scores,
                               risk_groups,
                               score_types = c("StromalScore", "ImmuneScore", 
                                             "ESTIMATEScore", "TumorPurity"),
                               save_path = NULL) {
  
  # Prepare data
  plot_data <- merge(
    immune_scores[, c("ID", score_types)],
    risk_groups,
    by.x = "ID",
    by.y = "row.names"
  )
  
  # Reshape for plotting
  plot_data_long <- reshape2::melt(
    plot_data,
    id.vars = c("ID", "Risk"),
    variable.name = "Score_Type",
    value.name = "Score"
  )
  
  # Create plot
  p <- ggplot(plot_data_long, aes(x = Risk, y = Score, fill = Risk)) +
    geom_boxplot(alpha = 0.7) +
    geom_point(position = position_jitter(width = 0.2), alpha = 0.3) +
    facet_wrap(~ Score_Type, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c("High" = "#E57164", "Low" = "#A184BC")) +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 12, face = "bold")
    ) +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",
      label.y.npc = 0.95
    ) +
    labs(
      x = "",
      y = "Score",
      title = "Immune Scores by Risk Group"
    )
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 10, height = 8, dpi = 300)
    cat("Plot saved to:", save_path, "\n")
  }
  
  return(p)
}
