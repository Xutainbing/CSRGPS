################################################################################
# Script: 09_immunotherapy_analysis.R
# Description: Immunotherapy response prediction and analysis
################################################################################

#' Load immunotherapy dataset
#' @param dataset_name Name of the dataset
#' @param expr_path Path to expression file
#' @param clin_path Path to clinical file
#' @return Processed dataset
load_immunotherapy_dataset <- function(dataset_name,
                                      expr_path,
                                      clin_path) {
  
  cat("Loading", dataset_name, "dataset\n")
  
  # Load expression data
  if (grepl("\\.rds$", expr_path, ignore.case = TRUE)) {
    expr_data <- readRDS(expr_path)
    if ("GENE_SYMBOL" %in% colnames(expr_data)) {
      rownames(expr_data) <- expr_data$GENE_SYMBOL
      expr_data <- expr_data[, -1]
    }
  } else {
    expr_data <- read.table(expr_path, header = TRUE, row.names = 1)
  }
  
  # Log transform if needed
  if (max(expr_data, na.rm = TRUE) > 100) {
    expr_data <- log2(expr_data + 1)
  }
  
  # Load clinical data
  if (grepl("\\.tsv$", clin_path, ignore.case = TRUE)) {
    clin_data <- read.table(clin_path, sep = "\t", header = TRUE)
  } else {
    clin_data <- read.csv(clin_path)
  }
  
  # Process clinical data
  if ("sample_id" %in% colnames(clin_data)) {
    rownames(clin_data) <- clin_data$sample_id
  }
  
  # Extract key columns
  clin_subset <- clin_data[, c("OS", "OS.time")]
  
  # Add response if available
  if ("response_NR" %in% colnames(clin_data)) {
    clin_subset$Response <- clin_data$response_NR
  } else if ("Best.Confirmed.Overall.Response" %in% colnames(clin_data)) {
    clin_subset$Response <- ifelse(
      clin_data$Best.Confirmed.Overall.Response %in% c("CR", "PR"),
      "R", "N"
    )
  }
  
  # Match samples
  common_samples <- intersect(rownames(clin_subset), colnames(expr_data))
  expr_data <- expr_data[, common_samples]
  clin_subset <- clin_subset[common_samples, ]
  
  # Combine data
  combined_data <- cbind(
    ID = rownames(clin_subset),
    clin_subset,
    t(expr_data)
  )
  
  cat("Loaded", nrow(combined_data), "samples,", ncol(expr_data), "genes\n")
  
  return(combined_data)
}

#' Predict immunotherapy response using risk score
#' @param risk_scores Risk scores from the model
#' @param response_data Response data
#' @return Response prediction results
predict_immunotherapy_response <- function(risk_scores,
                                         response_data) {
  
  cat("\n========== Immunotherapy Response Prediction ==========\n")
  
  # Merge risk scores with response
  merged_data <- merge(
    risk_scores,
    response_data,
    by.x = "row.names",
    by.y = "row.names"
  )
  
  # Survival analysis by risk group
  surv_fit <- survfit(Surv(OS.time, OS) ~ Risk, data = merged_data)
  surv_diff <- survdiff(Surv(OS.time, OS) ~ Risk, data = merged_data)
  p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  
  cat("Survival difference p-value:", format.pval(p_value), "\n")
  
  # Response rate comparison if available
  if ("Response" %in% colnames(merged_data)) {
    response_table <- table(merged_data$Risk, merged_data$Response)
    response_rate <- prop.table(response_table, 1)
    
    cat("\nResponse rates:\n")
    print(response_rate)
    
    # Chi-square test
    chi_test <- chisq.test(response_table)
    cat("Chi-square test p-value:", format.pval(chi_test$p.value), "\n")
  }
  
  cat("========================================================\n\n")
  
  return(list(
    merged_data = merged_data,
    survival_fit = surv_fit,
    p_value = p_value,
    response_analysis = if ("Response" %in% colnames(merged_data)) {
      list(
        table = response_table,
        rates = response_rate,
        chi_test = chi_test
      )
    } else NULL
  ))
}

#' Calculate TIDE scores
#' @param expr_matrix Expression matrix
#' @param output_file Output file path
#' @return TIDE prediction results
calculate_tide <- function(expr_matrix,
                          output_file = "results/tide_input.txt") {
  
  cat("Preparing data for TIDE analysis\n")
  
  # Normalize expression data
  expr_normalized <- t(apply(expr_matrix, 1, function(x) {x - mean(x)}))
  
  # Write normalized data
  write.table(
    expr_normalized,
    output_file,
    sep = "\t",
    quote = FALSE,
    col.names = TRUE,
    row.names = TRUE
  )
  
  cat("TIDE input file saved to:", output_file, "\n")
  cat("Please submit to: http://tide.dfci.harvard.edu/\n")
  
  return(output_file)
}

#' Analyze TIDE results
#' @param tide_results Path to TIDE result file
#' @param risk_groups Risk group assignments
#' @return TIDE analysis results
analyze_tide_results <- function(tide_results,
                                risk_groups) {
  
  cat("Analyzing TIDE results\n")
  
  # Read TIDE results
  tide_data <- read.csv(tide_results, row.names = 1)
  
  # Add risk groups
  tide_data$Risk <- risk_groups[match(rownames(tide_data), rownames(risk_groups)), "Risk"]
  
  # Compare TIDE scores
  tide_comparison <- list()
  
  # TIDE score
  tide_comparison$TIDE <- wilcox.test(TIDE ~ Risk, data = tide_data)
  
  # Dysfunction score
  tide_comparison$Dysfunction <- wilcox.test(Dysfunction ~ Risk, data = tide_data)
  
  # Exclusion score
  tide_comparison$Exclusion <- wilcox.test(Exclusion ~ Risk, data = tide_data)
  
  # MSI score
  tide_comparison$MSI <- wilcox.test(MSI.Score ~ Risk, data = tide_data)
  
  # Response prediction
  response_table <- table(tide_data$Risk, tide_data$Responder)
  tide_comparison$Response_chi <- chisq.test(response_table)
  
  # Print results
  cat("\nTIDE Score comparison:\n")
  cat("P-value:", format.pval(tide_comparison$TIDE$p.value), "\n")
  
  cat("\nResponse prediction:\n")
  print(prop.table(response_table, 1))
  cat("Chi-square p-value:", format.pval(tide_comparison$Response_chi$p.value), "\n")
  
  return(list(
    tide_data = tide_data,
    comparisons = tide_comparison,
    response_table = response_table
  ))
}

#' Perform SubMap analysis
#' @param expr_data Expression matrix for query dataset
#' @param query_labels Query dataset labels
#' @param reference_data Expression matrix for reference dataset
#' @param reference_labels Reference dataset labels
#' @param output_dir Output directory
#' @return SubMap results
perform_submap <- function(expr_data,
                          query_labels,
                          reference_data,
                          reference_labels,
                          output_dir = "results/submap/") {
  
  cat("Preparing data for SubMap analysis\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Write query dataset
  query_file <- file.path(output_dir, "query_expr.txt")
  write.table(expr_data, query_file, sep = "\t", quote = FALSE)
  
  # Write query labels
  query_label_file <- file.path(output_dir, "query_labels.txt")
  write.table(
    t(as.numeric(factor(query_labels)) - 1),
    query_label_file,
    sep = "\t",
    quote = FALSE,
    col.names = colnames(expr_data),
    row.names = FALSE
  )
  
  # Write reference dataset
  ref_file <- file.path(output_dir, "reference_expr.txt")
  write.table(reference_data, ref_file, sep = "\t", quote = FALSE)
  
  # Write reference labels
  ref_label_file <- file.path(output_dir, "reference_labels.txt")
  write.table(
    reference_labels,
    ref_label_file,
    sep = "\t",
    quote = FALSE,
    row.names = TRUE
  )
  
  cat("SubMap input files saved to:", output_dir, "\n")
  cat("Please run SubMap analysis using GenePattern\n")
  
  return(list(
    query_file = query_file,
    query_label_file = query_label_file,
    reference_file = ref_file,
    reference_label_file = ref_label_file
  ))
}

#' Visualize immunotherapy response
#' @param response_data Response prediction data
#' @param plot_type Type of plot
#' @param save_path Path to save plot
#' @return Plot object
plot_immunotherapy_response <- function(response_data,
                                       plot_type = c("survival", "response_rate", "both"),
                                       save_path = NULL) {
  
  plot_type <- match.arg(plot_type)
  
  plots <- list()
  
  # Survival plot
  if (plot_type %in% c("survival", "both")) {
    surv_fit <- survfit(Surv(OS.time, OS) ~ Risk, data = response_data)
    
    plots$survival <- ggsurvplot(
      surv_fit,
      data = response_data,
      palette = c("#982b2b", "#0074b3"),
      pval = TRUE,
      conf.int = TRUE,
      risk.table = TRUE,
      risk.table.height = 0.25,
      legend.labs = c("High Risk", "Low Risk"),
      title = "Overall Survival - Immunotherapy Cohort"
    )
  }
  
  # Response rate plot
  if (plot_type %in% c("response_rate", "both") && "Response" %in% colnames(response_data)) {
    response_summary <- response_data %>%
      group_by(Risk, Response) %>%
      summarise(Count = n()) %>%
      group_by(Risk) %>%
      mutate(Proportion = Count / sum(Count))
    
    plots$response <- ggplot(response_summary, aes(x = Risk, y = Proportion, fill = Response)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("N" = "#A184BC", "R" = "#E57164")) +
      geom_text(aes(label = scales::percent(Proportion)), 
               position = position_stack(vjust = 0.5),
               color = "white", size = 5) +
      theme_bw() +
      labs(
        x = "",
        y = "Proportion",
        title = "Response Rate by Risk Group"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.position = "right"
      )
  }
  
  # Combine plots if both
  if (plot_type == "both" && length(plots) == 2) {
    combined_plot <- cowplot::plot_grid(
      plots$survival$plot,
      plots$response,
      ncol = 2,
      rel_widths = c(1.2, 1)
    )
    
    if (!is.null(save_path)) {
      ggsave(save_path, combined_plot, width = 14, height = 6, dpi = 300)
      cat("Plot saved to:", save_path, "\n")
    }
    
    return(combined_plot)
  } else if (length(plots) == 1) {
    if (!is.null(save_path)) {
      ggsave(save_path, plots[[1]], width = 8, height = 6, dpi = 300)
      cat("Plot saved to:", save_path, "\n")
    }
    
    return(plots[[1]])
  }
}

#' Analyze multiple immunotherapy cohorts
#' @param ml_results Machine learning model results
#' @param immunotherapy_datasets List of immunotherapy datasets
#' @param model_name Model to use
#' @return Combined analysis results
analyze_immunotherapy_cohorts <- function(ml_results,
                                         immunotherapy_datasets,
                                         model_name = "StepCox[both] + Ridge") {
  
  cat("\n========== Multi-Cohort Immunotherapy Analysis ==========\n")
  
  results <- list()
  
  for (dataset_name in names(immunotherapy_datasets)) {
    cat("\nAnalyzing:", dataset_name, "\n")
    cat("----------------------------------------\n")
    
    # Get risk scores
    risk_scores <- ml_results$riskscore[[model_name]][[dataset_name]]
    
    # Add risk groups
    risk_scores$Risk <- ifelse(
      risk_scores$RS > median(risk_scores$RS),
      "High", "Low"
    )
    
    # Analyze response
    if ("Response" %in% colnames(immunotherapy_datasets[[dataset_name]])) {
      response_analysis <- predict_immunotherapy_response(
        risk_scores,
        immunotherapy_datasets[[dataset_name]][, c("Response", "OS", "OS.time")]
      )
      
      results[[dataset_name]] <- response_analysis
    } else {
      # Just survival analysis
      surv_result <- survival_analysis(risk_scores)
      results[[dataset_name]] <- surv_result
    }
  }
  
  cat("==========================================================\n\n")
  
  return(results)
}
