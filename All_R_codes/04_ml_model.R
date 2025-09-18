################################################################################
# Script: 04_ml_model.R
# Description: Machine learning prognostic model development
################################################################################

#' Prepare datasets for machine learning
#' @param datasets List of datasets
#' @param gene_list Gene list to use
#' @return Prepared datasets with common genes
prepare_ml_datasets <- function(datasets, gene_list) {
  cat("\n========== Preparing ML Datasets ==========\n")
  
  # Find common genes across all datasets
  all_genes <- lapply(datasets, function(x) {
    colnames(x)[!colnames(x) %in% c("ID", "OS.time", "OS")]
  })
  
  common_genes <- Reduce(intersect, all_genes)
  common_genes <- intersect(common_genes, gene_list)
  
  cat("Number of common genes:", length(common_genes), "\n")
  
  # Keep only common genes and survival columns
  prepared_datasets <- lapply(datasets, function(x) {
    cols_to_keep <- c("ID", "OS.time", "OS", common_genes)
    cols_to_keep <- cols_to_keep[cols_to_keep %in% colnames(x)]
    return(x[, cols_to_keep])
  })
  
  cat("============================================\n\n")
  
  return(prepared_datasets)
}

#' Perform univariate Cox regression
#' @param data Dataset with survival and gene expression
#' @param p_cutoff P-value cutoff
#' @return Significant genes
univariate_cox_analysis <- function(data, p_cutoff = 0.05) {
  cat("Performing univariate Cox regression\n")
  
  # Get gene columns
  gene_cols <- colnames(data)[!colnames(data) %in% c("ID", "OS.time", "OS")]
  
  # Initialize results
  cox_results <- data.frame()
  
  # Progress counter
  pb <- txtProgressBar(min = 0, max = length(gene_cols), style = 3)
  
  for (i in seq_along(gene_cols)) {
    gene <- gene_cols[i]
    
    # Cox regression
    cox_formula <- as.formula(paste("Surv(OS.time, OS) ~", gene))
    cox_model <- coxph(cox_formula, data = data)
    cox_summary <- summary(cox_model)
    
    # Extract results
    result <- data.frame(
      Gene = gene,
      HR = cox_summary$conf.int[, "exp(coef)"],
      HR_lower = cox_summary$conf.int[, "lower .95"],
      HR_upper = cox_summary$conf.int[, "upper .95"],
      P_value = cox_summary$coefficients[, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
    
    cox_results <- rbind(cox_results, result)
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  # Filter significant genes
  significant_genes <- cox_results$Gene[cox_results$P_value < p_cutoff]
  
  cat("\nSignificant genes:", length(significant_genes), "/", length(gene_cols), "\n")
  
  return(list(
    results = cox_results,
    significant_genes = significant_genes
  ))
}

#' Build machine learning models using Mime1
#' @param train_data Training dataset
#' @param all_datasets List of all datasets for validation
#' @param gene_list Candidate gene list
#' @param mode Model building mode
#' @param seed Random seed
#' @return ML model results
build_ml_models <- function(train_data,
                           all_datasets,
                           gene_list,
                           mode = "all",
                           seed = 123) {
  
  cat("\n========== Building ML Models ==========\n")
  cat("Mode:", mode, "\n")
  cat("Training dataset samples:", nrow(train_data), "\n")
  cat("Candidate genes:", length(gene_list), "\n")
  
  # Build models
  ml_results <- ML.Dev.Prog.Sig(
    train_data = train_data,
    list_train_vali_Data = all_datasets,
    unicox.filter.for.candi = TRUE,
    unicox_p_cutoff = 0.05,
    candidate_genes = gene_list,
    mode = mode,
    nodesize = 5,
    seed = seed
  )
  
  cat("=========================================\n\n")
  
  return(ml_results)
}

#' Extract risk scores from ML results
#' @param ml_results ML model results
#' @param model_name Specific model to extract
#' @return Risk scores for each dataset
extract_risk_scores <- function(ml_results, 
                               model_name = "StepCox[both] + Ridge") {
  
  cat("Extracting risk scores for model:", model_name, "\n")
  
  risk_scores <- ml_results[["riskscore"]][[model_name]]
  
  # Add risk groups based on median
  for (dataset in names(risk_scores)) {
    risk_scores[[dataset]]$Risk <- ifelse(
      risk_scores[[dataset]]$RS > median(risk_scores[[dataset]]$RS),
      "High",
      "Low"
    )
  }
  
  return(risk_scores)
}

#' Calculate C-index for model validation
#' @param risk_scores Risk score results
#' @return C-index for each dataset
calculate_cindex <- function(risk_scores) {
  cat("Calculating C-index for validation\n")
  
  cindex_results <- data.frame()
  
  for (dataset in names(risk_scores)) {
    # Cox model with risk score
    cox_model <- coxph(
      Surv(OS.time, OS) ~ RS,
      data = risk_scores[[dataset]]
    )
    
    # Get C-index
    cindex <- summary(cox_model)$concordance[1]
    
    result <- data.frame(
      Dataset = dataset,
      C_index = cindex,
      N = nrow(risk_scores[[dataset]]),
      stringsAsFactors = FALSE
    )
    
    cindex_results <- rbind(cindex_results, result)
  }
  
  print(cindex_results)
  
  return(cindex_results)
}

#' Compare model performance
#' @param ml_results ML model results
#' @param datasets List of datasets
#' @return Performance comparison results
compare_model_performance <- function(ml_results, datasets) {
  cat("\n========== Model Performance Comparison ==========\n")
  
  # Get all model names
  model_names <- names(ml_results$riskscore)
  
  # Calculate C-index for each model
  performance_results <- list()
  
  for (model in model_names) {
    risk_scores <- extract_risk_scores(ml_results, model)
    cindex <- calculate_cindex(risk_scores)
    cindex$Model <- model
    performance_results[[model]] <- cindex
  }
  
  # Combine results
  combined_performance <- do.call(rbind, performance_results)
  
  # Calculate mean C-index across datasets
  mean_cindex <- aggregate(
    C_index ~ Model,
    data = combined_performance,
    FUN = mean
  )
  mean_cindex <- mean_cindex[order(mean_cindex$C_index, decreasing = TRUE), ]
  
  cat("\nMean C-index across datasets:\n")
  print(mean_cindex)
  
  cat("==================================================\n\n")
  
  return(list(
    detailed = combined_performance,
    summary = mean_cindex
  ))
}
