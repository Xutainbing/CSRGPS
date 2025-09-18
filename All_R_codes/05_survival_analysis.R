################################################################################
# Script: 05_survival_analysis.R
# Description: Survival analysis functions
################################################################################

#' Perform survival analysis
#' @param data Data frame with survival data and risk groups
#' @param group_var Grouping variable name
#' @param time_var Time variable name
#' @param event_var Event variable name
#' @param time_unit Time unit for display
#' @return Survival analysis results
survival_analysis <- function(data, 
                             group_var = "Risk",
                             time_var = "OS.time",
                             event_var = "OS",
                             time_unit = "days") {
  
  cat("Performing survival analysis\n")
  
  # Convert time if needed
  if (time_unit == "years") {
    data[[time_var]] <- data[[time_var]] / 365
  } else if (time_unit == "months") {
    data[[time_var]] <- data[[time_var]] / 30
  }
  
  # Create survival formula
  surv_formula <- as.formula(paste("Surv(", time_var, ",", event_var, ") ~", group_var))
  
  # Fit survival model
  fit <- survfit(surv_formula, data = data)
  
  # Log-rank test
  survdiff_result <- survdiff(surv_formula, data = data)
  p_value <- 1 - pchisq(survdiff_result$chisq, length(survdiff_result$n) - 1)
  
  # Median survival
  median_survival <- summary(fit)$table[, "median"]
  
  cat("Log-rank test p-value:", format.pval(p_value), "\n")
  cat("Median survival:\n")
  print(median_survival)
  
  return(list(
    fit = fit,
    p_value = p_value,
    median_survival = median_survival,
    data = data
  ))
}

#' Create survival plot
#' @param surv_result Survival analysis results
#' @param title Plot title
#' @param palette Color palette
#' @param save_path Path to save plot
#' @return ggsurvplot object
create_survival_plot <- function(surv_result,
                                title = "",
                                palette = c("#982b2b", "#0074b3"),
                                save_path = NULL) {
  
  # Count samples per group
  group_counts <- table(surv_result$data$Risk)
  
  # Create survival plot
  p <- ggsurvplot(
    surv_result$fit,
    data = surv_result$data,
    palette = palette,
    risk.table = TRUE,
    risk.table.col = "strata",
    risk.table.height = 0.25,
    risk.table.title = "Number at risk",
    risk.table.y.text = FALSE,
    conf.int = TRUE,
    pval = TRUE,
    pval.size = 5,
    pval.coord = c(0.15, 0.1),
    legend = c(0.75, 0.85),
    legend.title = "",
    legend.labs = paste0(names(group_counts), " (n=", group_counts, ")"),
    xlab = "Time",
    ylab = "Overall Survival Probability",
    title = title,
    surv.median.line = "hv",
    ggtheme = theme_bw() + 
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
      ),
    font.x = 12,
    font.y = 12,
    font.legend = 10,
    font.tickslab = 10
  )
  
  if (!is.null(save_path)) {
    pdf(save_path, width = 8, height = 6)
    print(p)
    dev.off()
    cat("Survival plot saved to:", save_path, "\n")
  }
  
  return(p)
}

#' Calculate time-dependent ROC
#' @param data Data with risk scores and survival
#' @param time_points Time points for ROC (in years)
#' @param marker_var Marker variable name
#' @return Time-dependent ROC results
calculate_time_roc <- function(data,
                              time_points = c(1, 3, 5),
                              marker_var = "RS") {
  
  cat("Calculating time-dependent ROC\n")
  
  roc_results <- list()
  auc_values <- numeric()
  
  for (t in time_points) {
    roc <- survivalROC(
      Stime = data$OS.time,
      status = data$OS,
      marker = data[[marker_var]],
      predict.time = 365 * t,
      method = "KM"
    )
    
    roc_results[[paste0(t, "_year")]] <- roc
    auc_values <- c(auc_values, roc$AUC)
    
    cat(t, "year AUC:", round(roc$AUC, 3), "\n")
  }
  
  return(list(
    roc_objects = roc_results,
    auc_values = auc_values,
    time_points = time_points
  ))
}

#' Plot time-dependent ROC curves
#' @param roc_results Results from calculate_time_roc
#' @param save_path Path to save plot
#' @return ggplot object
plot_time_roc <- function(roc_results, save_path = NULL) {
  
  # Prepare data for plotting
  plot_data <- data.frame()
  
  for (i in seq_along(roc_results$time_points)) {
    t <- roc_results$time_points[i]
    roc_name <- paste0(t, "_year")
    
    df <- data.frame(
      FPR = as.numeric(roc_results$roc_objects[[roc_name]]$FP),
      TPR = as.numeric(roc_results$roc_objects[[roc_name]]$TP),
      Time = factor(paste0(t, " Year")),
      stringsAsFactors = FALSE
    )
    
    plot_data <- rbind(plot_data, df)
  }
  
  # Create labels with AUC
  auc_labels <- paste0(
    roc_results$time_points, " Year (AUC = ",
    round(roc_results$auc_values, 3), ")"
  )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = FPR, y = TPR, color = Time)) +
    geom_line(size = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    scale_color_manual(
      values = c("#8AD879", "#FA9F42", "#F3533A"),
      labels = auc_labels
    ) +
    theme_bw() +
    labs(
      x = "False Positive Rate",
      y = "True Positive Rate",
      title = "Time-Dependent ROC Curves"
    ) +
    theme(
      panel.grid = element_blank(),
      legend.position = c(0.7, 0.3),
      legend.background = element_rect(fill = "white", color = "black"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 8, height = 6, dpi = 300)
    cat("ROC plot saved to:", save_path, "\n")
  }
  
  return(p)
}

#' Perform subgroup survival analysis
#' @param data Data with risk scores and clinical variables
#' @param subgroup_var Subgroup variable name
#' @param risk_var Risk variable name
#' @return Subgroup analysis results
subgroup_survival_analysis <- function(data,
                                      subgroup_var,
                                      risk_var = "Risk") {
  
  cat("\nSubgroup survival analysis by", subgroup_var, "\n")
  
  subgroups <- unique(data[[subgroup_var]])
  subgroups <- subgroups[!is.na(subgroups)]
  
  results <- list()
  
  for (subgroup in subgroups) {
    cat("\nAnalyzing subgroup:", subgroup, "\n")
    
    # Filter data
    subgroup_data <- data[data[[subgroup_var]] == subgroup, ]
    
    if (nrow(subgroup_data) < 10) {
      cat("Insufficient samples (n =", nrow(subgroup_data), "), skipping\n")
      next
    }
    
    # Survival analysis
    surv_result <- survival_analysis(subgroup_data, group_var = risk_var)
    
    # Store results
    results[[as.character(subgroup)]] <- surv_result
  }
  
  return(results)
}

#' Forest plot for survival HR
#' @param data Data for analysis
#' @param variables Variables to include
#' @param reference_levels Reference levels for categorical variables
#' @return Forest plot
create_forest_plot <- function(data,
                              variables,
                              reference_levels = NULL) {
  
  cat("Creating forest plot for survival analysis\n")
  
  # Univariate Cox for each variable
  forest_data <- data.frame()
  
  for (var in variables) {
    # Skip if too many NAs
    if (sum(is.na(data[[var]])) > nrow(data) * 0.5) {
      next
    }
    
    # Cox regression
    cox_formula <- as.formula(paste("Surv(OS.time, OS) ~", var))
    cox_model <- coxph(cox_formula, data = data)
    cox_summary <- summary(cox_model)
    
    # Extract results
    result <- data.frame(
      Variable = var,
      N = sum(!is.na(data[[var]])),
      HR = cox_summary$conf.int[1, "exp(coef)"],
      HR_lower = cox_summary$conf.int[1, "lower .95"],
      HR_upper = cox_summary$conf.int[1, "upper .95"],
      P_value = cox_summary$coefficients[1, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
    
    forest_data <- rbind(forest_data, result)
  }
  
  # Order by HR
  forest_data <- forest_data[order(forest_data$HR), ]
  forest_data$Variable <- factor(forest_data$Variable, levels = forest_data$Variable)
  
  # Create plot
  p <- ggplot(forest_data, aes(x = HR, y = Variable)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
    scale_x_log10() +
    theme_bw() +
    labs(
      x = "Hazard Ratio (95% CI)",
      y = "",
      title = "Forest Plot - Univariate Cox Regression"
    ) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
  
  return(list(
    plot = p,
    data = forest_data
  ))
}
