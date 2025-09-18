################################################################################
# Script: 10_nomogram_analysis.R
# Description: Nomogram construction and validation
################################################################################

library(rms)
library(regplot)
library(ggDCA)
library(rmda)

#' Build nomogram model
#' @param data Clinical data with risk scores
#' @param variables Variables to include
#' @param time_points Time points for survival probability
#' @return Nomogram model and plot
build_nomogram <- function(data,
                          variables = c("Age", "Type", "IDH_mutation", "CSRGPS"),
                          time_points = c(1, 3, 5)) {
  
  cat("\n========== Building Nomogram ==========\n")
  cat("Variables:", paste(variables, collapse = ", "), "\n")
  cat("Time points:", paste(time_points, "years"), "\n")
  
  # Prepare data
  nomogram_data <- data[, c("OS.time", "OS", variables)]
  
  # Convert time to years if needed
  if (max(nomogram_data$OS.time, na.rm = TRUE) > 100) {
    nomogram_data$OS.time <- nomogram_data$OS.time / 365
  }
  
  # Convert categorical variables to factors
  for (var in variables) {
    if (is.character(nomogram_data[[var]])) {
      nomogram_data[[var]] <- factor(nomogram_data[[var]])
    }
  }
  
  # Remove missing values
  nomogram_data <- na.omit(nomogram_data)
  cat("Samples used:", nrow(nomogram_data), "\n")
  
  # Set up data distribution
  dd <- datadist(nomogram_data)
  options(datadist = "dd")
  
  # Build Cox model
  formula_str <- paste("Surv(OS.time, OS) ~", paste(variables, collapse = " + "))
  cox_model <- cph(
    as.formula(formula_str),
    data = nomogram_data,
    x = TRUE,
    y = TRUE,
    surv = TRUE
  )
  
  # Print model summary
  print(cox_model)
  
  # Create survival functions
  surv <- Survival(cox_model)
  surv_functions <- list()
  for (t in time_points) {
    surv_functions[[paste0(t, "_year")]] <- function(x, time = t) surv(time, x)
  }
  
  # Create nomogram
  nom <- nomogram(
    cox_model,
    fun = surv_functions,
    lp = FALSE,
    funlabel = paste0(time_points, "-year survival"),
    fun.at = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
  )
  
  # Plot nomogram
  plot(nom)
  
  cat("========================================\n\n")
  
  return(list(
    model = cox_model,
    nomogram = nom,
    data = nomogram_data
  ))
}

#' Create interactive nomogram
#' @param cox_model Cox model from nomogram
#' @param data Nomogram data
#' @param observation_index Index of observation to highlight
#' @param save_path Path to save plot
#' @return Interactive nomogram plot
create_interactive_nomogram <- function(cox_model,
                                      data,
                                      observation_index = 1,
                                      save_path = NULL) {
  
  library(regplot)
  
  # Create interactive plot
  p <- regplot(
    cox_model,
    plots = c("density", "boxes"),
    observation = data[observation_index, ],
    points = TRUE,
    droplines = TRUE,
    title = "Interactive Survival Nomogram",
    dencol = "#EA5455",
    boxcol = "#002B5B",
    prfail = TRUE,
    failtime = c(1, 3, 5),
    clickable = FALSE  # Set to TRUE for interactive version
  )
  
  if (!is.null(save_path)) {
    pdf(save_path, width = 12, height = 10)
    print(p)
    dev.off()
    cat("Interactive nomogram saved to:", save_path, "\n")
  }
  
  return(p)
}

#' Validate nomogram with calibration curves
#' @param nomogram_results Nomogram results
#' @param time_points Time points to validate
#' @param method Validation method
#' @param save_path Path to save plot
#' @return Calibration results
validate_nomogram <- function(nomogram_results,
                            time_points = c(1, 3, 5),
                            method = "boot",
                            save_path = NULL) {
  
  cat("\n========== Nomogram Validation ==========\n")
  
  calibration_results <- list()
  
  # Set time units
  units(nomogram_results$data$OS.time) <- "Year"
  
  for (t in time_points) {
    cat("Calibrating for", t, "year survival\n")
    
    # Rebuild model for specific time point
    cox_t <- cph(
      Surv(OS.time, OS) ~ .,
      x = TRUE,
      y = TRUE,
      surv = TRUE,
      time.inc = t,
      data = nomogram_results$data[, -which(names(nomogram_results$data) %in% c("OS.time", "OS"))]
    )
    
    # Calibration
    cal <- calibrate(
      cox_t,
      cmethod = "KM",
      method = method,
      u = t,
      m = round(nrow(nomogram_results$data) / 3),
      B = 1000
    )
    
    calibration_results[[paste0(t, "_year")]] <- cal
  }
  
  # Plot calibration curves
  if (!is.null(save_path)) {
    pdf(save_path, width = 8, height = 8)
    
    # Set up colors
    colors <- c("#8AD879", "#FA9F42", "#F3533A")
    
    for (i in seq_along(time_points)) {
      cal <- calibration_results[[paste0(time_points[i], "_year")]]
      
      if (i == 1) {
        plot(cal,
             lty = 1,
             lwd = 1.5,
             xlim = c(0, 1),
             ylim = c(0, 1),
             xlab = "Predicted Probability",
             ylab = "Actual Probability",
             errbar.col = colors[i],
             col = colors[i],
             subtitles = FALSE)
      } else {
        plot(cal,
             lty = 1,
             lwd = 1.5,
             errbar.col = colors[i],
             col = colors[i],
             subtitles = FALSE,
             add = TRUE)
      }
    }
    
    # Add diagonal line
    abline(0, 1, lty = 2, lwd = 1, col = "gray")
    
    # Add legend
    legend("bottomright",
           legend = paste0(time_points, "-Year"),
           col = colors[1:length(time_points)],
           lty = 1,
           lwd = 1.5)
    
    dev.off()
    cat("Calibration curves saved to:", save_path, "\n")
  }
  
  cat("==========================================\n\n")
  
  return(calibration_results)
}

#' Calculate nomogram discrimination (C-index and time-ROC)
#' @param nomogram_results Nomogram results
#' @param time_points Time points for ROC
#' @return Discrimination metrics
calculate_nomogram_discrimination <- function(nomogram_results,
                                            time_points = c(1, 3, 5)) {
  
  cat("Calculating discrimination metrics\n")
  
  # C-index
  cox_summary <- summary(nomogram_results$model)
  c_index <- cox_summary$concordance[1]
  cat("C-index:", round(c_index, 3), "\n")
  
  # Time-dependent ROC
  roc_results <- list()
  auc_values <- numeric()
  
  for (t in time_points) {
    # Get predicted probabilities
    surv_prob <- 1 - summary(
      survfit(nomogram_results$model, newdata = nomogram_results$data),
      times = t
    )$surv
    
    # Calculate ROC
    roc <- survivalROC(
      Stime = nomogram_results$data$OS.time,
      status = nomogram_results$data$OS,
      marker = -surv_prob,  # Negative because higher risk = lower survival
      predict.time = t,
      method = "KM"
    )
    
    roc_results[[paste0(t, "_year")]] <- roc
    auc_values <- c(auc_values, roc$AUC)
    
    cat(t, "year AUC:", round(roc$AUC, 3), "\n")
  }
  
  return(list(
    c_index = c_index,
    roc_results = roc_results,
    auc_values = auc_values
  ))
}

#' Perform decision curve analysis
#' @param nomogram_results Nomogram results
#' @param time_point Time point for analysis
#' @param save_path Path to save plot
#' @return DCA results
decision_curve_analysis <- function(nomogram_results,
                                   time_point = 3,
                                   save_path = NULL) {
  
  cat("\n========== Decision Curve Analysis ==========\n")
  cat("Time point:", time_point, "years\n")
  
  library(ggDCA)
  
  # Prepare models for comparison
  models <- list()
  
  # Full model
  models$nomogram <- nomogram_results$model
  
  # Individual variable models
  for (var in names(nomogram_results$data)[-c(1, 2)]) {
    formula_str <- paste("Surv(OS.time, OS) ~", var)
    models[[var]] <- cph(
      as.formula(formula_str),
      data = nomogram_results$data,
      x = TRUE,
      y = TRUE,
      surv = TRUE
    )
  }
  
  # Perform DCA
  dca_result <- dca(
    models$nomogram,
    models$CSRGPS,
    models$Age,
    models$Type,
    model.names = names(models),
    times = time_point
  )
  
  # Plot DCA
  p <- ggplot(dca_result) +
    geom_line(aes(x = thresholds, y = NB, color = model, linetype = model), size = 1) +
    scale_color_brewer(palette = "Set2") +
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "twodash")) +
    labs(
      x = "Threshold Probability",
      y = "Net Benefit",
      title = paste("Decision Curve Analysis -", time_point, "Year Survival")
    ) +
    theme_bw() +
    theme(
      legend.position = c(0.8, 0.8),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 10, height = 8, dpi = 300)
    cat("DCA plot saved to:", save_path, "\n")
  }
  
  cat("==============================================\n\n")
  
  return(list(
    dca_result = dca_result,
    plot = p
  ))
}

#' Create clinical impact curve
#' @param nomogram_results Nomogram results
#' @param population_size Population size for impact calculation
#' @param time_point Time point for analysis
#' @param save_path Path to save plot
#' @return Clinical impact results
clinical_impact_curve <- function(nomogram_results,
                                 population_size = 1000,
                                 time_point = 3,
                                 save_path = NULL) {
  
  cat("Creating clinical impact curve\n")
  
  library(rmda)
  
  # Prepare data for rmda
  decision_data <- nomogram_results$data
  decision_data$OS_binary <- decision_data$OS
  
  # Decision curve
  baseline_model <- decision_curve(
    OS_binary ~ Age + Type + IDH_mutation + CSRGPS,
    data = decision_data,
    thresholds = seq(0, 1, by = 0.01),
    bootstraps = 1000
  )
  
  # Plot clinical impact
  p <- plot_clinical_impact(
    baseline_model,
    population.size = population_size,
    cost.benefit.axis = TRUE,
    n.cost.benefits = 8,
    col = c('red', 'green'),
    confidence.intervals = TRUE,
    ylim = c(0, population_size),
    legend.position = "bottomleft"
  )
  
  if (!is.null(save_path)) {
    pdf(save_path, width = 10, height = 8)
    print(p)
    dev.off()
    cat("Clinical impact curve saved to:", save_path, "\n")
  }
  
  return(p)
}
