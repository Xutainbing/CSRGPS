################################################################################
# Script: 06_main_pipeline.R
# Description: Main execution pipeline for CSRGPS analysis
################################################################################

# Source all function scripts
source("00_setup.R")
source("01_data_loading.R")
source("02_consensus_clustering.R")
source("03_wgcna_analysis.R")
source("04_ml_model.R")
source("05_survival_analysis.R")

# Create log file
log_file <- file.path("results", paste0("analysis_log_", Sys.Date(), ".txt"))
sink(log_file, split = TRUE)

cat("================================================================================\n")
cat("CSRGPS Analysis Pipeline\n")
cat("Started at:", format(Sys.time()), "\n")
cat("================================================================================\n\n")

# ==============================================================================
# STEP 1: DATA LOADING
# ==============================================================================

cat("\n>>> STEP 1: DATA LOADING\n")
cat("----------------------------------------\n")

# Load TCGA data
tcga_data <- load_tcga_data()
tpm <- tcga_data$expression
cl <- tcga_data$clinical

# Load senescence-related genes
srg_genes <- load_srg_genes()

# Filter expression matrix for SRG genes
tpm_srg <- tpm[rownames(tpm) %in% srg_genes, ]
cat("SRG expression matrix:", nrow(tpm_srg), "genes,", ncol(tpm_srg), "samples\n")

# ==============================================================================
# STEP 2: CONSENSUS CLUSTERING
# ==============================================================================

cat("\n>>> STEP 2: CONSENSUS CLUSTERING\n")
cat("----------------------------------------\n")

# Normalize expression for clustering
tpm_srg_norm <- normalize_for_clustering(tpm_srg)

# Perform consensus clustering
clustering_results <- perform_consensus_clustering(
  tpm_srg_norm,
  max_k = 9,
  output_dir = "results/consensus/"
)

# Add cluster assignments to clinical data
cl <- merge(cl, clustering_results$cluster_assignments,
            by.x = "ID", by.y = "Sample", all.x = TRUE)

# Plot PAC scores
pac_plot <- plot_pac_scores(
  clustering_results$pac_scores,
  "figures/pac_scores.pdf"
)

# Plot consensus heatmap
consensus_heatmap <- plot_consensus_heatmap(
  clustering_results$consensus_matrix,
  clustering_results$cluster_assignments,
  "figures/consensus_heatmap.pdf"
)

# PCA visualization
pca_plot <- perform_cluster_pca(
  tpm_srg_norm,
  clustering_results$cluster_assignments,
  "figures/cluster_pca.pdf"
)

# Survival analysis by cluster
cluster_surv <- survival_analysis(
  cl,
  group_var = "Cluster",
  time_var = "OS.time",
  event_var = "OS"
)

cluster_surv_plot <- create_survival_plot(
  cluster_surv,
  title = "Survival by Consensus Cluster",
  palette = c("#db6968", "#2E9FDF"),
  save_path = "figures/cluster_survival.pdf"
)

# ==============================================================================
# STEP 3: WGCNA ANALYSIS
# ==============================================================================

cat("\n>>> STEP 3: WGCNA ANALYSIS\n")
cat("----------------------------------------\n")

# Prepare trait data
trait_data <- cl[, c("Age", "Gender", "Type", "Cluster")]

# Run WGCNA pipeline
wgcna_results <- run_wgcna_pipeline(
  tpm_srg,
  trait_data,
  output_dir = "results/wgcna/"
)

# Extract hub genes from significant modules
hub_genes_blue <- extract_hub_genes(
  wgcna_results$data_mat,
  wgcna_results$network,
  "blue",
  trait_data[, "Cluster", drop = FALSE]
)

hub_genes_brown <- extract_hub_genes(
  wgcna_results$data_mat,
  wgcna_results$network,
  "brown",
  trait_data[, "Cluster", drop = FALSE]
)

# Combine WGCNA genes
wgcna_genes <- unique(c(hub_genes_blue$hub_genes, hub_genes_brown$hub_genes))
cat("Total WGCNA hub genes:", length(wgcna_genes), "\n")

# ==============================================================================
# STEP 4: GENE SELECTION
# ==============================================================================

cat("\n>>> STEP 4: GENE SELECTION\n")
cat("----------------------------------------\n")

# Univariate Cox analysis
cox_results <- univariate_cox_analysis(
  cbind(cl[, c("ID", "OS.time", "OS")], t(tpm_srg)),
  p_cutoff = 0.05
)

# Combine gene lists (intersection of WGCNA and Cox significant genes)
selected_genes <- intersect(wgcna_genes, cox_results$significant_genes)
cat("Final selected genes:", length(selected_genes), "\n")

# Save gene list
write.csv(
  data.frame(Gene = selected_genes),
  "results/selected_genes.csv",
  row.names = FALSE
)

# ==============================================================================
# STEP 5: MACHINE LEARNING MODEL
# ==============================================================================

cat("\n>>> STEP 5: MACHINE LEARNING MODEL\n")
cat("----------------------------------------\n")

# Load all datasets
all_datasets <- load_all_datasets()

# Prepare datasets for ML
ml_datasets <- prepare_ml_datasets(all_datasets, selected_genes)

# Build ML models
ml_results <- build_ml_models(
  train_data = ml_datasets$TCGA,
  all_datasets = ml_datasets,
  gene_list = selected_genes,
  mode = "all",
  seed = 123
)

# Save ML results
saveRDS(ml_results, "results/ml_models/ml_results.rds")

# Compare model performance
model_performance <- compare_model_performance(ml_results, ml_datasets)

# Extract risk scores for best model
best_model <- "StepCox[both] + Ridge"
risk_scores <- extract_risk_scores(ml_results, best_model)

# Calculate C-index
cindex_results <- calculate_cindex(risk_scores)

# ==============================================================================
# STEP 6: SURVIVAL ANALYSIS
# ==============================================================================

cat("\n>>> STEP 6: SURVIVAL ANALYSIS\n")
cat("----------------------------------------\n")

# Survival analysis for each dataset
for (dataset_name in names(risk_scores)) {
  cat("\nDataset:", dataset_name, "\n")
  
  # Survival analysis
  surv_result <- survival_analysis(
    risk_scores[[dataset_name]],
    group_var = "Risk"
  )
  
  # Create survival plot
  surv_plot <- create_survival_plot(
    surv_result,
    title = paste(dataset_name, "- CSRGPS Risk Groups"),
    save_path = paste0("figures/survival_", dataset_name, ".pdf")
  )
  
  # Time-dependent ROC
  roc_result <- calculate_time_roc(
    risk_scores[[dataset_name]],
    time_points = c(1, 3, 5)
  )
  
  # Plot ROC
  roc_plot <- plot_time_roc(
    roc_result,
    save_path = paste0("figures/roc_", dataset_name, ".pdf")
  )
}

# ==============================================================================
# STEP 7: META-ANALYSIS
# ==============================================================================

cat("\n>>> STEP 7: META-ANALYSIS\n")
cat("----------------------------------------\n")

# Combine all risk scores for meta-analysis
meta_data <- data.frame()

for (dataset_name in names(risk_scores)) {
  dataset_risk <- risk_scores[[dataset_name]]
  dataset_risk$Dataset <- dataset_name
  meta_data <- rbind(meta_data, dataset_risk)
}

# Overall survival analysis
meta_surv <- survival_analysis(meta_data, group_var = "Risk")

meta_surv_plot <- create_survival_plot(
  meta_surv,
  title = "Meta-Analysis - CSRGPS Risk Groups",
  save_path = "figures/survival_meta.pdf"
)

# ==============================================================================
# SAVE FINAL RESULTS
# ==============================================================================

cat("\n>>> SAVING FINAL RESULTS\n")
cat("----------------------------------------\n")

# Save all results
final_results <- list(
  clustering = clustering_results,
  wgcna = wgcna_results,
  selected_genes = selected_genes,
  ml_model = ml_results,
  risk_scores = risk_scores,
  cindex = cindex_results,
  model_performance = model_performance
)

saveRDS(final_results, "results/final_csrgps_results.rds")

# ==============================================================================
# COMPLETION
# ==============================================================================

cat("\n================================================================================\n")
cat("CSRGPS Analysis Pipeline Completed\n")
cat("Ended at:", format(Sys.time()), "\n")
cat("Results saved in: results/\n")
cat("Figures saved in: figures/\n")
cat("================================================================================\n")

# Close log
sink()
