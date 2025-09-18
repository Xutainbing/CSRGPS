################################################################################
# Script: 12_main_complete_pipeline.R
# Description: Complete CSRGPS analysis pipeline with all modules
################################################################################

# Source all scripts
source("00_setup.R")
source("01_data_loading.R")
source("02_consensus_clustering.R")
source("03_wgcna_analysis.R")
source("04_ml_model.R")
source("05_survival_analysis.R")
source("06_main_pipeline.R")  # Basic pipeline
source("07_immune_infiltration.R")
source("08_single_cell_analysis.R")
source("09_immunotherapy_analysis.R")
source("10_nomogram_analysis.R")
source("11_enrichment_analysis.R")

# Start logging
log_file <- file.path("results", paste0("complete_analysis_log_", Sys.Date(), ".txt"))
sink(log_file, split = TRUE)

cat("================================================================================\n")
cat("Complete CSRGPS Analysis Pipeline\n")
cat("Started at:", format(Sys.time()), "\n")
cat("================================================================================\n\n")

# ==============================================================================
# Run basic pipeline first (Steps 1-7 from 06_main_pipeline.R)
# ==============================================================================

cat("\n>>> RUNNING BASIC PIPELINE\n")
cat("----------------------------------------\n")
source("06_main_pipeline.R")

# Load saved results
final_results <- readRDS("results/final_csrgps_results.rds")

# ==============================================================================
# STEP 8: IMMUNE INFILTRATION ANALYSIS
# ==============================================================================

cat("\n>>> STEP 8: IMMUNE INFILTRATION ANALYSIS\n")
cat("----------------------------------------\n")

# Perform immune deconvolution
immune_results <- immune_deconvolution(
  tpm,
  methods = c("cibersort", "epic", "mcpcounter", "xcell", "estimate", "quantiseq")
)

# Compare between risk groups
immune_comparison <- compare_immune_by_risk(
  immune_results$combined,
  final_results$risk_scores$TCGA
)

# Immune checkpoint analysis
checkpoint_analysis <- analyze_immune_checkpoints(
  tpm,
  risk_groups = final_results$risk_scores$TCGA
)

# Calculate IPS
ips_scores <- calculate_ips(tpm)

# ssGSEA for immune signatures
immune_ssgsea <- immune_ssgsea(tpm)

# Visualizations
immune_heatmap <- plot_immune_heatmap(
  immune_results$combined,
  cl,
  save_path = "figures/immune_heatmap.pdf"
)

immune_boxplot <- plot_immune_boxplot(
  immune_results$individual$estimate,
  final_results$risk_scores$TCGA,
  save_path = "figures/immune_scores.pdf"
)

# ==============================================================================
# STEP 9: SINGLE-CELL ANALYSIS
# ==============================================================================

cat("\n>>> STEP 9: SINGLE-CELL ANALYSIS\n")
cat("----------------------------------------\n")

# Load single-cell data
sc_data <- load_single_cell_data("data/single_cell_glioma.rds")

# Integrate single-cell data
sc_integrated <- integrate_single_cell(sc_data)

# Perform Scissor analysis
scissor_results <- perform_scissor(
  tpm,
  sc_integrated,
  cl$Risk,
  alpha_range = seq(0.001, 0.01, 0.002)
)

# Visualize Scissor results
scissor_plot <- plot_scissor_results(
  scissor_results$sc_data,
  save_path = "figures/scissor_results.pdf"
)

# Trajectory analysis
trajectory_cds <- trajectory_analysis(scissor_results$sc_data)

# Find trajectory genes
trajectory_genes <- find_trajectory_genes(trajectory_cds)

# CytoTRACE analysis
cytotrace_results <- perform_cytotrace(scissor_results$sc_data)

# InferCNV analysis (if reference cells available)
if ("celltype" %in% colnames(sc_data@meta.data)) {
  reference_cells <- rownames(sc_data@meta.data)[
    sc_data@meta.data$celltype %in% c("T_cells", "Endothelial")
  ]
  
  if (length(reference_cells) > 0) {
    infercnv_results <- perform_infercnv(
      sc_data,
      reference_cells,
      "data/gene_order_file.txt",
      output_dir = "results/infercnv/"
    )
    
    # Calculate CNV scores
    cnv_scores <- calculate_cnv_scores(
      "results/infercnv/infercnv.observations.txt"
    )
  }
}

# ==============================================================================
# STEP 10: IMMUNOTHERAPY ANALYSIS
# ==============================================================================

cat("\n>>> STEP 10: IMMUNOTHERAPY ANALYSIS\n")
cat("----------------------------------------\n")

# Load immunotherapy datasets
immunotherapy_datasets <- list()

# Example: Load multiple immunotherapy cohorts
immunotherapy_datasets$PRJNA482620 <- load_immunotherapy_dataset(
  "PRJNA482620",
  "data/immunotherapy/PRJNA482620_exp.rds",
  "data/immunotherapy/PRJNA482620_clin.tsv"
)

immunotherapy_datasets$IMvigor210 <- load_immunotherapy_dataset(
  "IMvigor210",
  "data/immunotherapy/IMvigor210_exp.rds",
  "data/immunotherapy/IMvigor210_clin.tsv"
)

# Analyze immunotherapy cohorts
immunotherapy_analysis <- analyze_immunotherapy_cohorts(
  final_results$ml_model,
  immunotherapy_datasets
)

# TIDE analysis
tide_input <- calculate_tide(tpm, "results/tide_input.txt")

# After getting TIDE results from web server
if (file.exists("results/tide_results.csv")) {
  tide_analysis <- analyze_tide_results(
    "results/tide_results.csv",
    final_results$risk_scores$TCGA
  )
}

# Visualize immunotherapy response
for (cohort in names(immunotherapy_analysis)) {
  plot_immunotherapy_response(
    immunotherapy_analysis[[cohort]]$merged_data,
    plot_type = "both",
    save_path = paste0("figures/immunotherapy_", cohort, ".pdf")
  )
}

# ==============================================================================
# STEP 11: NOMOGRAM CONSTRUCTION
# ==============================================================================

cat("\n>>> STEP 11: NOMOGRAM CONSTRUCTION\n")
cat("----------------------------------------\n")

# Prepare data for nomogram
nomogram_data <- merge(
  cl,
  final_results$risk_scores$TCGA,
  by.x = "ID",
  by.y = "row.names"
)

# Build nomogram
nomogram_results <- build_nomogram(
  nomogram_data,
  variables = c("Age", "Type", "IDH_mutation", "CSRGPS" = "RS"),
  time_points = c(1, 3, 5)
)

# Create interactive nomogram
interactive_nom <- create_interactive_nomogram(
  nomogram_results$model,
  nomogram_results$data,
  save_path = "figures/interactive_nomogram.pdf"
)

# Validate nomogram
calibration_results <- validate_nomogram(
  nomogram_results,
  save_path = "figures/nomogram_calibration.pdf"
)

# Calculate discrimination
discrimination_results <- calculate_nomogram_discrimination(nomogram_results)

# Decision curve analysis
dca_results <- decision_curve_analysis(
  nomogram_results,
  save_path = "figures/decision_curve.pdf"
)

# Clinical impact curve
clinical_impact <- clinical_impact_curve(
  nomogram_results,
  save_path = "figures/clinical_impact.pdf"
)

# ==============================================================================
# STEP 12: DIFFERENTIAL EXPRESSION AND ENRICHMENT
# ==============================================================================

cat("\n>>> STEP 12: DIFFERENTIAL EXPRESSION AND ENRICHMENT\n")
cat("----------------------------------------\n")

# Differential expression analysis
DEG_results <- differential_expression(
  tpm,
  final_results$risk_scores$TCGA$Risk
)

# GO enrichment
go_bp <- go_enrichment(DEG_results, ont = "BP")
go_cc <- go_enrichment(DEG_results, ont = "CC")
go_mf <- go_enrichment(DEG_results, ont = "MF")

# KEGG enrichment
kegg_results <- kegg_enrichment(DEG_results)

# GSEA analysis
gsea_results <- perform_gsea(DEG_results)

# GSVA analysis
gsva_results <- perform_gsva(
  tpm,
  group_factor = final_results$risk_scores$TCGA$Risk
)

# Visualizations
plot_enrichment(go_bp, "dotplot", save_path = "figures/go_bp_dotplot.pdf")
plot_enrichment(kegg_results, "barplot", save_path = "figures/kegg_barplot.pdf")

# GSEA plot for top pathways
top_pathways <- gsea_results@result$ID[1:5]
plot_gsea(gsea_results, top_pathways, save_path = "figures/gsea_top_pathways.pdf")

# ==============================================================================
# STEP 13: SAVE ALL RESULTS
# ==============================================================================

cat("\n>>> SAVING ALL RESULTS\n")
cat("----------------------------------------\n")

# Compile all results
complete_results <- list(
  # Basic pipeline results
  clustering = final_results$clustering,
  wgcna = final_results$wgcna,
  selected_genes = final_results$selected_genes,
  ml_model = final_results$ml_model,
  risk_scores = final_results$risk_scores,
  
  # Extended analyses
  immune = list(
    deconvolution = immune_results,
    comparison = immune_comparison,
    checkpoints = checkpoint_analysis,
    ips = ips_scores,
    ssgsea = immune_ssgsea
  ),
  
  single_cell = list(
    scissor = scissor_results,
    trajectory = trajectory_cds,
    trajectory_genes = trajectory_genes,
    cytotrace = cytotrace_results
  ),
  
  immunotherapy = immunotherapy_analysis,
  
  nomogram = list(
    model = nomogram_results,
    calibration = calibration_results,
    discrimination = discrimination_results,
    dca = dca_results
  ),
  
  enrichment = list(
    DEG = DEG_results,
    GO = list(BP = go_bp, CC = go_cc, MF = go_mf),
    KEGG = kegg_results,
    GSEA = gsea_results,
    GSVA = gsva_results
  )
)

# Save complete results
saveRDS(complete_results, "results/complete_csrgps_analysis.rds")

# ==============================================================================
# GENERATE SUMMARY REPORT
# ==============================================================================

cat("\n>>> GENERATING SUMMARY REPORT\n")
cat("----------------------------------------\n")

# Create summary statistics
summary_stats <- data.frame(
  Analysis = c(
    "Total samples",
    "Senescence genes analyzed",
    "Optimal clusters",
    "WGCNA modules",
    "Selected prognostic genes",
    "ML models tested",
    "Best model C-index",
    "Immune cell types analyzed",
    "Single cells analyzed",
    "Scissor+ cells",
    "Scissor- cells",
    "DEGs identified",
    "Enriched GO terms",
    "Enriched KEGG pathways"
  ),
  Value = c(
    ncol(tpm),
    length(srg_genes),
    final_results$clustering$optimal_k,
    length(unique(final_results$wgcna$module_colors)),
    length(final_results$selected_genes),
    length(final_results$ml_model$riskscore),
    round(max(final_results$model_performance$summary$C_index), 3),
    ncol(immune_results$individual$cibersort) - 1,
    ncol(sc_integrated),
    length(scissor_results$selected_cells$positive),
    length(scissor_results$selected_cells$negative),
    sum(DEG_results$change != "NOT"),
    nrow(go_bp@result),
    nrow(kegg_results@result)
  )
)

write.csv(summary_stats, "results/analysis_summary.csv", row.names = FALSE)

cat("\nAnalysis Summary:\n")
print(summary_stats)

# ==============================================================================
# COMPLETION
# ==============================================================================

cat("\n================================================================================\n")
cat("Complete CSRGPS Analysis Pipeline Finished\n")
cat("Ended at:", format(Sys.time()), "\n")
cat("All results saved in: results/\n")
cat("All figures saved in: figures/\n")
cat("================================================================================\n")

# Close log
sink()
