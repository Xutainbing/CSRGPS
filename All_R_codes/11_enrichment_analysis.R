################################################################################
# Script: 11_enrichment_analysis.R
# Description: Differential expression and pathway enrichment analysis
################################################################################

library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA)
library(msigdbr)
library(enrichplot)
library(GseaVis)

#' Perform differential expression analysis
#' @param expr_matrix Expression matrix
#' @param group_factor Group factor
#' @param adjust_method P-value adjustment method
#' @param logFC_cutoff Log fold change cutoff
#' @param p_cutoff P-value cutoff
#' @return DEG results
differential_expression <- function(expr_matrix,
                                  group_factor,
                                  adjust_method = "BH",
                                  logFC_cutoff = 1,
                                  p_cutoff = 0.05) {
  
  cat("\n========== Differential Expression Analysis ==========\n")
  cat("Groups:", unique(group_factor), "\n")
  cat("Comparison:", paste(unique(group_factor), collapse = " vs "), "\n")
  
  # Create design matrix
  design <- model.matrix(~0 + factor(group_factor))
  colnames(design) <- levels(factor(group_factor))
  rownames(design) <- colnames(expr_matrix)
  
  # Create contrast matrix
  contrast_matrix <- makeContrasts(
    paste0(unique(group_factor), collapse = "-"),
    levels = design
  )
  
  # Fit linear model
  fit <- lmFit(expr_matrix, design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  # Get all results
  DEG <- topTable(fit2, coef = 1, n = Inf, adjust.method = adjust_method)
  DEG <- na.omit(DEG)
  
  # Add gene symbols
  DEG$symbol <- rownames(DEG)
  
  # Add change status
  DEG$change <- "NOT"
  DEG$change[DEG$P.Value < p_cutoff & DEG$logFC > logFC_cutoff] <- "UP"
  DEG$change[DEG$P.Value < p_cutoff & DEG$logFC < -logFC_cutoff] <- "DOWN"
  DEG$change <- factor(DEG$change, levels = c("UP", "DOWN", "NOT"))
  
  # Summary
  cat("\nDEG Summary:\n")
  cat("Up-regulated:", sum(DEG$change == "UP"), "\n")
  cat("Down-regulated:", sum(DEG$change == "DOWN"), "\n")
  cat("Total significant:", sum(DEG$change != "NOT"), "\n")
  cat("======================================================\n\n")
  
  return(DEG)
}

#' Convert gene symbols to Entrez IDs
#' @param gene_symbols Gene symbols
#' @param org_db Organism database
#' @return Conversion results
convert_to_entrez <- function(gene_symbols,
                             org_db = org.Hs.eg.db) {
  
  cat("Converting gene symbols to Entrez IDs\n")
  
  conversion <- bitr(
    unique(gene_symbols),
    fromType = "SYMBOL",
    toType = c("ENTREZID", "ENSEMBL"),
    OrgDb = org_db
  )
  
  cat("Converted:", nrow(conversion), "/", length(unique(gene_symbols)), "genes\n")
  
  return(conversion)
}

#' Perform GO enrichment analysis
#' @param DEG_results DEG results with Entrez IDs
#' @param ont GO ontology (BP, CC, MF, or ALL)
#' @param pvalue_cutoff P-value cutoff
#' @param qvalue_cutoff Q-value cutoff
#' @return GO enrichment results
go_enrichment <- function(DEG_results,
                         ont = "BP",
                         pvalue_cutoff = 0.05,
                         qvalue_cutoff = 0.05) {
  
  cat("Performing GO enrichment analysis (", ont, ")\n")
  
  # Add Entrez IDs if not present
  if (!"ENTREZID" %in% colnames(DEG_results)) {
    entrez_conversion <- convert_to_entrez(DEG_results$symbol)
    DEG_results <- merge(DEG_results, entrez_conversion,
                        by.x = "symbol", by.y = "SYMBOL")
  }
  
  # Get gene lists
  gene_up <- DEG_results$ENTREZID[DEG_results$change == "UP"]
  gene_down <- DEG_results$ENTREZID[DEG_results$change == "DOWN"]
  gene_diff <- c(gene_up, gene_down)
  gene_all <- DEG_results$ENTREZID
  
  # GO enrichment
  go_results <- enrichGO(
    gene = gene_diff,
    universe = gene_all,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff
  )
  
  cat("Enriched GO terms:", nrow(go_results@result), "\n")
  
  return(go_results)
}

#' Perform KEGG enrichment analysis
#' @param DEG_results DEG results with Entrez IDs
#' @param organism Organism code
#' @param pvalue_cutoff P-value cutoff
#' @param qvalue_cutoff Q-value cutoff
#' @return KEGG enrichment results
kegg_enrichment <- function(DEG_results,
                           organism = "hsa",
                           pvalue_cutoff = 0.05,
                           qvalue_cutoff = 0.05) {
  
  cat("Performing KEGG enrichment analysis\n")
  
  # Add Entrez IDs if not present
  if (!"ENTREZID" %in% colnames(DEG_results)) {
    entrez_conversion <- convert_to_entrez(DEG_results$symbol)
    DEG_results <- merge(DEG_results, entrez_conversion,
                        by.x = "symbol", by.y = "SYMBOL")
  }
  
  # Get gene lists
  gene_diff <- DEG_results$ENTREZID[DEG_results$change != "NOT"]
  gene_all <- DEG_results$ENTREZID
  
  # KEGG enrichment
  kegg_results <- enrichKEGG(
    gene = gene_diff,
    universe = gene_all,
    organism = organism,
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff,
    pAdjustMethod = "BH"
  )
  
  cat("Enriched KEGG pathways:", nrow(kegg_results@result), "\n")
  
  return(kegg_results)
}

#' Perform GSEA analysis
#' @param DEG_results DEG results
#' @param gene_sets Gene sets for GSEA
#' @param nperm Number of permutations
#' @return GSEA results
perform_gsea <- function(DEG_results,
                        gene_sets = NULL,
                        nperm = 1000) {
  
  cat("\n========== GSEA Analysis ==========\n")
  
  # Create ranked gene list
  gene_list <- DEG_results$logFC
  names(gene_list) <- toupper(DEG_results$symbol)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Load gene sets if not provided
  if (is.null(gene_sets)) {
    # Load Hallmark gene sets
    m_df <- msigdbr(species = "Homo sapiens", category = "H")
    gene_sets <- split(m_df$gene_symbol, m_df$gs_name)
  }
  
  # Run GSEA
  gsea_results <- GSEA(
    gene_list,
    TERM2GENE = data.frame(
      term = rep(names(gene_sets), lengths(gene_sets)),
      gene = unlist(gene_sets)
    ),
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    nPerm = nperm
  )
  
  # Get significant pathways
  up_pathways <- gsea_results@result[gsea_results@result$NES > 1 & 
                                    gsea_results@result$p.adjust < 0.05, ]
  down_pathways <- gsea_results@result[gsea_results@result$NES < -1 & 
                                      gsea_results@result$p.adjust < 0.05, ]
  
  cat("Up-regulated pathways:", nrow(up_pathways), "\n")
  cat("Down-regulated pathways:", nrow(down_pathways), "\n")
  cat("====================================\n\n")
  
  return(gsea_results)
}

#' Perform GSVA analysis
#' @param expr_matrix Expression matrix
#' @param gene_sets Gene sets
#' @param group_factor Group factor
#' @param method GSVA method
#' @return GSVA results
perform_gsva <- function(expr_matrix,
                        gene_sets = NULL,
                        group_factor,
                        method = "gsva") {
  
  cat("Performing GSVA analysis\n")
  
  # Load gene sets if not provided
  if (is.null(gene_sets)) {
    m_df <- msigdbr(species = "Homo sapiens", category = "H")
    gene_sets <- split(m_df$gene_symbol, m_df$gs_name)
  }
  
  # Run GSVA
  gsva_scores <- gsva(
    expr_matrix,
    gene_sets,
    method = method,
    kcdf = "Gaussian",
    parallel.sz = 1
  )
  
  # Differential analysis of GSVA scores
  design <- model.matrix(~0 + factor(group_factor))
  colnames(design) <- levels(factor(group_factor))
  contrast_matrix <- makeContrasts(
    paste0(unique(group_factor), collapse = "-"),
    levels = design
  )
  
  fit <- lmFit(gsva_scores, design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  diff_pathways <- topTable(fit2, coef = 1, n = Inf)
  
  return(list(
    scores = gsva_scores,
    differential = diff_pathways
  ))
}

#' Visualize enrichment results
#' @param enrichment_results Enrichment results
#' @param plot_type Type of plot
#' @param top_n Number of top terms to show
#' @param save_path Path to save plot
#' @return Plot object
plot_enrichment <- function(enrichment_results,
                          plot_type = c("barplot", "dotplot", "network", "heatmap"),
                          top_n = 20,
                          save_path = NULL) {
  
  plot_type <- match.arg(plot_type)
  
  # Select top terms
  if (nrow(enrichment_results@result) > top_n) {
    enrichment_results@result <- enrichment_results@result[1:top_n, ]
  }
  
  # Create plot based on type
  if (plot_type == "barplot") {
    p <- barplot(enrichment_results, showCategory = top_n) +
      theme(axis.text.y = element_text(size = 10))
  } else if (plot_type == "dotplot") {
    p <- dotplot(enrichment_results, showCategory = top_n) +
      theme(axis.text.y = element_text(size = 10))
  } else if (plot_type == "network") {
    p <- cnetplot(enrichment_results, categorySize = "pvalue",
                 foldChange = TRUE, circular = TRUE)
  } else if (plot_type == "heatmap") {
    p <- heatplot(enrichment_results, foldChange = TRUE)
  }
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 10, height = 8, dpi = 300)
    cat("Plot saved to:", save_path, "\n")
  }
  
  return(p)
}

#' Create GSEA plot
#' @param gsea_results GSEA results
#' @param pathway_ids Pathway IDs to plot
#' @param save_path Path to save plot
#' @return GSEA plot
plot_gsea <- function(gsea_results,
                     pathway_ids,
                     save_path = NULL) {
  
  library(GseaVis)
  
  # Create GSEA plot
  p <- gseaNb(
    object = gsea_results,
    geneSetID = pathway_ids,
    subPlot = 2,
    termWidth = 30,
    rmHt = TRUE,
    curveCol = brewer.pal(length(pathway_ids), "Set1")
  )
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 12, height = 8 * length(pathway_ids), dpi = 300)
    cat("GSEA plot saved to:", save_path, "\n")
  }
  
  return(p)
}
