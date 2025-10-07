#!/usr/bin/env Rscript

#' Wrapper Script for Running Separate Analyses
#'
#' This script reads config_separate_analyses.yaml and runs independent
#' analyses for different trimester comparisons, then compares the results
#'
#' @author Expression Integration Pipeline
#' @date 2025-10-05

cat("\n=== Separate Analyses Pipeline ===\n\n")

# Load required libraries
library(yaml)

# Source all modules
source("R/config.R")
source("R/utils.R")
source("R/01_data_merging.R")
source("R/02_batch_correction.R")
source("R/03_pca_analysis.R")
source("R/04_differential_expression.R")
source("R/05_network_clustering.R")
source("R/06_enrichment_analysis.R")
source("R/07_visualization.R")

# Load configuration
config <- load_config("config/config_separate_analyses.yaml")

# Initialize STRING database (shared across analyses)
cat("\n=== Initializing STRING Database (shared) ===\n")
string_db <- initialize_stringdb(
  version = config$stringdb$version,
  species = config$stringdb$species,
  score_threshold = config$stringdb$score_threshold,
  use_cache = config$stringdb$use_cache,
  cache_dir = config$stringdb$cache_dir
)

# Storage for results
all_results <- list()

# Loop through each separate analysis
for (analysis_id in names(config$separate_analyses)) {

  analysis <- config$separate_analyses[[analysis_id]]

  cat("\n\n")
  cat("=====================================================\n")
  cat("  Running Analysis:", analysis$name, "\n")
  cat("  ID:", analysis_id, "\n")
  cat("=====================================================\n\n")

  # Create output directory for this analysis
  output_base <- file.path(config$paths$output, analysis$output_subdir)
  output_dirs <- create_output_dirs(output_base)

  # ===== 1. DATA MERGING =====
  cat("\n--- Step 1: Data Merging ---\n")

  merged_exprs <- merge_expression_data(
    mapped_path = config$paths$mapped_data,
    pattern = config$merging$pattern
  )

  # Load phenodata
  pdata <- read.csv(config$paths$phenodata, stringsAsFactors = FALSE)

  # Align phenodata
  pdata_aligned <- align_phenodata(
    mrgd = merged_exprs,
    pdata = pdata,
    sample_col = config$merging$sample_column
  )

  # ===== 2. FILTERING =====
  cat("\n--- Step 2: Filtering (Analysis-Specific) ---\n")

  filtered <- filter_samples(
    mrgd = merged_exprs,
    pdata = pdata_aligned,
    filters = analysis$filtering
  )

  merged_exprs <- filtered$exprs
  pdata_aligned <- filtered$pdata

  cat("Filtered to:", ncol(merged_exprs), "samples\n")

  # ===== 3. BATCH CORRECTION =====
  cat("\n--- Step 3: Batch Correction ---\n")

  # Add trimester variable
  pdata_aligned <- add_trimester_variable(pdata_aligned)

  # Apply ComBat
  exprs_corrected <- apply_combat_correction(
    exprs = merged_exprs,
    pdata = pdata_aligned,
    batch_col = analysis$batch_correction$batch_column,
    mod_formula = as.formula(analysis$batch_correction$model_formula),
    par_prior = analysis$batch_correction$par_prior,
    prior_plots = analysis$batch_correction$prior_plots
  )

  # Save corrected expression
  save_data(exprs_corrected, file.path(output_base, "exprs_corrected.tsv"))
  save_data(pdata_aligned, file.path(output_base, "pdata_aligned.tsv"))

  # ===== 4. PCA =====
  cat("\n--- Step 4: PCA Analysis ---\n")

  pca <- perform_pca(exprs_corrected, center = config$pca$center, scale = config$pca$scale)

  plot_pca_metadata(
    pca = pca,
    pdata = pdata_aligned,
    variables = config$pca$color_by,
    output_dir = output_dirs$plots,
    prefix = paste0("pca_", analysis_id)
  )

  # ===== 5. DIFFERENTIAL EXPRESSION =====
  cat("\n--- Step 5: Differential Expression ---\n")

  de_results <- differential_expression(
    exprs = exprs_corrected,
    pdata = pdata_aligned,
    design_formula = as.formula(analysis$differential_expression$design_formula),
    contrasts = analysis$differential_expression$contrasts,
    method = analysis$differential_expression$decide_method
  )

  # Extract results
  difexp_all <- extract_difexp_results(
    efit = de_results$efit,
    results = de_results$results,
    exprs = exprs_corrected,
    pdata = pdata_aligned,
    species = "human"
  )

  # Filter by thresholds
  difexp_filtered <- filter_difexp(
    difexp = difexp_all,
    logfc_threshold = analysis$differential_expression$logfc_threshold,
    pvalue_threshold = analysis$differential_expression$pvalue_threshold
  )

  # Save DE results
  save_data(difexp_all, file.path(output_dirs$difexp, "difexp_all.csv"))
  save_data(difexp_filtered, file.path(output_dirs$difexp, "difexp_filtered.csv"))

  cat("Significant genes:", nrow(difexp_filtered), "\n")

  # ===== 6. NETWORK CLUSTERING =====
  cat("\n--- Step 6: Network Clustering ---\n")

  if (nrow(difexp_filtered) > 0) {
    # Map to STRING
    mapped <- map_to_stringdb(
      difexp = difexp_filtered,
      background_exprs = exprs_corrected,
      string_db = string_db,
      gene_col = "SYMBOL"
    )

    difexp_mapped <- mapped$difexp

    # Build network
    G <- build_ppi_network(difexp_mapped$STRING_id, string_db)

    # Cluster
    if (vcount(G) > 0) {
      fgreedy <- fastgreedy_clustering(G)
      difexp_clustered <- add_cluster_column(difexp_mapped, fgreedy, G)

      # Recursive clustering
      if (config$clustering$recursive) {
        cluster_recursive(
          difexp = difexp_clustered,
          string_db = string_db,
          output_dir = output_base,
          min_cluster_size = config$clustering$min_cluster_size
        )
      }
    } else {
      difexp_clustered <- difexp_mapped
      difexp_clustered$cluster <- NA
    }
  } else {
    cat("No significant genes for clustering\n")
    difexp_clustered <- difexp_filtered
    difexp_clustered$cluster <- NA
  }

  # ===== 7. ENRICHMENT =====
  cat("\n--- Step 7: Enrichment Analysis ---\n")

  if (nrow(difexp_clustered) > 0 && !all(is.na(difexp_clustered$cluster))) {
    enrichment_summary <- cluster_enrichment(
      difexp = difexp_clustered,
      string_db = string_db,
      output_dir = output_base,
      max_pvalue = config$enrichment$max_pvalue
    )
  } else {
    cat("No clusters for enrichment\n")
    enrichment_summary <- data.frame()
  }

  # Save final results
  save_data(difexp_clustered, file.path(output_base, "difexp_final.csv"))

  # Store for comparison
  all_results[[analysis_id]] <- list(
    name = analysis$name,
    difexp_all = difexp_all,
    difexp_filtered = difexp_filtered,
    difexp_final = difexp_clustered,
    n_significant = nrow(difexp_filtered),
    enrichment_summary = enrichment_summary
  )

  cat("\n✓ Analysis", analysis$name, "completed\n")
}

# ===== 8. COMPARISON ANALYSIS =====
if (config$comparison$enabled && length(all_results) >= 2) {

  cat("\n\n")
  cat("=====================================================\n")
  cat("  Comparison Analysis\n")
  cat("=====================================================\n\n")

  comparison_dir <- file.path(config$paths$output, config$comparison$output_subdir)
  if (!dir.exists(comparison_dir)) {
    dir.create(comparison_dir, recursive = TRUE)
  }

  # Get analyses to compare
  analyses_to_compare <- config$comparison$analyses_to_compare

  if (length(analyses_to_compare) != 2) {
    stop("Comparison requires exactly 2 analyses")
  }

  res1 <- all_results[[analyses_to_compare[1]]]
  res2 <- all_results[[analyses_to_compare[2]]]

  cat("Comparing:\n")
  cat("  Analysis 1:", res1$name, "(", res1$n_significant, "sig genes )\n")
  cat("  Analysis 2:", res2$name, "(", res2$n_significant, "sig genes )\n\n")

  # Get significant gene lists
  sig1 <- res1$difexp_filtered$SYMBOL
  sig2 <- res2$difexp_filtered$SYMBOL

  # Venn diagram counts
  common <- intersect(sig1, sig2)
  only1 <- setdiff(sig1, sig2)
  only2 <- setdiff(sig2, sig1)

  cat("Venn diagram:\n")
  cat("  Only in", res1$name, ":", length(only1), "\n")
  cat("  Common:", length(common), "\n")
  cat("  Only in", res2$name, ":", length(only2), "\n\n")

  # Merge results for common genes
  merged_comparison <- merge(
    res1$difexp_filtered[, c("SYMBOL", "ENTREZID", "logFC", "adj.P.Val", "AveExpr")],
    res2$difexp_filtered[, c("SYMBOL", "logFC", "adj.P.Val", "AveExpr")],
    by = "SYMBOL",
    suffixes = c(paste0("_", analyses_to_compare[1]), paste0("_", analyses_to_compare[2]))
  )

  # Add categories
  merged_comparison$category <- "common"
  merged_comparison$category[merged_comparison$SYMBOL %in% only1] <- paste0("only_", analyses_to_compare[1])
  merged_comparison$category[merged_comparison$SYMBOL %in% only2] <- paste0("only_", analyses_to_compare[2])

  # Direction comparison for common genes
  if (nrow(merged_comparison) > 0) {
    logfc1_col <- paste0("logFC_", analyses_to_compare[1])
    logfc2_col <- paste0("logFC_", analyses_to_compare[2])

    merged_comparison$same_direction <- sign(merged_comparison[[logfc1_col]]) ==
                                         sign(merged_comparison[[logfc2_col]])

    cat("Direction comparison for common genes:\n")
    cat("  Same direction:", sum(merged_comparison$same_direction & merged_comparison$category == "common"), "\n")
    cat("  Opposite direction:", sum(!merged_comparison$same_direction & merged_comparison$category == "common"), "\n\n")

    # Fold change correlation
    if (config$comparison$fold_change_correlation) {
      cor_val <- cor(merged_comparison[[logfc1_col]], merged_comparison[[logfc2_col]],
                      use = "complete.obs")
      cat("Fold change correlation:", round(cor_val, 3), "\n\n")

      # Scatter plot
      png(file.path(comparison_dir, "fold_change_scatter.png"), width = 800, height = 800)
      plot(merged_comparison[[logfc1_col]], merged_comparison[[logfc2_col]],
           xlab = paste("logFC", res1$name),
           ylab = paste("logFC", res2$name),
           main = paste("Fold Change Comparison (r =", round(cor_val, 3), ")"),
           pch = 16, col = rgb(0, 0, 0, 0.3))
      abline(h = 0, v = 0, lty = 2, col = "gray")
      abline(0, 1, col = "red", lty = 2)
      dev.off()
      cat("Saved scatter plot\n")
    }
  }

  # Save comparison results
  save_data(merged_comparison, file.path(comparison_dir, "comparison_results.csv"))

  # Summary table
  summary_df <- data.frame(
    analysis = c(res1$name, res2$name),
    n_significant = c(res1$n_significant, res2$n_significant),
    n_unique = c(length(only1), length(only2)),
    n_common = c(length(common), length(common))
  )

  save_data(summary_df, file.path(comparison_dir, "comparison_summary.csv"))

  cat("\n✓ Comparison analysis completed\n")
}

cat("\n\n=== All Analyses Complete ===\n")
cat("Results saved to:", config$paths$output, "\n\n")
