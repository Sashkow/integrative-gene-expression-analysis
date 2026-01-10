#!/usr/bin/env Rscript

#' Test which dataset gives most DE genes when added to baseline
#'
#' Baseline: Old baseline datasets (GSE122214, GSE22490, GSE37901, GSE9984)
#' Test: Add each additional dataset and compare genes and DEGs
#'
#' @author Expression Integration Pipeline
#' @date 2025-10-26

cat("\n=== Testing Additional Datasets Addition to Old Baseline ===\n\n")

# Load required libraries
library(yaml)
library(ggplot2)

# Source all modules
source("R/config.R")
source("R/utils.R")
source("R/01_data_merging.R")
source("R/02_batch_correction.R")
source("R/03_pca_analysis.R")
source("R/04_differential_expression.R")
source("R/05_network_clustering.R")

#' Filter differential expression results to only include genes in STRING DB
#'
#' @param difexp Differential expression results data frame
#' @param string_db STRINGdb object
#' @param gene_col Column name containing gene symbols (default: "SYMBOL")
#' @return Filtered difexp data frame with only STRING-mappable genes
filter_difexp_by_stringdb <- function(difexp, string_db, gene_col = "SYMBOL") {

  if (!gene_col %in% colnames(difexp)) {
    stop("Column ", gene_col, " not found in difexp data")
  }

  # Map to STRING DB
  difexp_mapped <- string_db$map(difexp, gene_col,
                                  removeUnmappedRows = TRUE,
                                  takeFirst = TRUE)

  n_before <- nrow(difexp)
  n_after <- nrow(difexp_mapped)
  n_lost <- n_before - n_after
  pct_retained <- round(100 * n_after / n_before, 1)

  cat("  STRING DB filtering: ", n_before, " -> ", n_after,
      " genes (", pct_retained, "% retained, ", n_lost, " lost)\n", sep = "")

  return(difexp_mapped)
}

#' Plot logFC-FDR comparison between baseline and current iteration
#'
#' @param baseline_difexp Data frame with baseline differential expression results
#' @param current_difexp Data frame with current iteration differential expression results
#' @param output_file Path to save the plot
#' @param title Title for the plot
plot_logfc_fdr_comparison <- function(baseline_difexp, current_difexp, output_file, title = "logFC-FDR Comparison") {

  # Merge baseline and current by gene identifier (using SYMBOL or ENTREZID)
  # Assuming both have SYMBOL column
  gene_col <- if ("SYMBOL" %in% colnames(baseline_difexp)) "SYMBOL" else "ENTREZID"

  baseline_data <- baseline_difexp[, c(gene_col, "logFC", "adj.P.Val")]
  colnames(baseline_data) <- c("gene", "logFC_baseline", "FDR_baseline")

  current_data <- current_difexp[, c(gene_col, "logFC", "adj.P.Val")]
  colnames(current_data) <- c("gene", "logFC_current", "FDR_current")

  # Merge on gene
  merged_data <- merge(baseline_data, current_data, by = "gene", all = FALSE)

  # Convert FDR to -log10(FDR) for better visualization
  merged_data$neglog10FDR_baseline <- -log10(merged_data$FDR_baseline + 1e-300)
  merged_data$neglog10FDR_current <- -log10(merged_data$FDR_current + 1e-300)

  # Cap at reasonable values for visualization
  merged_data$neglog10FDR_baseline <- pmin(merged_data$neglog10FDR_baseline, 50)
  merged_data$neglog10FDR_current <- pmin(merged_data$neglog10FDR_current, 50)

  # Determine point significance
  merged_data$significance <- "Not significant"
  merged_data$significance[merged_data$FDR_baseline < 0.05 | merged_data$FDR_current < 0.05] <- "Significant in one"
  merged_data$significance[merged_data$FDR_baseline < 0.05 & merged_data$FDR_current < 0.05] <- "Significant in both"

  # Calculate change magnitude for coloring lines
  merged_data$logFC_change <- abs(merged_data$logFC_current - merged_data$logFC_baseline)
  merged_data$FDR_change <- abs(merged_data$neglog10FDR_current - merged_data$neglog10FDR_baseline)
  merged_data$total_change <- sqrt(merged_data$logFC_change^2 + merged_data$FDR_change^2)

  # Create the plot
  p <- ggplot(merged_data) +
    # Draw lines connecting baseline and current for each gene
    geom_segment(aes(x = logFC_baseline, y = neglog10FDR_baseline,
                     xend = logFC_current, yend = neglog10FDR_current,
                     color = total_change),
                 alpha = 0.3, linewidth = 0.3) +
    # Baseline points
    geom_point(aes(x = logFC_baseline, y = neglog10FDR_baseline, shape = "Baseline"),
               color = "blue", alpha = 0.5, size = 1.5) +
    # Current iteration points
    geom_point(aes(x = logFC_current, y = neglog10FDR_current, shape = "Current"),
               color = "red", alpha = 0.5, size = 1.5) +
    # Add significance threshold lines
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50", alpha = 0.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50", alpha = 0.5) +
    # Color scale for lines
    scale_color_gradient(low = "lightgray", high = "purple", name = "Change\nMagnitude") +
    # Shape scale
    scale_shape_manual(values = c("Baseline" = 16, "Current" = 17), name = "Dataset") +
    # Labels and theme
    labs(
      title = title,
      subtitle = paste0("N genes = ", nrow(merged_data)),
      x = "logFC",
      y = "-log10(FDR)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )

  # Save plot using png device to avoid graphics API issues
  png(output_file, width = 10, height = 8, units = "in", res = 300)
  print(p)
  dev.off()

  cat("✓ Saved comparison plot to:", output_file, "\n")

  # Calculate detailed statistics
  # FDR changes
  genes_fdr_increased <- sum(merged_data$FDR_current > merged_data$FDR_baseline)
  genes_fdr_decreased <- sum(merged_data$FDR_current < merged_data$FDR_baseline)

  # Significance changes
  genes_lost_significance <- sum(merged_data$FDR_baseline < 0.05 & merged_data$FDR_current >= 0.05)
  genes_gained_significance <- sum(merged_data$FDR_baseline >= 0.05 & merged_data$FDR_current < 0.05)

  # LogFC changes
  genes_logfc_change_gt1 <- sum(merged_data$logFC_change > 1)

  # LogFC threshold crossings (considering absolute values)
  baseline_above1 <- abs(merged_data$logFC_baseline) > 1
  current_above1 <- abs(merged_data$logFC_current) > 1
  genes_logfc_fell_below1 <- sum(baseline_above1 & !current_above1)
  genes_logfc_rose_above1 <- sum(!baseline_above1 & current_above1)

  # Return summary statistics
  stats <- list(
    n_genes = nrow(merged_data),
    n_sig_baseline = sum(merged_data$FDR_baseline < 0.05),
    n_sig_current = sum(merged_data$FDR_current < 0.05),
    n_sig_both = sum(merged_data$FDR_baseline < 0.05 & merged_data$FDR_current < 0.05),
    mean_logFC_change = mean(merged_data$logFC_change),
    median_logFC_change = median(merged_data$logFC_change),
    genes_fdr_increased = genes_fdr_increased,
    genes_fdr_decreased = genes_fdr_decreased,
    genes_lost_significance = genes_lost_significance,
    genes_gained_significance = genes_gained_significance,
    genes_logfc_change_gt1 = genes_logfc_change_gt1,
    genes_logfc_fell_below1 = genes_logfc_fell_below1,
    genes_logfc_rose_above1 = genes_logfc_rose_above1
  )

  return(stats)
}

# Load configuration
config <- load_config("config/config_test_first_datasets.yaml")

# Create output directory for test iterations
test_output_dir <- "output/test_first_datasets"
dir.create(test_output_dir, showWarnings = FALSE, recursive = TRUE)
cat("Created output directory:", test_output_dir, "\n\n")

# Check if STRING DB filtering is enabled
filter_by_stringdb <- config$differential_expression$filter_by_stringdb
if (is.null(filter_by_stringdb)) {
  filter_by_stringdb <- FALSE  # Default to FALSE if not specified
}

# Initialize STRING database if filtering is enabled
string_db <- NULL
if (filter_by_stringdb) {
  cat("STRING DB filtering: ENABLED\n")
  cat("Initializing STRING database...\n")
  string_db <- initialize_stringdb(
    version = "11",
    species = 9606,
    score_threshold = 400,
    use_cache = TRUE,
    cache_dir = "data/stringdb_cache"
  )
  cat("✓ STRING database initialized\n\n")
} else {
  cat("STRING DB filtering: DISABLED\n")
  cat("All genes will be included in output\n\n")
}

# Extract datasets from config
baseline_datasets <- config$baseline$datasets
addon_datasets <- config$first_trimester_datasets$datasets

cat("Baseline datasets:", paste(baseline_datasets, collapse = ", "), "\n")
cat("Additional datasets to test:", paste(addon_datasets, collapse = ", "), "\n\n")

# Load full phenodata once
pdata_full <- read.csv(config$paths$phenodata, stringsAsFactors = FALSE)

# Define batch column
batch_col <- config$batch_correction$batch_column

# Storage for results
results_summary <- data.frame(
  dataset = character(),
  type = character(),
  n_first_trim = integer(),
  n_second_trim = integer(),
  total_samples = integer(),
  n_genes_merged = integer(),
  n_degs_all = integer(),
  n_degs_logfc1 = integer(),
  genes_fdr_increased = integer(),
  genes_fdr_decreased = integer(),
  genes_lost_significance = integer(),
  genes_gained_significance = integer(),
  genes_logfc_change_gt1 = integer(),
  genes_logfc_fell_below1 = integer(),
  genes_logfc_rose_above1 = integer(),
  stringsAsFactors = FALSE
)

# Storage for baseline comparison
baseline_genes_merged <- NULL
baseline_degs <- NULL
baseline_difexp <- NULL
baseline_difexp_all <- NULL

# ===================================================================
# STEP 1: RUN BASELINE ANALYSIS (Old baseline: GSE122214, GSE22490, GSE37901, GSE9984)
# ===================================================================

cat("\n")
cat("==================================================================\n")
cat("           BASELINE ANALYSIS                                     \n")
cat("==================================================================\n\n")

# Exclude all addon datasets - we only want the baseline datasets
baseline_exclusions <- addon_datasets

cat("Excluding:", paste(baseline_exclusions, collapse = ", "), "\n\n")

# Filter phenodata
pdata_filtered <- filter_phenodata(
  pdata = pdata_full,
  filters = config$filtering
)

pdata_filtered <- pdata_filtered[!(pdata_filtered[[batch_col]] %in% baseline_exclusions), ]

# Count samples
n_first <- sum(pdata_filtered$Gestational.Age.Category == "First Trimester", na.rm = TRUE)
n_second <- sum(pdata_filtered$Gestational.Age.Category == "Second Trimester", na.rm = TRUE)

cat("Samples: First Trimester =", n_first, ", Second Trimester =", n_second, "\n")

# Identify required datasets
required_datasets <- identify_datasets(
  pdata_filtered = pdata_filtered,
  batch_col = batch_col
)

cat("Datasets included:", paste(required_datasets, collapse = ", "), "\n")

# Run baseline analysis
tryCatch({
  # Merge expression data
  merged_exprs <- merge_selected_datasets(
    mapped_path = config$paths$mapped_data,
    datasets = required_datasets,
    pattern = config$merging$pattern
  )

  # Align phenodata
  pdata_aligned <- pdata_filtered[
    make.names(pdata_filtered[[config$merging$sample_column]]) %in% colnames(merged_exprs),
  ]

  # Subset expression
  merged_exprs <- subset_expression_by_phenodata(
    mrgd = merged_exprs,
    pdata_filtered = pdata_aligned,
    sample_col = config$merging$sample_column
  )

  cat("Expression matrix:", nrow(merged_exprs), "genes x", ncol(merged_exprs), "samples\n")

  # Save baseline genes
  baseline_genes_merged <- rownames(merged_exprs)

  # Add trimester variable
  pdata_aligned <- add_trimester_variable(pdata_aligned)

  # Batch correction
  cat("Applying ComBat batch correction...\n")
  exprs_corrected <- apply_combat_correction(
    exprs = merged_exprs,
    pdata = pdata_aligned,
    batch_col = batch_col,
    mod_formula = as.formula(config$batch_correction$model_formula),
    par_prior = config$batch_correction$par_prior,
    prior_plots = FALSE
  )

  # Differential expression
  cat("Running differential expression...\n")
  de_results <- differential_expression(
    exprs = exprs_corrected,
    pdata = pdata_aligned,
    design_formula = as.formula(config$differential_expression$design_formula),
    contrasts = config$differential_expression$contrasts,
    method = config$differential_expression$decide_method
  )

  # Extract results
  difexp_all <- extract_difexp_results(
    efit = de_results$efit,
    results = de_results$results,
    exprs = exprs_corrected,
    pdata = pdata_aligned,
    species = "human"
  )

  # Count DEGs
  difexp_all_sig <- filter_difexp(
    difexp = difexp_all,
    logfc_threshold = 0,
    pvalue_threshold = config$differential_expression$pvalue_threshold
  )

  difexp_logfc1 <- filter_difexp(
    difexp = difexp_all,
    logfc_threshold = 1,
    pvalue_threshold = config$differential_expression$pvalue_threshold
  )

  n_degs_all <- nrow(difexp_all_sig)
  n_degs_logfc1 <- nrow(difexp_logfc1)

  cat("✓ DEGs (p < 0.05):", n_degs_all, "\n")
  cat("✓ DEGs (p < 0.05, |logFC| > 1):", n_degs_logfc1, "\n")

  # Conditionally filter by STRING DB
  if (filter_by_stringdb) {
    cat("Filtering for STRING DB proteins...\n")
    difexp_logfc1_final <- filter_difexp_by_stringdb(difexp_logfc1, string_db, gene_col = "SYMBOL")
    difexp_all_final <- filter_difexp_by_stringdb(difexp_all, string_db, gene_col = "SYMBOL")

    n_degs_logfc1_final <- nrow(difexp_logfc1_final)
    cat("✓ DEGs with STRING proteins (|logFC| > 1):", n_degs_logfc1_final, "\n")
  } else {
    # No STRING filtering - use original results
    difexp_logfc1_final <- difexp_logfc1
    difexp_all_final <- difexp_all
    n_degs_logfc1_final <- n_degs_logfc1
  }

  # Save baseline DEGs and unfiltered results for comparison
  baseline_degs <- difexp_logfc1_final$ENTREZID
  baseline_difexp <- difexp_logfc1_final
  baseline_difexp_all <- difexp_all_final

  # Store baseline results (no comparison stats for baseline)
  baseline_name <- paste0("BASELINE (", paste(baseline_datasets, collapse = ", "), ")")
  results_summary <- rbind(results_summary, data.frame(
    dataset = baseline_name,
    type = "baseline",
    n_first_trim = n_first,
    n_second_trim = n_second,
    total_samples = ncol(merged_exprs),
    n_genes_merged = nrow(merged_exprs),
    n_degs_all = n_degs_all,
    n_degs_logfc1 = n_degs_logfc1,
    genes_fdr_increased = NA,
    genes_fdr_decreased = NA,
    genes_lost_significance = NA,
    genes_gained_significance = NA,
    genes_logfc_change_gt1 = NA,
    genes_logfc_fell_below1 = NA,
    genes_logfc_rose_above1 = NA,
    stringsAsFactors = FALSE
  ))

  # Save baseline iteration files
  baseline_dir <- file.path(test_output_dir, "BASELINE")
  dir.create(baseline_dir, showWarnings = FALSE, recursive = TRUE)

  # Save expression matrix (before ComBat)
  write.table(merged_exprs,
              file.path(baseline_dir, "merged_exprs_before_combat.tsv"),
              sep = "\t", quote = FALSE)

  # Save expression matrix (after ComBat)
  write.table(exprs_corrected,
              file.path(baseline_dir, "merged_exprs_after_combat.tsv"),
              sep = "\t", quote = FALSE)

  # Save phenodata
  write.csv(pdata_aligned,
            file.path(baseline_dir, "phenodata.csv"),
            row.names = FALSE)

  # Save differential expression results (filtered by logFC)
  write.csv(difexp_logfc1_final,
            file.path(baseline_dir, "difexp_genes.csv"),
            row.names = FALSE)

  # Save unfiltered differential expression results (all genes)
  write.csv(difexp_all_final,
            file.path(baseline_dir, "difexp_genes_unfiltered.csv"),
            row.names = FALSE)

  # Save original results before STRING filtering (only if filtering was applied)
  if (filter_by_stringdb) {
    write.csv(difexp_logfc1,
              file.path(baseline_dir, "difexp_genes_before_stringdb.csv"),
              row.names = FALSE)

    write.csv(difexp_all,
              file.path(baseline_dir, "difexp_genes_unfiltered_before_stringdb.csv"),
              row.names = FALSE)
  }

  cat("✓ Saved baseline files to:", baseline_dir, "\n")
  cat("✓ Baseline analysis complete\n")

}, error = function(e) {
  cat("✗ Error in baseline analysis:", conditionMessage(e), "\n")
})

# ===================================================================
# STEP 2: TEST ADDITIONAL DATASETS
# ===================================================================

cat("\n\n")
cat("==================================================================\n")
cat("           TESTING ADDITIONAL DATASETS                           \n")
cat("==================================================================\n\n")

# Determine which datasets to test based on config
test_all_subsets <- config$testing$test_all_subsets

if (test_all_subsets) {
  # Generate all non-empty subsets of addon_datasets
  n <- length(addon_datasets)
  all_combinations <- list()

  for (i in 1:(2^n - 1)) {
    subset_idx <- as.logical(intToBits(i)[1:n])
    combination <- addon_datasets[subset_idx]
    all_combinations[[length(all_combinations) + 1]] <- combination
  }

  cat("Testing mode: ALL SUBSETS\n")
  cat("Total combinations to test:", length(all_combinations), "\n\n")
  datasets_to_test <- all_combinations
} else {
  cat("Testing mode: INDIVIDUAL DATASETS\n")
  cat("Datasets to test:", paste(addon_datasets, collapse = ", "), "\n\n")
  datasets_to_test <- as.list(addon_datasets)
}

for (test_datasets in datasets_to_test) {

  # Convert to character vector if single dataset
  if (length(test_datasets) == 1) {
    test_name <- test_datasets
    first_ds <- test_datasets
  } else {
    test_name <- paste(test_datasets, collapse = "+")
    first_ds <- test_datasets
  }

  cat("\n")
  cat("========================================\n")
  cat("Testing with additional dataset(s):", test_name, "\n")
  cat("========================================\n\n")

  # Create exclusion list (exclude all addon datasets EXCEPT the ones we're testing)
  other_addon_datasets <- setdiff(addon_datasets, test_datasets)
  exclusions <- other_addon_datasets

  cat("Excluding:", paste(exclusions, collapse = ", "), "\n\n")

  # Filter phenodata
  pdata_filtered <- filter_phenodata(
    pdata = pdata_full,
    filters = config$filtering
  )

  # Apply exclusions
  pdata_filtered <- pdata_filtered[!(pdata_filtered[[batch_col]] %in% exclusions), ]

  # Count samples
  n_first <- sum(pdata_filtered$Gestational.Age.Category == "First Trimester", na.rm = TRUE)
  n_second <- sum(pdata_filtered$Gestational.Age.Category == "Second Trimester", na.rm = TRUE)

  cat("Samples: First Trimester =", n_first, ", Second Trimester =", n_second, "\n")

  # Skip if no First samples
  if (n_first == 0) {
    cat("⚠ No First Trimester samples for", test_name, "- skipping\n")
    next
  }

  # Identify required datasets
  required_datasets <- identify_datasets(
    pdata_filtered = pdata_filtered,
    batch_col = batch_col
  )

  cat("Datasets included:", paste(required_datasets, collapse = ", "), "\n")

  # Run analysis
  tryCatch({
    # Merge expression data
    merged_exprs <- merge_selected_datasets(
      mapped_path = config$paths$mapped_data,
      datasets = required_datasets,
      pattern = config$merging$pattern
    )

    # Align phenodata
    pdata_aligned <- pdata_filtered[
      make.names(pdata_filtered[[config$merging$sample_column]]) %in% colnames(merged_exprs),
    ]

    # Subset expression
    merged_exprs <- subset_expression_by_phenodata(
      mrgd = merged_exprs,
      pdata_filtered = pdata_aligned,
      sample_col = config$merging$sample_column
    )

    cat("Expression matrix:", nrow(merged_exprs), "genes x", ncol(merged_exprs), "samples\n")

    # Compare genes to baseline
    if (!is.null(baseline_genes_merged)) {
      current_genes <- rownames(merged_exprs)
      genes_gained_vec <- setdiff(current_genes, baseline_genes_merged)
      genes_lost_vec <- setdiff(baseline_genes_merged, current_genes)
      genes_common <- length(intersect(current_genes, baseline_genes_merged))

      cat("Genes comparison to baseline:\n")
      cat("  Common genes:", genes_common, "\n")
      cat("  Genes gained:", length(genes_gained_vec), "\n")
      cat("  Genes lost:", length(genes_lost_vec), "\n")

      # Print gained genes
      if (length(genes_gained_vec) > 0) {
        cat("\n  Gained genes (symbols):\n")
        cat("  ", paste(head(genes_gained_vec, 50), collapse = ", "), "\n")
        if (length(genes_gained_vec) > 50) {
          cat("  ... and", length(genes_gained_vec) - 50, "more\n")
        }
        # Save to file
        gained_file <- file.path(test_output_dir, paste0("first_", gsub("\\+", "_", test_name), "_genes_gained.txt"))
        writeLines(genes_gained_vec, gained_file)
        cat("  Saved to:", gained_file, "\n")
      }

      # Print lost genes
      if (length(genes_lost_vec) > 0) {
        cat("\n  Lost genes (symbols):\n")
        cat("  ", paste(head(genes_lost_vec, 50), collapse = ", "), "\n")
        if (length(genes_lost_vec) > 50) {
          cat("  ... and", length(genes_lost_vec) - 50, "more\n")
        }
        # Save to file
        lost_file <- file.path(test_output_dir, paste0("first_", gsub("\\+", "_", test_name), "_genes_lost.txt"))
        writeLines(genes_lost_vec, lost_file)
        cat("  Saved to:", lost_file, "\n")
      }
      cat("\n")
    }

    # Add trimester variable
    pdata_aligned <- add_trimester_variable(pdata_aligned)

    # Batch correction
    cat("Applying ComBat batch correction...\n")
    exprs_corrected <- apply_combat_correction(
      exprs = merged_exprs,
      pdata = pdata_aligned,
      batch_col = batch_col,
      mod_formula = as.formula(config$batch_correction$model_formula),
      par_prior = config$batch_correction$par_prior,
      prior_plots = FALSE
    )

    # Differential expression
    cat("Running differential expression...\n")
    de_results <- differential_expression(
      exprs = exprs_corrected,
      pdata = pdata_aligned,
      design_formula = as.formula(config$differential_expression$design_formula),
      contrasts = config$differential_expression$contrasts,
      method = config$differential_expression$decide_method
    )

    # Extract results
    difexp_all <- extract_difexp_results(
      efit = de_results$efit,
      results = de_results$results,
      exprs = exprs_corrected,
      pdata = pdata_aligned,
      species = "human"
    )

    # Count DEGs
    difexp_all_sig <- filter_difexp(
      difexp = difexp_all,
      logfc_threshold = 0,
      pvalue_threshold = config$differential_expression$pvalue_threshold
    )

    difexp_logfc1 <- filter_difexp(
      difexp = difexp_all,
      logfc_threshold = 1,
      pvalue_threshold = config$differential_expression$pvalue_threshold
    )

    n_degs_all <- nrow(difexp_all_sig)
    n_degs_logfc1 <- nrow(difexp_logfc1)

    cat("✓ DEGs (p < 0.05):", n_degs_all, "\n")
    cat("✓ DEGs (p < 0.05, |logFC| > 1):", n_degs_logfc1, "\n")

    # Conditionally filter by STRING DB
    if (filter_by_stringdb) {
      cat("Filtering for STRING DB proteins...\n")
      difexp_logfc1_final <- filter_difexp_by_stringdb(difexp_logfc1, string_db, gene_col = "SYMBOL")
      difexp_all_final <- filter_difexp_by_stringdb(difexp_all, string_db, gene_col = "SYMBOL")

      n_degs_logfc1_final <- nrow(difexp_logfc1_final)
      cat("✓ DEGs with STRING proteins (|logFC| > 1):", n_degs_logfc1_final, "\n")
    } else {
      # No STRING filtering - use original results
      difexp_logfc1_final <- difexp_logfc1
      difexp_all_final <- difexp_all
      n_degs_logfc1_final <- n_degs_logfc1
    }

    # Compare DEGs to baseline
    if (!is.null(baseline_degs)) {
      current_degs <- difexp_logfc1_final$ENTREZID
      degs_common <- length(intersect(current_degs, baseline_degs))
      degs_gained <- length(setdiff(current_degs, baseline_degs))
      degs_lost <- length(setdiff(baseline_degs, current_degs))

      cat("DEGs comparison to baseline:\n")
      cat("  Common DEGs:", degs_common, "\n")
      cat("  DEGs gained:", degs_gained, "\n")
      cat("  DEGs lost:", degs_lost, "\n")
      cat("  Net change:", degs_gained - degs_lost, "\n")
    }

    # Calculate detailed comparison statistics with baseline
    comparison_stats <- NULL
    if (!is.null(baseline_difexp_all)) {
      # Merge baseline and current by gene
      gene_col <- if ("SYMBOL" %in% colnames(baseline_difexp_all)) "SYMBOL" else "ENTREZID"
      baseline_data <- baseline_difexp_all[, c(gene_col, "logFC", "adj.P.Val", "P.Value")]
      colnames(baseline_data) <- c("gene", "logFC_baseline", "FDR_baseline", "P_baseline")
      current_data <- difexp_all_final[, c(gene_col, "logFC", "adj.P.Val", "P.Value")]
      colnames(current_data) <- c("gene", "logFC_current", "FDR_current", "P_current")
      merged_data <- merge(baseline_data, current_data, by = "gene", all = FALSE)

      # Check if FDR normalization is enabled
      normalize_fdr <- ifelse(is.null(config$differential_expression$normalize_fdr_to_baseline),
                              FALSE,
                              config$differential_expression$normalize_fdr_to_baseline)

      cat("FDR normalization setting: normalize_fdr =", normalize_fdr, "\n")
      cat("baseline_genes_merged is null:", is.null(baseline_genes_merged), "\n")

      # Recalculate FDR for current data as if it had baseline gene count
      # This makes FDR comparable across iterations with different gene counts
      if (normalize_fdr && !is.null(baseline_genes_merged)) {
        baseline_n_genes <- length(baseline_genes_merged)
        current_n_genes <- nrow(difexp_all_final)

        cat("Recalculating FDR with baseline gene count for fair comparison...\n")
        cat("  Baseline genes:", baseline_n_genes, "\n")
        cat("  Current genes:", current_n_genes, "\n")

        # Adjust FDR: multiply p-values by ratio of gene counts, then apply BH correction
        # If current has fewer genes, adjusted FDR will be MORE stringent (higher values)
        merged_data$FDR_current_adjusted <- p.adjust(merged_data$P_current, method = "BH", n = baseline_n_genes)

        cat("  Using adjusted FDR for comparisons (normalized to", baseline_n_genes, "genes)\n")
      } else {
        # No adjustment - use original FDR
        merged_data$FDR_current_adjusted <- merged_data$FDR_current
      }

      # Calculate statistics using adjusted FDR
      comparison_stats <- list(
        genes_fdr_increased = sum(merged_data$FDR_current_adjusted > merged_data$FDR_baseline),
        genes_fdr_decreased = sum(merged_data$FDR_current_adjusted < merged_data$FDR_baseline),
        genes_lost_significance = sum(merged_data$FDR_baseline < 0.05 & merged_data$FDR_current_adjusted >= 0.05),
        genes_gained_significance = sum(merged_data$FDR_baseline >= 0.05 & merged_data$FDR_current_adjusted < 0.05),
        genes_logfc_change_gt1 = sum(abs(merged_data$logFC_current - merged_data$logFC_baseline) > 1),
        genes_logfc_fell_below1 = sum(abs(merged_data$logFC_baseline) > 1 & abs(merged_data$logFC_current) <= 1),
        genes_logfc_rose_above1 = sum(abs(merged_data$logFC_baseline) <= 1 & abs(merged_data$logFC_current) > 1)
      )
    }

    # Store results
    results_summary <- rbind(results_summary, data.frame(
      dataset = test_name,
      type = "with_first",
      n_first_trim = n_first,
      n_second_trim = n_second,
      total_samples = ncol(merged_exprs),
      n_genes_merged = nrow(merged_exprs),
      n_degs_all = n_degs_all,
      n_degs_logfc1 = n_degs_logfc1,
      genes_fdr_increased = if(!is.null(comparison_stats)) comparison_stats$genes_fdr_increased else NA,
      genes_fdr_decreased = if(!is.null(comparison_stats)) comparison_stats$genes_fdr_decreased else NA,
      genes_lost_significance = if(!is.null(comparison_stats)) comparison_stats$genes_lost_significance else NA,
      genes_gained_significance = if(!is.null(comparison_stats)) comparison_stats$genes_gained_significance else NA,
      genes_logfc_change_gt1 = if(!is.null(comparison_stats)) comparison_stats$genes_logfc_change_gt1 else NA,
      genes_logfc_fell_below1 = if(!is.null(comparison_stats)) comparison_stats$genes_logfc_fell_below1 else NA,
      genes_logfc_rose_above1 = if(!is.null(comparison_stats)) comparison_stats$genes_logfc_rose_above1 else NA,
      stringsAsFactors = FALSE
    ))

    # Save iteration files
    iteration_dir <- file.path(test_output_dir, gsub("\\+", "_", test_name))
    dir.create(iteration_dir, showWarnings = FALSE, recursive = TRUE)

    # Save expression matrix (before ComBat)
    write.table(merged_exprs,
                file.path(iteration_dir, "merged_exprs_before_combat.tsv"),
                sep = "\t", quote = FALSE)

    # Save expression matrix (after ComBat)
    write.table(exprs_corrected,
                file.path(iteration_dir, "merged_exprs_after_combat.tsv"),
                sep = "\t", quote = FALSE)

    # Save phenodata
    write.csv(pdata_aligned,
              file.path(iteration_dir, "phenodata.csv"),
              row.names = FALSE)

    # Save differential expression results (filtered by logFC)
    write.csv(difexp_logfc1_final,
              file.path(iteration_dir, "difexp_genes.csv"),
              row.names = FALSE)

    # Save unfiltered differential expression results (all genes)
    write.csv(difexp_all_final,
              file.path(iteration_dir, "difexp_genes_unfiltered.csv"),
              row.names = FALSE)

    # Save original results before STRING filtering (only if filtering was applied)
    if (filter_by_stringdb) {
      write.csv(difexp_logfc1,
                file.path(iteration_dir, "difexp_genes_before_stringdb.csv"),
                row.names = FALSE)

      write.csv(difexp_all,
                file.path(iteration_dir, "difexp_genes_unfiltered_before_stringdb.csv"),
                row.names = FALSE)
    }

    # Create logFC-FDR comparison plot with baseline
    if (!is.null(baseline_difexp_all)) {
      cat("\nGenerating logFC-FDR comparison plot...\n")
      plot_file <- file.path(iteration_dir, "logFC_FDR_comparison.png")
      plot_title <- paste0("Baseline vs ", test_name)

      tryCatch({
        plot_stats <- plot_logfc_fdr_comparison(
          baseline_difexp = baseline_difexp_all,
          current_difexp = difexp_all_final,
          output_file = plot_file,
          title = plot_title
        )

        cat("  Common genes:", plot_stats$n_genes, "\n")
        cat("  Significant in baseline:", plot_stats$n_sig_baseline, "\n")
        cat("  Significant in current:", plot_stats$n_sig_current, "\n")
        cat("  Significant in both:", plot_stats$n_sig_both, "\n")
        cat("  Mean |logFC| change:", round(plot_stats$mean_logFC_change, 3), "\n")

      }, error = function(e) {
        cat("⚠ Could not generate comparison plot:", conditionMessage(e), "\n")
      })
    }

    cat("✓ Saved iteration files to:", iteration_dir, "\n")
    cat("✓ Analysis complete\n")

  }, error = function(e) {
    cat("✗ Error with", test_name, ":", conditionMessage(e), "\n")
    cat("Error details:", as.character(e), "\n")
  })
}

# ===================================================================
# STEP 3: SINGLE DATASET ANALYSES
# ===================================================================

cat("\n\n")
cat("==================================================================\n")
cat("           SINGLE DATASET ANALYSES                               \n")
cat("==================================================================\n\n")

# Identify datasets with both First Trimester and Second Trimester samples
pdata_for_singles <- filter_phenodata(
  pdata = pdata_full,
  filters = config$filtering
)

# Count samples per dataset
dataset_counts <- table(pdata_for_singles[[batch_col]], pdata_for_singles$Gestational.Age.Category)

# Find datasets with both trimesters (at least 3 samples each for power)
datasets_with_both <- rownames(dataset_counts)[
  dataset_counts[, "First Trimester"] >= 3 & dataset_counts[, "Second Trimester"] >= 3
]

cat("Datasets with both First Trimester and Second Trimester (>=3 each):\n")
cat(paste(datasets_with_both, collapse = ", "), "\n\n")

# Storage for single dataset results
single_results <- data.frame(
  dataset = character(),
  n_first_trim = integer(),
  n_second_trim = integer(),
  total_samples = integer(),
  n_genes = integer(),
  n_degs_all = integer(),
  n_degs_logfc1 = integer(),
  stringsAsFactors = FALSE
)

# Analyze each dataset individually
for (single_ds in datasets_with_both) {

  cat("\n")
  cat("========================================\n")
  cat("Analyzing dataset:", single_ds, "\n")
  cat("========================================\n\n")

  # Filter to only this dataset
  pdata_single <- pdata_for_singles[pdata_for_singles[[batch_col]] == single_ds, ]

  # Count samples
  n_first <- sum(pdata_single$Gestational.Age.Category == "First Trimester", na.rm = TRUE)
  n_second <- sum(pdata_single$Gestational.Age.Category == "Second Trimester", na.rm = TRUE)

  cat("Samples: First Trimester =", n_first, ", Second Trimester =", n_second, "\n")

  tryCatch({
    # Load expression data for this dataset only
    merged_exprs <- merge_selected_datasets(
      mapped_path = config$paths$mapped_data,
      datasets = single_ds,
      pattern = config$merging$pattern
    )

    # Align phenodata
    pdata_aligned <- pdata_single[
      make.names(pdata_single[[config$merging$sample_column]]) %in% colnames(merged_exprs),
    ]

    # Subset expression
    merged_exprs <- subset_expression_by_phenodata(
      mrgd = merged_exprs,
      pdata_filtered = pdata_aligned,
      sample_col = config$merging$sample_column
    )

    cat("Expression matrix:", nrow(merged_exprs), "genes x", ncol(merged_exprs), "samples\n")

    # Add trimester variable
    pdata_aligned <- add_trimester_variable(pdata_aligned)

    # NO batch correction for single dataset (only one batch)
    cat("Single dataset - no batch correction needed\n")
    exprs_corrected <- merged_exprs

    # Differential expression
    cat("Running differential expression...\n")
    de_results <- differential_expression(
      exprs = exprs_corrected,
      pdata = pdata_aligned,
      design_formula = as.formula(config$differential_expression$design_formula),
      contrasts = config$differential_expression$contrasts,
      method = config$differential_expression$decide_method
    )

    # Extract results
    difexp_all <- extract_difexp_results(
      efit = de_results$efit,
      results = de_results$results,
      exprs = exprs_corrected,
      pdata = pdata_aligned,
      species = "human"
    )

    # Count DEGs
    difexp_all_sig <- filter_difexp(
      difexp = difexp_all,
      logfc_threshold = 0,
      pvalue_threshold = config$differential_expression$pvalue_threshold
    )

    difexp_logfc1 <- filter_difexp(
      difexp = difexp_all,
      logfc_threshold = 1,
      pvalue_threshold = config$differential_expression$pvalue_threshold
    )

    n_degs_all <- nrow(difexp_all_sig)
    n_degs_logfc1 <- nrow(difexp_logfc1)

    cat("✓ DEGs (p < 0.05):", n_degs_all, "\n")
    cat("✓ DEGs (p < 0.05, |logFC| > 1):", n_degs_logfc1, "\n")

    # Save differential expression results (filtered)
    difexp_file <- file.path(test_output_dir, paste0("first_", single_ds, "_single_difexp.csv"))
    write.csv(difexp_logfc1, difexp_file, row.names = FALSE)
    cat("Saved DEG results to:", difexp_file, "\n")

    # Save unfiltered differential expression results (all genes)
    difexp_file_unfiltered <- file.path(test_output_dir, paste0("first_", single_ds, "_single_difexp_unfiltered.csv"))
    write.csv(difexp_all, difexp_file_unfiltered, row.names = FALSE)
    cat("Saved unfiltered results to:", difexp_file_unfiltered, "\n")

    # Store results
    single_results <- rbind(single_results, data.frame(
      dataset = single_ds,
      n_first_trim = n_first,
      n_second_trim = n_second,
      total_samples = ncol(merged_exprs),
      n_genes = nrow(merged_exprs),
      n_degs_all = n_degs_all,
      n_degs_logfc1 = n_degs_logfc1,
      stringsAsFactors = FALSE
    ))

    cat("✓ Single dataset analysis complete\n")

  }, error = function(e) {
    cat("✗ Error with", single_ds, ":", conditionMessage(e), "\n")
  })
}

# Summary of single dataset analyses
if (nrow(single_results) > 0) {
  cat("\n\n")
  cat("==================================================================\n")
  cat("           SINGLE DATASET SUMMARY                                \n")
  cat("==================================================================\n\n")

  print(single_results, row.names = FALSE)

  # Save summary
  single_summary_file <- file.path(test_output_dir, "first_single_dataset_summary.csv")
  write.csv(single_results, single_summary_file, row.names = FALSE)
  cat("\nSingle dataset summary saved to:", single_summary_file, "\n")
}

# ===================================================================
# STEP 4: SUMMARY AND DETAILED REPORT
# ===================================================================

cat("\n\n")
cat("==================================================================\n")
cat("                        SUMMARY RESULTS                           \n")
cat("==================================================================\n\n")

# Print dataset composition summary
cat("BASELINE DATASETS:\n")
pdata_baseline_info <- filter_phenodata(
  pdata = pdata_full,
  filters = config$filtering
)
pdata_baseline_info <- pdata_baseline_info[!(pdata_baseline_info[[batch_col]] %in% baseline_exclusions), ]
baseline_table <- table(pdata_baseline_info[[batch_col]], pdata_baseline_info$Gestational.Age.Category)
print(baseline_table)
cat("\nBaseline total: First Trimester =", sum(baseline_table[, "First Trimester"]),
    ", Second Trimester =", sum(baseline_table[, "Second Trimester"]), "\n\n")

cat("ADDITIONAL DATASETS TESTED:\n")
cat(paste(addon_datasets, collapse = ", "), "\n\n")

cat("ALL DATASETS BY TRIMESTER:\n")
pdata_all_info <- filter_phenodata(
  pdata = pdata_full,
  filters = config$filtering
)
all_table <- table(pdata_all_info[[batch_col]], pdata_all_info$Gestational.Age.Category)
print(all_table)
cat("\nTotal: First Trimester =", sum(all_table[, "First Trimester"]),
    ", Second Trimester =", sum(all_table[, "Second Trimester"]), "\n\n")

cat("------------------------------------------------------------------\n")
cat("                    COMPARISON RESULTS                            \n")
cat("------------------------------------------------------------------\n\n")

if (nrow(results_summary) > 0) {
  print(results_summary, row.names = FALSE)

  # Save summary
  comparison_summary_file <- file.path(test_output_dir, "first_dataset_comparison_summary.csv")
  write.csv(results_summary, comparison_summary_file, row.names = FALSE)
  cat("\nSummary saved to:", comparison_summary_file, "\n\n")

  # Find best dataset combination
  with_addons <- results_summary[results_summary$type == "with_first", ]
  if (nrow(with_addons) > 0) {
    with_addons <- with_addons[order(-with_addons$n_degs_logfc1), ]
    cat("\n==================================================================\n")
    cat("BEST ADDITIONAL DATASET(S):", with_addons$dataset[1], "\n")
    cat("DEGs (|logFC| > 1):", with_addons$n_degs_logfc1[1], "\n")
    cat("==================================================================\n\n")
  }

} else {
  cat("⚠ No results\n\n")
}

cat("\n")
cat("==================================================================\n")
cat("                    ANALYSIS COMPLETE                             \n")
cat("==================================================================\n\n")
