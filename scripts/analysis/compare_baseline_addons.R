#!/usr/bin/env Rscript

#' Core functions for comparing baseline vs baseline+addon dataset combinations
#'
#' Shared logic for test_first_datasets.R and test_term_datasets.R
#'
#' @author Expression Integration Pipeline

library(yaml)
library(ggplot2)

#' Filter differential expression results to only include genes in STRING DB
#'
#' @param difexp Differential expression results data frame
#' @param string_db STRINGdb object
#' @param gene_col Column name containing gene symbols (default: "SYMBOL")
#' @param log_excluded If TRUE, returns list with filtered data and excluded genes
#' @return Filtered difexp data frame with only STRING-mappable genes
filter_difexp_by_stringdb <- function(difexp, string_db, gene_col = "SYMBOL",
                                       log_excluded = FALSE) {
  if (!gene_col %in% colnames(difexp)) {
    stop("Column ", gene_col, " not found in difexp data")
  }

  difexp_with_mapping <- string_db$map(difexp, gene_col,
                                        removeUnmappedRows = FALSE,
                                        takeFirst = TRUE)

  unmapped_mask <- is.na(difexp_with_mapping$STRING_id)
  excluded_genes <- difexp_with_mapping[[gene_col]][unmapped_mask]
  difexp_mapped <- difexp_with_mapping[!unmapped_mask, ]

  n_before <- nrow(difexp)
  n_after <- nrow(difexp_mapped)
  n_lost <- n_before - n_after
  pct_retained <- round(100 * n_after / n_before, 1)

  cat("  STRING DB filtering: ", n_before, " -> ", n_after,
      " genes (", pct_retained, "% retained, ", n_lost, " lost)\n", sep = "")

  if (log_excluded && n_lost > 0) {
    cat("  Excluded genes: ", paste(head(excluded_genes, 20), collapse = ", "),
        if (n_lost > 20) paste0(" ... and ", n_lost - 20, " more") else "",
        "\n", sep = "")
  }

  if (log_excluded) {
    return(list(filtered = difexp_mapped, excluded = excluded_genes))
  }
  return(difexp_mapped)
}

#' Plot logFC-FDR comparison between baseline and current iteration
#'
#' @param baseline_difexp Data frame with baseline differential expression results
#' @param current_difexp Data frame with current iteration differential expression results
#' @param output_file Path to save the plot
#' @param title Title for the plot
#' @return List of comparison statistics
plot_logfc_fdr_comparison <- function(baseline_difexp, current_difexp, output_file,
                                       title = "logFC-FDR Comparison") {
  gene_col <- if ("SYMBOL" %in% colnames(baseline_difexp)) "SYMBOL" else "ENTREZID"

  baseline_data <- baseline_difexp[, c(gene_col, "logFC", "adj.P.Val")]
  colnames(baseline_data) <- c("gene", "logFC_baseline", "FDR_baseline")

  current_data <- current_difexp[, c(gene_col, "logFC", "adj.P.Val")]
  colnames(current_data) <- c("gene", "logFC_current", "FDR_current")

  merged_data <- merge(baseline_data, current_data, by = "gene", all = FALSE)

  merged_data$neglog10FDR_baseline <- -log10(merged_data$FDR_baseline + 1e-300)
  merged_data$neglog10FDR_current <- -log10(merged_data$FDR_current + 1e-300)
  merged_data$neglog10FDR_baseline <- pmin(merged_data$neglog10FDR_baseline, 50)
  merged_data$neglog10FDR_current <- pmin(merged_data$neglog10FDR_current, 50)

  merged_data$significance <- "Not significant"
  merged_data$significance[merged_data$FDR_baseline < 0.05 | merged_data$FDR_current < 0.05] <- "Significant in one"
  merged_data$significance[merged_data$FDR_baseline < 0.05 & merged_data$FDR_current < 0.05] <- "Significant in both"

  merged_data$logFC_change <- abs(merged_data$logFC_current - merged_data$logFC_baseline)
  merged_data$FDR_change <- abs(merged_data$neglog10FDR_current - merged_data$neglog10FDR_baseline)
  merged_data$total_change <- sqrt(merged_data$logFC_change^2 + merged_data$FDR_change^2)

  p <- ggplot(merged_data) +
    geom_segment(aes(x = logFC_baseline, y = neglog10FDR_baseline,
                     xend = logFC_current, yend = neglog10FDR_current,
                     color = total_change),
                 alpha = 0.3, linewidth = 0.3) +
    geom_point(aes(x = logFC_baseline, y = neglog10FDR_baseline, shape = "Baseline"),
               color = "blue", alpha = 0.5, size = 1.5) +
    geom_point(aes(x = logFC_current, y = neglog10FDR_current, shape = "Current"),
               color = "red", alpha = 0.5, size = 1.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50", alpha = 0.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50", alpha = 0.5) +
    scale_color_gradient(low = "lightgray", high = "purple", name = "Change\nMagnitude") +
    scale_shape_manual(values = c("Baseline" = 16, "Current" = 17), name = "Dataset") +
    labs(title = title, subtitle = paste0("N genes = ", nrow(merged_data)),
         x = "logFC", y = "-log10(FDR)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          legend.position = "right", panel.grid.minor = element_blank())

  png(output_file, width = 10, height = 8, units = "in", res = 300)
  print(p)
  dev.off()
  cat("Saved comparison plot to:", output_file, "\n")

  baseline_above1 <- abs(merged_data$logFC_baseline) > 1
  current_above1 <- abs(merged_data$logFC_current) > 1

  list(
    n_genes = nrow(merged_data),
    n_sig_baseline = sum(merged_data$FDR_baseline < 0.05),
    n_sig_current = sum(merged_data$FDR_current < 0.05),
    n_sig_both = sum(merged_data$FDR_baseline < 0.05 & merged_data$FDR_current < 0.05),
    mean_logFC_change = mean(merged_data$logFC_change),
    median_logFC_change = median(merged_data$logFC_change),
    genes_fdr_increased = sum(merged_data$FDR_current > merged_data$FDR_baseline),
    genes_fdr_decreased = sum(merged_data$FDR_current < merged_data$FDR_baseline),
    genes_lost_significance = sum(merged_data$FDR_baseline < 0.05 & merged_data$FDR_current >= 0.05),
    genes_gained_significance = sum(merged_data$FDR_baseline >= 0.05 & merged_data$FDR_current < 0.05),
    genes_logfc_change_gt1 = sum(merged_data$logFC_change > 1),
    genes_logfc_fell_below1 = sum(baseline_above1 & !current_above1),
    genes_logfc_rose_above1 = sum(!baseline_above1 & current_above1)
  )
}

# Source shared logging utilities
source("scripts/utils/logging_utils.R")

#' Initialize STRING database if filtering is enabled
initialize_string_db_if_needed <- function(filter_by_stringdb) {
  if (is.null(filter_by_stringdb)) filter_by_stringdb <- FALSE

  if (filter_by_stringdb) {
    cat("STRING DB filtering: ENABLED\n")
    cat("Initializing STRING database...\n")
    string_db <- initialize_stringdb(
      version = "11", species = 9606, score_threshold = 400,
      use_cache = TRUE, cache_dir = "data/stringdb_cache"
    )
    cat("STRING database initialized\n\n")
    return(string_db)
  }

  cat("STRING DB filtering: DISABLED\n")
  cat("All genes will be included in output\n\n")
  NULL
}

#' Calculate comparison statistics between baseline and current analysis
calculate_comparison_stats <- function(baseline_difexp_all, difexp_all_final,
                                        baseline_genes_merged, normalize_fdr) {
  if (is.null(baseline_difexp_all)) return(NULL)

  gene_col <- if ("SYMBOL" %in% colnames(baseline_difexp_all)) "SYMBOL" else "ENTREZID"
  baseline_data <- baseline_difexp_all[, c(gene_col, "logFC", "adj.P.Val", "P.Value")]
  colnames(baseline_data) <- c("gene", "logFC_baseline", "FDR_baseline", "P_baseline")
  current_data <- difexp_all_final[, c(gene_col, "logFC", "adj.P.Val", "P.Value")]
  colnames(current_data) <- c("gene", "logFC_current", "FDR_current", "P_current")
  merged_data <- merge(baseline_data, current_data, by = "gene", all = FALSE)

  cat("FDR normalization setting: normalize_fdr =", normalize_fdr, "\n")

  if (normalize_fdr && !is.null(baseline_genes_merged)) {
    baseline_n_genes <- length(baseline_genes_merged)
    cat("Recalculating FDR with baseline gene count for fair comparison...\n")
    cat("  Baseline genes:", baseline_n_genes, ", Current genes:", nrow(difexp_all_final), "\n")
    merged_data$FDR_current_adjusted <- p.adjust(merged_data$P_current, method = "BH", n = baseline_n_genes)
  } else {
    merged_data$FDR_current_adjusted <- merged_data$FDR_current
  }

  list(
    genes_fdr_increased = sum(merged_data$FDR_current_adjusted > merged_data$FDR_baseline),
    genes_fdr_decreased = sum(merged_data$FDR_current_adjusted < merged_data$FDR_baseline),
    genes_lost_significance = sum(merged_data$FDR_baseline < 0.05 & merged_data$FDR_current_adjusted >= 0.05),
    genes_gained_significance = sum(merged_data$FDR_baseline >= 0.05 & merged_data$FDR_current_adjusted < 0.05),
    genes_logfc_change_gt1 = sum(abs(merged_data$logFC_current - merged_data$logFC_baseline) > 1),
    genes_logfc_fell_below1 = sum(abs(merged_data$logFC_baseline) > 1 & abs(merged_data$logFC_current) <= 1),
    genes_logfc_rose_above1 = sum(abs(merged_data$logFC_baseline) <= 1 & abs(merged_data$logFC_current) > 1)
  )
}

#' Save iteration output files
save_iteration_files <- function(iteration_dir, merged_exprs, exprs_corrected,
                                  pdata_aligned, difexp_logfc1_final, difexp_all_final,
                                  difexp_logfc1 = NULL, difexp_all = NULL,
                                  filter_by_stringdb = FALSE) {
  dir.create(iteration_dir, showWarnings = FALSE, recursive = TRUE)

  write.table(merged_exprs, file.path(iteration_dir, "merged_exprs_before_combat.tsv"),
              sep = "\t", quote = FALSE)
  write.table(exprs_corrected, file.path(iteration_dir, "merged_exprs_after_combat.tsv"),
              sep = "\t", quote = FALSE)
  write.csv(pdata_aligned, file.path(iteration_dir, "phenodata.csv"), row.names = FALSE)
  write.csv(difexp_logfc1_final, file.path(iteration_dir, "difexp_genes.csv"), row.names = FALSE)
  write.csv(difexp_all_final, file.path(iteration_dir, "difexp_genes_unfiltered.csv"), row.names = FALSE)

  if (filter_by_stringdb && !is.null(difexp_logfc1) && !is.null(difexp_all)) {
    write.csv(difexp_logfc1, file.path(iteration_dir, "difexp_genes_before_stringdb.csv"), row.names = FALSE)
    write.csv(difexp_all, file.path(iteration_dir, "difexp_genes_unfiltered_before_stringdb.csv"), row.names = FALSE)
  }

  cat("Saved iteration files to:", iteration_dir, "\n")
}

#' Run a single analysis iteration (baseline or with addon datasets)
run_analysis_iteration <- function(config, pdata_full, exclusions, batch_col,
                                    trimester_col_1, trimester_col_2,
                                    string_db, filter_by_stringdb,
                                    baseline_genes_merged = NULL,
                                    baseline_degs = NULL,
                                    baseline_difexp_all = NULL) {

  pdata_filtered <- filter_phenodata(pdata = pdata_full, filters = config$filtering)
  pdata_filtered <- pdata_filtered[!(pdata_filtered[[batch_col]] %in% exclusions), ]

  n_trim_1 <- sum(pdata_filtered$Gestational.Age.Category == trimester_col_1, na.rm = TRUE)
  n_trim_2 <- sum(pdata_filtered$Gestational.Age.Category == trimester_col_2, na.rm = TRUE)
  cat("Samples:", trimester_col_1, "=", n_trim_1, ",", trimester_col_2, "=", n_trim_2, "\n")

  required_datasets <- identify_datasets(pdata_filtered = pdata_filtered, batch_col = batch_col)
  cat("Datasets included:", paste(required_datasets, collapse = ", "), "\n")

  merged_exprs <- merge_selected_datasets(
    mapped_path = config$paths$mapped_data,
    datasets = required_datasets,
    pattern = config$merging$pattern
  )

  pdata_aligned <- pdata_filtered[
    make.names(pdata_filtered[[config$merging$sample_column]]) %in% colnames(merged_exprs),
  ]

  merged_exprs <- subset_expression_by_phenodata(
    mrgd = merged_exprs,
    pdata_filtered = pdata_aligned,
    sample_col = config$merging$sample_column
  )

  cat("Expression matrix:", nrow(merged_exprs), "genes x", ncol(merged_exprs), "samples\n")

  n_trim_1 <- sum(pdata_aligned$Gestational.Age.Category == trimester_col_1, na.rm = TRUE)
  n_trim_2 <- sum(pdata_aligned$Gestational.Age.Category == trimester_col_2, na.rm = TRUE)

  genes_gained_vec <- genes_lost_vec <- NULL
  if (!is.null(baseline_genes_merged)) {
    current_genes <- rownames(merged_exprs)
    genes_gained_vec <- setdiff(current_genes, baseline_genes_merged)
    genes_lost_vec <- setdiff(baseline_genes_merged, current_genes)
    cat("Genes comparison: Common =", length(intersect(current_genes, baseline_genes_merged)),
        ", Gained =", length(genes_gained_vec), ", Lost =", length(genes_lost_vec), "\n")
  }

  pdata_aligned <- add_trimester_variable(pdata_aligned)

  cat("Applying ComBat batch correction...\n")
  exprs_corrected <- apply_combat_correction(
    exprs = merged_exprs, pdata = pdata_aligned, batch_col = batch_col,
    mod_formula = as.formula(config$batch_correction$model_formula),
    par_prior = config$batch_correction$par_prior, prior_plots = FALSE
  )

  cat("Running differential expression...\n")
  de_results <- differential_expression(
    exprs = exprs_corrected, pdata = pdata_aligned,
    design_formula = as.formula(config$differential_expression$design_formula),
    contrasts = config$differential_expression$contrasts,
    method = config$differential_expression$decide_method
  )

  difexp_all <- extract_difexp_results(
    efit = de_results$efit, results = de_results$results,
    exprs = exprs_corrected, pdata = pdata_aligned, species = "human"
  )

  difexp_all_sig <- filter_difexp(difexp = difexp_all, logfc_threshold = 0,
                                   pvalue_threshold = config$differential_expression$pvalue_threshold)
  difexp_logfc1 <- filter_difexp(difexp = difexp_all, logfc_threshold = 1,
                                  pvalue_threshold = config$differential_expression$pvalue_threshold)

  n_degs_all <- nrow(difexp_all_sig)
  n_degs_logfc1 <- nrow(difexp_logfc1)
  cat("DEGs (p < 0.05):", n_degs_all, ", DEGs (|logFC| > 1):", n_degs_logfc1, "\n")

  if (filter_by_stringdb && !is.null(string_db)) {
    cat("Filtering for STRING DB proteins...\n")
    difexp_logfc1_final <- filter_difexp_by_stringdb(difexp_logfc1, string_db, "SYMBOL", TRUE)$filtered
    difexp_all_final <- filter_difexp_by_stringdb(difexp_all, string_db, "SYMBOL")
    cat("DEGs with STRING proteins (|logFC| > 1):", nrow(difexp_logfc1_final), "\n")
  } else {
    difexp_logfc1_final <- difexp_logfc1
    difexp_all_final <- difexp_all
  }

  if (!is.null(baseline_degs)) {
    current_degs <- difexp_logfc1_final$ENTREZID
    cat("DEGs comparison: Common =", length(intersect(current_degs, baseline_degs)),
        ", Gained =", length(setdiff(current_degs, baseline_degs)),
        ", Lost =", length(setdiff(baseline_degs, current_degs)), "\n")
  }

  normalize_fdr <- isTRUE(config$differential_expression$normalize_fdr_to_baseline)
  comparison_stats <- calculate_comparison_stats(baseline_difexp_all, difexp_all_final,
                                                  baseline_genes_merged, normalize_fdr)

  # Calculate combined score: sum of |logFC| * -log10(adj.P.Val) for DEGs
  deg_score <- sum(
    abs(difexp_logfc1_final$logFC) * -log10(difexp_logfc1_final$adj.P.Val + 1e-300)
  )

  list(
    merged_exprs = merged_exprs, exprs_corrected = exprs_corrected,
    pdata_aligned = pdata_aligned, difexp_all = difexp_all, difexp_logfc1 = difexp_logfc1,
    difexp_all_final = difexp_all_final, difexp_logfc1_final = difexp_logfc1_final,
    n_trim_1 = n_trim_1, n_trim_2 = n_trim_2,
    n_degs_all = n_degs_all, n_degs_logfc1 = n_degs_logfc1,
    genes_merged = rownames(merged_exprs), degs = difexp_logfc1_final$ENTREZID,
    comparison_stats = comparison_stats,
    genes_gained = genes_gained_vec, genes_lost = genes_lost_vec,
    deg_score = deg_score
  )
}

#' Run single dataset analysis (no batch correction)
run_single_dataset_analysis <- function(config, pdata_single, dataset_name,
                                         trimester_col_1, trimester_col_2) {
  n_trim_1 <- sum(pdata_single$Gestational.Age.Category == trimester_col_1, na.rm = TRUE)
  n_trim_2 <- sum(pdata_single$Gestational.Age.Category == trimester_col_2, na.rm = TRUE)
  cat("Samples:", trimester_col_1, "=", n_trim_1, ",", trimester_col_2, "=", n_trim_2, "\n")

  merged_exprs <- merge_selected_datasets(
    mapped_path = config$paths$mapped_data, datasets = dataset_name, pattern = config$merging$pattern
  )

  pdata_aligned <- pdata_single[
    make.names(pdata_single[[config$merging$sample_column]]) %in% colnames(merged_exprs),
  ]

  merged_exprs <- subset_expression_by_phenodata(
    mrgd = merged_exprs, pdata_filtered = pdata_aligned, sample_col = config$merging$sample_column
  )

  cat("Expression matrix:", nrow(merged_exprs), "genes x", ncol(merged_exprs), "samples\n")

  n_trim_1 <- sum(pdata_aligned$Gestational.Age.Category == trimester_col_1, na.rm = TRUE)
  n_trim_2 <- sum(pdata_aligned$Gestational.Age.Category == trimester_col_2, na.rm = TRUE)

  pdata_aligned <- add_trimester_variable(pdata_aligned)
  cat("Single dataset - no batch correction needed\n")

  de_results <- differential_expression(
    exprs = merged_exprs, pdata = pdata_aligned,
    design_formula = as.formula(config$differential_expression$design_formula),
    contrasts = config$differential_expression$contrasts,
    method = config$differential_expression$decide_method
  )

  difexp_all <- extract_difexp_results(
    efit = de_results$efit, results = de_results$results,
    exprs = merged_exprs, pdata = pdata_aligned, species = "human"
  )

  difexp_all_sig <- filter_difexp(difexp = difexp_all, logfc_threshold = 0,
                                   pvalue_threshold = config$differential_expression$pvalue_threshold)
  difexp_logfc1 <- filter_difexp(difexp = difexp_all, logfc_threshold = 1,
                                  pvalue_threshold = config$differential_expression$pvalue_threshold)

  cat("DEGs (p < 0.05):", nrow(difexp_all_sig), ", DEGs (|logFC| > 1):", nrow(difexp_logfc1), "\n")

  list(
    difexp_all = difexp_all, difexp_logfc1 = difexp_logfc1,
    n_trim_1 = n_trim_1, n_trim_2 = n_trim_2,
    total_samples = ncol(merged_exprs), n_genes = nrow(merged_exprs),
    n_degs_all = nrow(difexp_all_sig), n_degs_logfc1 = nrow(difexp_logfc1)
  )
}

#' Main function to run dataset testing analysis
#'
#' @param config_file Path to configuration file
#' @param output_dir Output directory path
#' @param trimester_col_1 First trimester category (e.g., "First Trimester", "Second Trimester")
#' @param trimester_col_2 Second trimester category (e.g., "Second Trimester", "Term")
#' @param addon_config_key Config key for addon datasets (e.g., "first_trimester_datasets", "term_datasets")
#' @param file_prefix Prefix for output files (e.g., "first_", "term_")
#' @param type_label Type label for results (e.g., "with_first", "with_term")
#' @param enable_logging Whether to enable logging to file
#' @param use_global_exclusions Whether to use global exclusions from config
run_dataset_comparison <- function(config_file, output_dir,
                                    trimester_col_1, trimester_col_2,
                                    addon_config_key, file_prefix, type_label,
                                    enable_logging = FALSE,
                                    use_global_exclusions = FALSE) {

  if (enable_logging) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    log_file <- file.path(output_dir, paste0("run_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
    sink(log_file, split = TRUE)
    cat("Log file:", log_file, "\nStarted:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  }

  cat("\n=== Testing Datasets:", trimester_col_1, "vs", trimester_col_2, "===\n\n")

  source("scripts/integrative_analysis/integrative_analysis_disser_pipeline/config.R")
  source("scripts/integrative_analysis/integrative_analysis_disser_pipeline/utils.R")
  source("scripts/integrative_analysis/integrative_analysis_disser_pipeline/01_data_merging.R")
  source("scripts/integrative_analysis/integrative_analysis_disser_pipeline/02_batch_correction.R")
  source("scripts/integrative_analysis/integrative_analysis_disser_pipeline/03_pca_analysis.R")
  source("scripts/integrative_analysis/integrative_analysis_disser_pipeline/04_differential_expression.R")
  source("scripts/integrative_analysis/integrative_analysis_disser_pipeline/05_network_clustering.R")

  config <- load_config(config_file)
  archive_previous_results(output_dir)
  cat("Output directory:", output_dir, "\n\n")

  filter_by_stringdb <- isTRUE(config$differential_expression$filter_by_stringdb)
  string_db <- initialize_string_db_if_needed(filter_by_stringdb)

  baseline_datasets <- config$baseline$datasets
  addon_datasets <- config[[addon_config_key]]$datasets
  global_exclusions <- if (use_global_exclusions && !is.null(config$global_exclusions)) {
    config$global_exclusions$datasets
  } else {
    character(0)
  }

  cat("Baseline datasets:", paste(baseline_datasets, collapse = ", "), "\n")
  cat("Additional datasets to test:", paste(addon_datasets, collapse = ", "), "\n")
  if (length(global_exclusions) > 0) cat("Global exclusions:", paste(global_exclusions, collapse = ", "), "\n")
  cat("\n")

  pdata_full <- read.csv(config$paths$phenodata, stringsAsFactors = FALSE)
  batch_col <- config$batch_correction$batch_column

  # Validate dataset files and metadata matching before analysis
  source("tests/check_duplicate_dataset_files.R")
  all_datasets <- unique(c(baseline_datasets, addon_datasets))
  validation <- validate_dataset_files(
    mapped_path = config$paths$mapped_data,
    phenodata_path = config$paths$phenodata,
    datasets = all_datasets,
    pattern = config$merging$pattern,
    sample_col = config$merging$sample_column,
    batch_col = batch_col,
    verbose = TRUE
  )

  if (!validation$valid) {
    stop("Validation failed. Fix data issues before running analysis.")
  }

  col_1_name <- gsub(" ", "_", tolower(trimester_col_1))
  col_2_name <- gsub(" ", "_", tolower(trimester_col_2))

  results_summary <- data.frame(
    dataset = character(), type = character(),
    stringsAsFactors = FALSE
  )
  results_summary[[col_1_name]] <- integer()
  results_summary[[col_2_name]] <- integer()
  results_summary$total_samples <- integer()
  results_summary$n_genes_merged <- integer()
  results_summary$n_degs_all <- integer()
  results_summary$n_degs_logfc1 <- integer()
  results_summary$genes_fdr_increased <- integer()
  results_summary$genes_fdr_decreased <- integer()
  results_summary$genes_lost_significance <- integer()
  results_summary$genes_gained_significance <- integer()
  results_summary$genes_logfc_change_gt1 <- integer()
  results_summary$genes_logfc_fell_below1 <- integer()
  results_summary$genes_logfc_rose_above1 <- integer()
  results_summary$deg_score <- numeric()

  baseline_genes_merged <- baseline_degs <- baseline_difexp_all <- NULL

  # === BASELINE ANALYSIS ===
  cat("\n==================================================================\n")
  cat("           BASELINE ANALYSIS                                     \n")
  cat("==================================================================\n\n")

  baseline_exclusions <- c(addon_datasets, global_exclusions)
  cat("Excluding:", paste(baseline_exclusions, collapse = ", "), "\n\n")

  tryCatch({
    baseline_result <- run_analysis_iteration(
      config, pdata_full, baseline_exclusions, batch_col,
      trimester_col_1, trimester_col_2, string_db, filter_by_stringdb
    )

    baseline_genes_merged <- baseline_result$genes_merged
    baseline_degs <- baseline_result$degs
    baseline_difexp_all <- baseline_result$difexp_all_final

    baseline_name <- paste0("BASELINE (", paste(baseline_datasets, collapse = ", "), ")")
    new_row <- list(
      dataset = baseline_name, type = "baseline",
      total_samples = ncol(baseline_result$merged_exprs),
      n_genes_merged = nrow(baseline_result$merged_exprs),
      n_degs_all = baseline_result$n_degs_all, n_degs_logfc1 = baseline_result$n_degs_logfc1,
      genes_fdr_increased = NA, genes_fdr_decreased = NA,
      genes_lost_significance = NA, genes_gained_significance = NA,
      genes_logfc_change_gt1 = NA, genes_logfc_fell_below1 = NA, genes_logfc_rose_above1 = NA,
      deg_score = baseline_result$deg_score
    )
    new_row[[col_1_name]] <- baseline_result$n_trim_1
    new_row[[col_2_name]] <- baseline_result$n_trim_2
    results_summary <- rbind(results_summary, as.data.frame(new_row, stringsAsFactors = FALSE))

    save_iteration_files(
      file.path(output_dir, "BASELINE"), baseline_result$merged_exprs,
      baseline_result$exprs_corrected, baseline_result$pdata_aligned,
      baseline_result$difexp_logfc1_final, baseline_result$difexp_all_final,
      baseline_result$difexp_logfc1, baseline_result$difexp_all, filter_by_stringdb
    )
    cat("Baseline analysis complete\n")
  }, error = function(e) cat("Error in baseline analysis:", conditionMessage(e), "\n"))

  # === TEST ADDITIONAL DATASETS ===
  cat("\n\n==================================================================\n")
  cat("           TESTING ADDITIONAL DATASETS                           \n")
  cat("==================================================================\n\n")

  if (isTRUE(config$testing$test_all_subsets)) {
    n <- length(addon_datasets)
    datasets_to_test <- lapply(1:(2^n - 1), function(i) addon_datasets[as.logical(intToBits(i)[1:n])])
    cat("Testing mode: ALL SUBSETS (", length(datasets_to_test), " combinations)\n\n")
  } else {
    datasets_to_test <- as.list(addon_datasets)
    cat("Testing mode: INDIVIDUAL DATASETS\n\n")
  }

  for (test_datasets in datasets_to_test) {
    test_name <- paste(test_datasets, collapse = "+")
    cat("\n========================================\n")
    cat("Testing:", test_name, "\n")
    cat("========================================\n\n")

    exclusions <- c(setdiff(addon_datasets, test_datasets), global_exclusions)
    cat("Excluding:", paste(exclusions, collapse = ", "), "\n\n")

    tryCatch({
      iter_result <- run_analysis_iteration(
        config, pdata_full, exclusions, batch_col,
        trimester_col_1, trimester_col_2, string_db, filter_by_stringdb,
        baseline_genes_merged, baseline_degs, baseline_difexp_all
      )

      if (length(iter_result$genes_gained) > 0) {
        writeLines(iter_result$genes_gained,
                   file.path(output_dir, paste0(file_prefix, gsub("\\+", "_", test_name), "_genes_gained.txt")))
      }
      if (length(iter_result$genes_lost) > 0) {
        writeLines(iter_result$genes_lost,
                   file.path(output_dir, paste0(file_prefix, gsub("\\+", "_", test_name), "_genes_lost.txt")))
      }

      cs <- iter_result$comparison_stats
      new_row <- list(
        dataset = test_name, type = type_label,
        total_samples = ncol(iter_result$merged_exprs),
        n_genes_merged = nrow(iter_result$merged_exprs),
        n_degs_all = iter_result$n_degs_all, n_degs_logfc1 = iter_result$n_degs_logfc1,
        genes_fdr_increased = if (!is.null(cs)) cs$genes_fdr_increased else NA,
        genes_fdr_decreased = if (!is.null(cs)) cs$genes_fdr_decreased else NA,
        genes_lost_significance = if (!is.null(cs)) cs$genes_lost_significance else NA,
        genes_gained_significance = if (!is.null(cs)) cs$genes_gained_significance else NA,
        genes_logfc_change_gt1 = if (!is.null(cs)) cs$genes_logfc_change_gt1 else NA,
        genes_logfc_fell_below1 = if (!is.null(cs)) cs$genes_logfc_fell_below1 else NA,
        genes_logfc_rose_above1 = if (!is.null(cs)) cs$genes_logfc_rose_above1 else NA,
        deg_score = iter_result$deg_score
      )
      new_row[[col_1_name]] <- iter_result$n_trim_1
      new_row[[col_2_name]] <- iter_result$n_trim_2
      results_summary <- rbind(results_summary, as.data.frame(new_row, stringsAsFactors = FALSE))

      iteration_dir <- file.path(output_dir, gsub("\\+", "_", test_name))
      save_iteration_files(
        iteration_dir, iter_result$merged_exprs, iter_result$exprs_corrected,
        iter_result$pdata_aligned, iter_result$difexp_logfc1_final,
        iter_result$difexp_all_final, iter_result$difexp_logfc1,
        iter_result$difexp_all, filter_by_stringdb
      )

      if (!is.null(baseline_difexp_all)) {
        tryCatch({
          plot_stats <- plot_logfc_fdr_comparison(
            baseline_difexp_all, iter_result$difexp_all_final,
            file.path(iteration_dir, "logFC_FDR_comparison.png"),
            paste0("Baseline vs ", test_name)
          )
          cat("Plot stats - Common:", plot_stats$n_genes, ", Mean |logFC| change:",
              round(plot_stats$mean_logFC_change, 3), "\n")
        }, error = function(e) cat("Could not generate plot:", conditionMessage(e), "\n"))
      }
      cat("Analysis complete\n")
    }, error = function(e) cat("Error with", test_name, ":", conditionMessage(e), "\n"))
  }

  # === SINGLE DATASET ANALYSES ===
  cat("\n\n==================================================================\n")
  cat("           SINGLE DATASET ANALYSES                               \n")
  cat("==================================================================\n\n")

  pdata_for_singles <- filter_phenodata(pdata = pdata_full, filters = config$filtering)
  dataset_counts <- table(pdata_for_singles[[batch_col]], pdata_for_singles$Gestational.Age.Category)

  datasets_with_both <- character(0)
  if (trimester_col_1 %in% colnames(dataset_counts) && trimester_col_2 %in% colnames(dataset_counts)) {
    datasets_with_both <- rownames(dataset_counts)[
      dataset_counts[, trimester_col_1] >= 3 & dataset_counts[, trimester_col_2] >= 3
    ]
  }
  cat("Datasets with both categories (>=3 each):", paste(datasets_with_both, collapse = ", "), "\n\n")

  single_results <- data.frame(dataset = character(), stringsAsFactors = FALSE)
  single_results[[col_1_name]] <- integer()
  single_results[[col_2_name]] <- integer()
  single_results$total_samples <- integer()
  single_results$n_genes <- integer()
  single_results$n_degs_all <- integer()
  single_results$n_degs_logfc1 <- integer()

  for (single_ds in datasets_with_both) {
    cat("\n========================================\n")
    cat("Analyzing dataset:", single_ds, "\n")
    cat("========================================\n\n")

    pdata_single <- pdata_for_singles[pdata_for_singles[[batch_col]] == single_ds, ]

    tryCatch({
      single_result <- run_single_dataset_analysis(
        config, pdata_single, single_ds, trimester_col_1, trimester_col_2
      )

      write.csv(single_result$difexp_logfc1,
                file.path(output_dir, paste0(file_prefix, single_ds, "_single_difexp.csv")), row.names = FALSE)
      write.csv(single_result$difexp_all,
                file.path(output_dir, paste0(file_prefix, single_ds, "_single_difexp_unfiltered.csv")), row.names = FALSE)

      new_row <- list(dataset = single_ds, total_samples = single_result$total_samples,
                      n_genes = single_result$n_genes, n_degs_all = single_result$n_degs_all,
                      n_degs_logfc1 = single_result$n_degs_logfc1)
      new_row[[col_1_name]] <- single_result$n_trim_1
      new_row[[col_2_name]] <- single_result$n_trim_2
      single_results <- rbind(single_results, as.data.frame(new_row, stringsAsFactors = FALSE))
      cat("Single dataset analysis complete\n")
    }, error = function(e) cat("Error with", single_ds, ":", conditionMessage(e), "\n"))
  }

  if (nrow(single_results) > 0) {
    cat("\n\n==================================================================\n")
    cat("           SINGLE DATASET SUMMARY                                \n")
    cat("==================================================================\n\n")
    print(single_results, row.names = FALSE)
    write.csv(single_results, file.path(output_dir, paste0(file_prefix, "single_dataset_summary.csv")), row.names = FALSE)
  }

  # === SUMMARY ===
  cat("\n\n==================================================================\n")
  cat("                        SUMMARY RESULTS                           \n")
  cat("==================================================================\n\n")

  if (nrow(results_summary) > 0) {
    # Sort by deg_score (best to worst)
    results_summary <- results_summary[order(-results_summary$deg_score), ]
    print(results_summary, row.names = FALSE)
    summary_csv <- file.path(output_dir,
                             paste0(file_prefix, "dataset_comparison_summary.csv"))
    summary_xlsx <- file.path(output_dir,
                              paste0(file_prefix, "dataset_comparison_summary.xlsx"))
    write.csv(results_summary, summary_csv, row.names = FALSE)
    openxlsx::write.xlsx(results_summary, summary_xlsx)

    with_addons <- results_summary[results_summary$type == type_label, ]
    if (nrow(with_addons) > 0) {
      with_addons <- with_addons[order(-with_addons$deg_score), ]
      cat("\n==================================================================\n")
      cat("BEST ADDITIONAL DATASET(S):", with_addons$dataset[1], "\n")
      cat("DEGs (|logFC| > 1):", with_addons$n_degs_logfc1[1], "\n")
      cat("DEG score:", round(with_addons$deg_score[1], 2), "\n")
      cat("==================================================================\n\n")
    }
  }

  cat("\n==================================================================\n")
  cat("                    ANALYSIS COMPLETE                             \n")
  cat("==================================================================\n")

  if (enable_logging) {
    cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
    sink()
  }

  list(results_summary = results_summary, single_results = single_results)
}
