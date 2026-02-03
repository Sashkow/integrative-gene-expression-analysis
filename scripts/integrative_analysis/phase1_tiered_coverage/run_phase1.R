#!/usr/bin/env Rscript

#' Run Phase 1: Tiered Gene Coverage Analysis
#'
#' Usage:
#'   Rscript run_phase1.R [options]
#'
#' Options:
#'   --config=PATH         Path to config YAML (default: config_phase1.yaml)
#'   --threshold=NUM       Override coverage threshold 0-1
#'   --output_dir=PATH     Override output directory
#'   --report_only         Only generate report, don't merge
#'   --no_archive          Don't archive previous results
#'   --help                Show this help message
#'
#' Examples:
#'   Rscript run_phase1.R
#'   Rscript run_phase1.R --threshold=0.5
#'   Rscript run_phase1.R --config=my_config.yaml
#'   Rscript run_phase1.R --report_only

library(yaml)

# Null coalescing operator
`%||%` <- function(a, b) if (is.null(a)) b else a

# Source shared utilities
utils_path <- "scripts/utils/logging_utils.R"
if (!file.exists(utils_path)) {
  utils_path <- file.path("..", "..", "utils", "logging_utils.R")
}
source(utils_path)

# Get script directory
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  return("scripts/integrative_analysis/phase1_tiered_coverage")
}

script_dir <- get_script_dir()

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Defaults
config_path <- file.path(script_dir, "config_phase1.yaml")
override_threshold <- NULL
override_output <- NULL
override_report_only <- NULL
no_archive <- FALSE

# Parse arguments
for (arg in args) {
  if (arg == "--help" || arg == "-h") {
    cat("
Phase 1: Tiered Gene Coverage Analysis

Usage:
  Rscript run_phase1.R [options]

Options:
  --config=PATH         Path to config YAML (default: config_phase1.yaml)
  --threshold=NUM       Override coverage threshold 0-1
  --output_dir=PATH     Override output directory
  --report_only         Only generate report, don't merge
  --help                Show this help message

Examples:
  Rscript run_phase1.R
  Rscript run_phase1.R --threshold=0.5
  Rscript run_phase1.R --config=my_config.yaml
  Rscript run_phase1.R --report_only
")
    quit(status = 0)
  } else if (grepl("^--config=", arg)) {
    config_path <- sub("^--config=", "", arg)
  } else if (grepl("^--threshold=", arg)) {
    override_threshold <- as.numeric(sub("^--threshold=", "", arg))
  } else if (grepl("^--output_dir=", arg)) {
    override_output <- sub("^--output_dir=", "", arg)
  } else if (arg == "--report_only") {
    override_report_only <- TRUE
  } else if (arg == "--no_archive") {
    no_archive <- TRUE
  }
}

# Load config
if (!file.exists(config_path)) {
  stop("Config file not found: ", config_path)
}
config <- yaml::read_yaml(config_path)

# Apply overrides
if (!is.null(override_threshold)) config$coverage$threshold <- override_threshold
if (!is.null(override_output)) config$paths$output <- override_output
if (!is.null(override_report_only)) config$analysis$report_only <- override_report_only

# Extract config values
mapped_path <- config$paths$mapped_data
output_dir <- config$paths$output
threshold <- config$coverage$threshold
report_only <- config$analysis$report_only %||% FALSE
verbose <- config$logging$verbose %||% TRUE
pattern <- config$files$pattern %||% "\\.tsv$"
datasets_filter <- config$files$datasets

# Validate threshold
if (is.na(threshold) || threshold < 0 || threshold > 1) {
  stop("Invalid threshold. Must be a number between 0 and 1.")
}

cat("\n")
cat("============================================================\n")
cat("  Phase 1: Tiered Gene Coverage Analysis\n")
cat("============================================================\n\n")

cat("Configuration:\n")
cat("  Config file:  ", config_path, "\n")
cat("  Mapped path:  ", mapped_path, "\n")
cat("  Threshold:    ", threshold * 100, "%\n")
cat("  Output dir:   ", output_dir, "\n")
cat("  Report only:  ", report_only, "\n")
if (!is.null(datasets_filter)) {
  cat("  Datasets:     ", paste(datasets_filter, collapse = ", "), "\n")
}
cat("\n")

# Source the module
source(file.path(script_dir, "tiered_coverage.R"))

# Archive previous results (unless --no_archive)
if (!no_archive) {
  archive_previous_results(output_dir)
} else {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}

# Setup logging
log_file <- NULL
if (config$logging$save_log %||% TRUE) {
  log_file <- setup_logging(output_dir, prefix = "phase1_log")
}

# Load datasets
cat("Loading datasets...\n")
datasets <- load_datasets_as_list(mapped_path, pattern = pattern, datasets = datasets_filter)

# Apply gene type filter if configured
protein_coding_only <- config$gene_filter$protein_coding_only %||% FALSE
if (protein_coding_only) {
  cat("\nFiltering to protein-coding genes only...\n")
  datasets <- filter_protein_coding(datasets)
} else if (!is.null(config$gene_filter$include_types)) {
  gene_types <- config$gene_filter$include_types
  cat("\nFiltering to gene types:", paste(gene_types, collapse = ", "), "\n")
  datasets <- filter_genes_by_type(datasets, gene_types)
}

# Generate coverage report
cat("\nGenerating coverage report...\n")
report <- generate_coverage_report(datasets)
print_coverage_report(report)

# Compare thresholds if enabled
if (config$analysis$compare_thresholds %||% FALSE) {
  cat("\n=== Threshold Comparison ===\n\n")
  coverage <- report$gene_coverage
  comparison_thresholds <- config$analysis$comparison_thresholds %||% c(1.0, 0.75, 0.5)

  cat(sprintf("%-12s %8s %8s\n", "Threshold", "Genes", "Pct"))
  cat(paste(rep("-", 30), collapse = ""), "\n")
  for (thresh in comparison_thresholds) {
    n_genes <- sum(coverage >= thresh)
    pct <- round(n_genes / length(coverage) * 100, 1)
    n_datasets_req <- ceiling(thresh * report$summary$n_datasets)
    cat(sprintf("%3.0f%% (>=%2d ds) %7d %7.1f%%\n",
                thresh * 100, n_datasets_req, n_genes, pct))
  }
  cat("\n")
}

# Save outputs
if (config$output$save_coverage_report %||% TRUE) {
  coverage_df <- data.frame(
    gene_id = names(report$gene_coverage),
    coverage = report$gene_coverage,
    n_datasets = rowSums(report$presence_matrix),
    stringsAsFactors = FALSE
  )
  coverage_df <- coverage_df[order(-coverage_df$coverage, coverage_df$gene_id), ]

  coverage_file <- file.path(output_dir, "gene_coverage.tsv")
  write.table(coverage_df, coverage_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Saved gene coverage to:", coverage_file, "\n")
}

if (config$output$save_presence_matrix %||% TRUE) {
  presence_file <- file.path(output_dir, "presence_matrix.tsv")
  write.table(report$presence_matrix, presence_file, sep = "\t", quote = FALSE)
  cat("Saved presence matrix to:", presence_file, "\n")
}

if (config$output$save_summary %||% TRUE) {
  summary_file <- file.path(output_dir, "summary.txt")

  summary_conn <- file(summary_file, "w")
  writeLines(c(
    "Phase 1: Tiered Gene Coverage Summary",
    "=====================================",
    "",
    paste("Date:", Sys.time()),
    paste("Config:", config_path),
    "",
    paste("Datasets analyzed:", report$summary$n_datasets),
    paste("Total unique genes:", report$summary$n_total_genes),
    "",
    "Coverage Statistics:",
    paste("  Mean coverage:   ", sprintf("%.1f%%", report$summary$mean_coverage * 100)),
    paste("  Median coverage: ", sprintf("%.1f%%", report$summary$median_coverage * 100)),
    "",
    "Tiered Gene Counts:",
    paste("  Tier 1 (100%):", report$tier_counts$tier1_100pct),
    paste("  Tier 2 (75%): ", report$tier_counts$tier2_75pct),
    paste("  Tier 3 (50%): ", report$tier_counts$tier3_50pct)
  ), summary_conn)
  close(summary_conn)

  cat("Saved summary to:", summary_file, "\n")
}

# Perform merge if not report_only
merged_file <- NULL
if (!report_only && (config$output$save_merged %||% TRUE)) {
  cat("\nPerforming tiered merge at", threshold * 100, "% threshold...\n")

  merged <- merge_datasets_tiered(datasets, coverage_threshold = threshold)

  merged_file <- file.path(output_dir, sprintf("merged_%.0fpct.tsv", threshold * 100))
  write.table(merged, merged_file, sep = "\t", quote = FALSE)
  cat("Saved merged expression to:", merged_file, "\n")

  cat("\nMerged matrix dimensions:", nrow(merged), "genes x", ncol(merged), "samples\n")

  na_count <- sum(is.na(merged))
  na_pct <- round(na_count / length(as.matrix(merged)) * 100, 1)
  cat("Missing values:", na_count, "(", na_pct, "%)\n")
}

# Close log
if (!is.null(log_file)) {
  close_logging(log_file)
}

cat("\n============================================================\n")
cat("  Phase 1 Complete!\n")
cat("============================================================\n\n")

cat("Output files in", output_dir, ":\n")
for (f in list.files(output_dir)) {
  cat("  ", f, "\n")
}
cat("\n")
