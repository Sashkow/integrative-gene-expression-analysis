#' Phase 1: Tiered Gene Coverage Analysis
#'
#' Functions for analyzing and merging expression data with flexible
#' gene coverage thresholds instead of requiring genes in ALL datasets.
#'
#' This module extends the disser pipeline with tiered merging capabilities.
#'
#' @author Expression Integration Pipeline
#' @date 2025

# ============================================================================
# Dependencies - reuse disser pipeline base functions
# ============================================================================

# Get the directory of this script for relative sourcing
.get_script_dir <- function() {
  # Try multiple methods to find script location
  if (exists("ofile", envir = parent.frame())) {
    return(dirname(parent.frame()$ofile))
  }

  # For interactive use, try to find relative to working directory
  candidates <- c(
    "scripts/integrative_analysis/phase1_tiered_coverage",
    "phase1_tiered_coverage",
    "."
  )

  for (dir in candidates) {
    if (dir.exists(dir)) return(dir)
  }

  return(".")
}

# Source disser pipeline utilities if not already loaded
.source_disser_if_needed <- function() {
  # Check if key function exists
  if (!exists("merge_expression_data", mode = "function")) {
    disser_path <- "scripts/integrative_analysis/integrative_analysis_disser_pipeline"

    if (!dir.exists(disser_path)) {
      # Try parent directory
      disser_path <- file.path("..", "integrative_analysis_disser_pipeline")
    }

    if (dir.exists(disser_path)) {
      source(file.path(disser_path, "01_data_merging.R"))
      message("Loaded disser pipeline: 01_data_merging.R")
    }
  }
}

# Auto-load dependencies when this file is sourced
.source_disser_if_needed()


# ============================================================================
# Gene Filtering Functions
# ============================================================================

#' Filter genes by type (protein-coding, ncRNA, etc.)
#'
#' @param datasets Named list of expression data frames
#' @param gene_types Character vector of gene types to keep
#' @return Filtered datasets list
#' @export
filter_genes_by_type <- function(datasets, gene_types = "protein-coding") {

  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Package org.Hs.eg.db required for gene type filtering")
  }

  # Get all unique genes
 all_genes <- unique(unlist(lapply(datasets, rownames)))
  cat("Total genes before filtering:", length(all_genes), "\n")

  # Map to gene types
  mapped_types <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = all_genes,
    column = "GENETYPE",
    keytype = "ENTREZID",
    multiVals = "first"
  )

  # Filter to requested types
 genes_to_keep <- names(mapped_types[mapped_types %in% gene_types &
                                       !is.na(mapped_types)])
  cat("Genes matching types [", paste(gene_types, collapse = ", "), "]:",
      length(genes_to_keep), "\n")

  # Filter each dataset
  filtered_datasets <- lapply(datasets, function(ds) {
    keep_rows <- rownames(ds) %in% genes_to_keep
    ds[keep_rows, , drop = FALSE]
  })

  names(filtered_datasets) <- names(datasets)
  return(filtered_datasets)
}


#' Get protein-coding genes only
#'
#' @param datasets Named list of expression data frames
#' @return Filtered datasets with only protein-coding genes
#' @export
filter_protein_coding <- function(datasets) {
  filter_genes_by_type(datasets, gene_types = "protein-coding")
}


# ============================================================================
# Phase 1: Tiered Gene Coverage Functions
# ============================================================================

#' Get gene presence matrix across datasets
#'
#' Creates a binary matrix showing which genes are present in which datasets
#'
#' @param datasets Named list of expression data frames (genes as rownames)
#' @return Binary matrix with genes as rows, datasets as columns
#' @export
get_gene_presence_matrix <- function(datasets) {

  if (length(datasets) == 0) {
    stop("No datasets provided")
  }

  # Get all unique genes across all datasets
  all_genes <- unique(unlist(lapply(datasets, rownames)))

  # Create presence matrix
  presence_matrix <- matrix(
    0L,
    nrow = length(all_genes),
    ncol = length(datasets),
    dimnames = list(all_genes, names(datasets))
  )

  # Fill in presence (1) for each gene in each dataset
  for (dataset_name in names(datasets)) {
    genes_in_dataset <- rownames(datasets[[dataset_name]])
    presence_matrix[genes_in_dataset, dataset_name] <- 1L
  }

  return(presence_matrix)
}


#' Calculate gene coverage across datasets
#'
#' Computes what percentage of datasets contain each gene
#'
#' @param datasets Named list of expression data frames (genes as rownames)
#' @return Named numeric vector with coverage proportion for each gene
#' @export
calculate_gene_coverage <- function(datasets) {

  if (length(datasets) == 0) {
    stop("No datasets provided")
  }

  presence_matrix <- get_gene_presence_matrix(datasets)

  # Calculate proportion of datasets containing each gene
  coverage <- rowSums(presence_matrix) / ncol(presence_matrix)

  return(coverage)
}


#' Filter genes by coverage threshold
#'
#' Returns genes that meet the minimum coverage threshold
#'
#' @param coverage Named numeric vector of gene coverage proportions
#' @param threshold Minimum coverage threshold (0-1), default 1.0 (100%)
#' @return Character vector of gene names meeting threshold
#' @export
filter_genes_by_coverage <- function(coverage, threshold = 1.0) {

  passing_genes <- names(coverage[coverage >= threshold])

  return(passing_genes)
}


#' Generate comprehensive gene coverage report
#'
#' Creates a detailed report of gene coverage across datasets
#'
#' @param datasets Named list of expression data frames
#' @return List containing presence_matrix, gene_coverage, summary, and tier_counts
#' @export
generate_coverage_report <- function(datasets) {

  presence_matrix <- get_gene_presence_matrix(datasets)
  gene_coverage <- rowSums(presence_matrix) / ncol(presence_matrix)

  # Calculate tier counts
  n_datasets <- ncol(presence_matrix)
  tier_counts <- list(
    tier1_100pct = sum(gene_coverage >= 1.0),
    tier2_75pct = sum(gene_coverage >= 0.75),
    tier3_50pct = sum(gene_coverage >= 0.50),
    total_genes = length(gene_coverage)
  )

  # Summary statistics
  summary_stats <- list(
    n_datasets = n_datasets,
    n_total_genes = length(gene_coverage),
    n_universal_genes = tier_counts$tier1_100pct,
    pct_universal = round(tier_counts$tier1_100pct / length(gene_coverage) * 100, 1),
    mean_coverage = round(mean(gene_coverage), 3),
    median_coverage = round(median(gene_coverage), 3)
  )

  report <- list(
    presence_matrix = presence_matrix,
    gene_coverage = gene_coverage,
    summary = summary_stats,
    tier_counts = tier_counts
  )

  return(report)
}


#' Print gene coverage report
#'
#' @param report Coverage report from generate_coverage_report()
#' @export
print_coverage_report <- function(report) {

  cat("\n")
  cat("=======================================================\n")
  cat("           GENE COVERAGE REPORT\n")
  cat("=======================================================\n\n")

  cat("Datasets analyzed:", report$summary$n_datasets, "\n")
  cat("Total unique genes:", report$summary$n_total_genes, "\n\n")

  cat("Coverage Statistics:\n")
  cat("  Mean coverage:   ", sprintf("%.1f%%", report$summary$mean_coverage * 100), "\n")
  cat("  Median coverage: ", sprintf("%.1f%%", report$summary$median_coverage * 100), "\n\n")

  cat("Tiered Gene Counts:\n")
  cat("  Tier 1 (100% - all datasets):  ", report$tier_counts$tier1_100pct,
      sprintf(" (%.1f%%)", report$tier_counts$tier1_100pct / report$summary$n_total_genes * 100), "\n")
  cat("  Tier 2 (>=75% of datasets):    ", report$tier_counts$tier2_75pct,
      sprintf(" (%.1f%%)", report$tier_counts$tier2_75pct / report$summary$n_total_genes * 100), "\n")
  cat("  Tier 3 (>=50% of datasets):    ", report$tier_counts$tier3_50pct,
      sprintf(" (%.1f%%)", report$tier_counts$tier3_50pct / report$summary$n_total_genes * 100), "\n\n")

  cat("Gene Recovery Improvement:\n")
  cat("  Current method (Tier 1): ", report$tier_counts$tier1_100pct, " genes\n")
  cat("  Tier 2 recovers:         +", report$tier_counts$tier2_75pct - report$tier_counts$tier1_100pct,
      " genes (", sprintf("%.0f%%", (report$tier_counts$tier2_75pct / report$tier_counts$tier1_100pct - 1) * 100),
      " improvement)\n")
  cat("  Tier 3 recovers:         +", report$tier_counts$tier3_50pct - report$tier_counts$tier1_100pct,
      " genes (", sprintf("%.0f%%", (report$tier_counts$tier3_50pct / report$tier_counts$tier1_100pct - 1) * 100),
      " improvement)\n")

  cat("\n=======================================================\n\n")
}


#' Merge datasets with tiered gene coverage threshold
#'
#' Merges multiple expression datasets keeping genes that meet
#' the specified coverage threshold. Missing values are set to NA.
#'
#' @param datasets Named list of expression data frames (genes as rownames)
#' @param coverage_threshold Minimum proportion of datasets gene must be in (0-1)
#' @return Merged expression matrix with NA for missing values
#' @export
merge_datasets_tiered <- function(datasets, coverage_threshold = 1.0) {

  if (length(datasets) == 0) {
    stop("No datasets provided")
  }

  if (coverage_threshold < 0 || coverage_threshold > 1) {
    stop("coverage_threshold must be between 0 and 1")
  }

  # Calculate coverage
  coverage <- calculate_gene_coverage(datasets)

  # Filter genes by threshold
  genes_to_keep <- filter_genes_by_coverage(coverage, coverage_threshold)

  if (length(genes_to_keep) == 0) {
    stop("No genes meet the coverage threshold of ", coverage_threshold * 100, "%")
  }

  message("Tiered merge at ", coverage_threshold * 100, "% threshold:")
  message("  ", length(genes_to_keep), " genes meet threshold (of ", length(coverage), " total)")

  # Get all sample names
  all_samples <- unlist(lapply(datasets, colnames))

  # Create empty merged matrix
  merged <- matrix(
    NA_real_,
    nrow = length(genes_to_keep),
    ncol = length(all_samples),
    dimnames = list(genes_to_keep, all_samples)
  )

  # Fill in values from each dataset
  for (dataset_name in names(datasets)) {
    dataset <- datasets[[dataset_name]]
    samples_in_dataset <- colnames(dataset)
    genes_in_dataset <- intersect(genes_to_keep, rownames(dataset))

    if (length(genes_in_dataset) > 0) {
      merged[genes_in_dataset, samples_in_dataset] <- as.matrix(dataset[genes_in_dataset, , drop = FALSE])
    }
  }

  # Report NA statistics
  na_count <- sum(is.na(merged))
  na_pct <- round(na_count / length(merged) * 100, 1)
  message("  Missing values: ", na_count, " (", na_pct, "% of matrix)")

  return(as.data.frame(merged))
}


#' Merge expression files with tiered gene coverage
#'
#' File-based version of merge_datasets_tiered for use with mapped expression files
#'
#' @param mapped_path Path to directory containing mapped expression files
#' @param pattern File pattern to match (default: "*.tsv")
#' @param coverage_threshold Minimum proportion of datasets gene must be in (0-1)
#' @param datasets Optional: specific dataset names to include (NULL = all)
#' @return List with merged expression matrix and coverage report
#' @export
merge_expression_data_tiered <- function(mapped_path, pattern = "\\.tsv$",
                                         coverage_threshold = 1.0, datasets = NULL) {

  # Get list of expression files
  exprs_files <- list.files(mapped_path, pattern = pattern, full.names = FALSE)

  if (length(exprs_files) == 0) {
    stop("No expression files found in ", mapped_path)
  }

  # Filter to specific datasets if provided
  if (!is.null(datasets)) {
    selected_files <- character(0)
    for (dataset in datasets) {
      matching <- grep(dataset, exprs_files, value = TRUE, fixed = TRUE)
      selected_files <- c(selected_files, matching)
    }
    exprs_files <- unique(selected_files)
  }

  message("Loading ", length(exprs_files), " expression files...")

  # Load all datasets into a list
  datasets_list <- list()
  for (file_name in exprs_files) {
    dataset_name <- gsub("\\.tsv$", "", file_name)
    datasets_list[[dataset_name]] <- read.table(
      file.path(mapped_path, file_name),
      header = TRUE, sep = '\t', row.names = 1
    )
    message("  Loaded ", dataset_name, ": ", nrow(datasets_list[[dataset_name]]), " genes, ",
            ncol(datasets_list[[dataset_name]]), " samples")
  }

  # Generate coverage report
  report <- generate_coverage_report(datasets_list)
  print_coverage_report(report)

  # Perform tiered merge
  merged <- merge_datasets_tiered(datasets_list, coverage_threshold)

  return(list(
    exprs = merged,
    coverage_report = report
  ))
}


#' Load datasets from files into a list
#'
#' Utility function to load expression files into a named list for use
#' with tiered coverage functions.
#'
#' @param mapped_path Path to directory containing mapped expression files
#' @param pattern File pattern to match (default: "*.tsv")
#' @param datasets Optional: specific dataset names to include (NULL = all)
#' @return Named list of expression data frames
#' @export
load_datasets_as_list <- function(mapped_path, pattern = "\\.tsv$", datasets = NULL) {

  # Get list of expression files
  exprs_files <- list.files(mapped_path, pattern = pattern, full.names = FALSE)

  if (length(exprs_files) == 0) {
    stop("No expression files found in ", mapped_path)
  }

  # Filter to specific datasets if provided
  if (!is.null(datasets)) {
    selected_files <- character(0)
    for (dataset in datasets) {
      matching <- grep(dataset, exprs_files, value = TRUE, fixed = TRUE)
      selected_files <- c(selected_files, matching)
    }
    exprs_files <- unique(selected_files)
  }

  message("Loading ", length(exprs_files), " expression files...")

  # Load all datasets into a list
  datasets_list <- list()
  for (file_name in exprs_files) {
    dataset_name <- gsub("\\.tsv$", "", file_name)
    datasets_list[[dataset_name]] <- read.table(
      file.path(mapped_path, file_name),
      header = TRUE, sep = '\t', row.names = 1
    )
    message("  Loaded ", dataset_name, ": ", nrow(datasets_list[[dataset_name]]), " genes, ",
            ncol(datasets_list[[dataset_name]]), " samples")
  }

  return(datasets_list)
}
