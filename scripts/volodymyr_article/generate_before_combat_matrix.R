#!/usr/bin/env Rscript

#' Generate Before-ComBat Expression Matrix
#'
#' Generates the before-ComBat expression matrix for a specific dataset combination
#' This is needed when the original test run didn't save pre-ComBat matrices
#'
#' @author Expression Integration Pipeline
#' @date 2025-11-13

cat("\n=== Generating Before-ComBat Expression Matrix ===\n\n")

# Load required libraries
library(yaml)

# Source all modules
source("scripts/integrative_analysis/integrative_analysis_disser_pipeline/config.R")
source("scripts/integrative_analysis/integrative_analysis_disser_pipeline/utils.R")
source("scripts/integrative_analysis/integrative_analysis_disser_pipeline/01_data_merging.R")
source("scripts/integrative_analysis/integrative_analysis_disser_pipeline/02_batch_correction.R")

# Load configuration
config <- load_config("config/config_test_first_datasets.yaml")

# Define the specific dataset combination we want
target_dir <- "output/test_first_datasets/GSE55439_GSE93520_GSE28551_GSE100051"

cat("Target directory:", target_dir, "\n\n")

# Check if after-ComBat file exists
if (!file.exists(file.path(target_dir, "merged_exprs.tsv"))) {
  stop("Target directory does not contain merged_exprs.tsv")
}

# Check if before-ComBat file already exists
before_combat_file <- file.path(target_dir, "merged_exprs_before_combat.tsv")
if (file.exists(before_combat_file)) {
  cat("✓ Before-ComBat matrix already exists:", before_combat_file, "\n")
  cat("Delete it first if you want to regenerate.\n")
  quit(status = 0)
}

# Load the phenodata from the target directory
cat("Loading phenodata...\n")
pdata <- read.csv(file.path(target_dir, "phenodata.csv"), stringsAsFactors = FALSE)
cat("Samples:", nrow(pdata), "\n")

# Get unique datasets from phenodata
datasets <- unique(pdata$accession)
cat("Datasets:", paste(datasets, collapse = ", "), "\n\n")

# Define batch column
batch_col <- config$batch_correction$batch_column

# Identify required datasets
cat("Merging expression data...\n")
merged_exprs <- merge_selected_datasets(
  mapped_path = config$paths$mapped_data,
  datasets = datasets,
  pattern = config$merging$pattern
)

cat("Expression matrix:", nrow(merged_exprs), "genes x", ncol(merged_exprs), "samples\n")

# Subset expression to match phenodata
# Find the sample column that matches
sample_col <- if ("arraydatafile_exprscolumnnames" %in% colnames(pdata)) {
  "arraydatafile_exprscolumnnames"
} else if ("exprs_column_names" %in% colnames(pdata)) {
  "exprs_column_names"
} else {
  "Array.Data.File"
}

cat("Using sample column:", sample_col, "\n")

pdata_aligned <- pdata[
  make.names(pdata[[sample_col]]) %in% colnames(merged_exprs),
]

merged_exprs <- subset_expression_by_phenodata(
  mrgd = merged_exprs,
  pdata_filtered = pdata_aligned,
  sample_col = sample_col
)

cat("After subsetting:", nrow(merged_exprs), "genes x", ncol(merged_exprs), "samples\n")

# Verify sample order matches
if (!all(colnames(merged_exprs) == make.names(pdata_aligned[[sample_col]]))) {
  cat("⚠ Warning: Sample order mismatch, reordering...\n")
  merged_exprs <- merged_exprs[, make.names(pdata_aligned[[sample_col]])]
}

# Save before-ComBat matrix
cat("\nSaving before-ComBat expression matrix...\n")
write.table(merged_exprs,
            before_combat_file,
            sep = "\t", quote = FALSE)

cat("✓ Saved:", before_combat_file, "\n")

# Also rename the old merged_exprs.tsv to make it clear it's after ComBat
after_combat_file <- file.path(target_dir, "merged_exprs_after_combat.tsv")
if (!file.exists(after_combat_file)) {
  cat("\nRenaming merged_exprs.tsv to merged_exprs_after_combat.tsv...\n")
  file.rename(
    file.path(target_dir, "merged_exprs.tsv"),
    after_combat_file
  )
  cat("✓ Renamed to:", after_combat_file, "\n")
}

cat("\n")
cat("==================================================================\n")
cat("                    GENERATION COMPLETE                           \n")
cat("==================================================================\n\n")

cat("You can now run pca_analysis_combat_comparison.R to generate both before and after ComBat PCA plots.\n\n")
