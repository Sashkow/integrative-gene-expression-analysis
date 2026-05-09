#!/usr/bin/env Rscript

#' Test which dataset gives most DE genes when added to baseline
#'
#' Comparison: Second Trimester vs Term
#'
#' Baseline: Old baseline datasets (GSE122214, GSE22490, GSE37901, GSE9984)
#' Test: Add each additional term dataset and compare genes and DEGs
#'
#' @author Expression Integration Pipeline

cat("\n=== Testing Term Datasets for 2_3 Comparison ===\n\n")

source("scripts/analysis/compare_baseline_addons.R")

run_dataset_comparison(
  config_file = "config/dataset_testing/config_test_term_datasets.yaml",
  output_dir = "output/dataset_testing/test_term_datasets",
  trimester_col_1 = "Second Trimester",
  trimester_col_2 = "Term",
  addon_config_key = "term_datasets",
  file_prefix = "term_",
  type_label = "with_term",
  enable_logging = TRUE,
  use_global_exclusions = TRUE
)
