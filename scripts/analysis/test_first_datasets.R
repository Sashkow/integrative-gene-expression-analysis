#!/usr/bin/env Rscript

#' Test which dataset gives most DE genes when added to baseline
#'
#' Comparison: First Trimester vs Second Trimester
#'
#' Baseline: Old baseline datasets (GSE122214, GSE22490, GSE37901, GSE9984)
#' Test: Add each additional first trimester dataset and compare genes and DEGs
#'
#' @author Expression Integration Pipeline

cat("\n=== Testing Additional Datasets Addition to Old Baseline ===\n\n")

source("scripts/analysis/compare_baseline_addons.R")

run_dataset_comparison(
  config_file = "config/config_test_first_datasets.yaml",
  output_dir = "output/test_first_datasets",
  trimester_col_1 = "First Trimester",
  trimester_col_2 = "Second Trimester",
  addon_config_key = "first_trimester_datasets",
  file_prefix = "first_",
  type_label = "with_first",
  enable_logging = FALSE,
  use_global_exclusions = FALSE
)
