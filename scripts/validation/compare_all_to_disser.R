#!/usr/bin/env Rscript

#' Compare Both 1_2 and 2_3 Results to Dissertation Data
#'
#' Runs dissertation comparison for both trimester comparisons
#' using the shared comparison function.
#'
#' @author Expression Integration Pipeline
#' @date 2025-01-07

# Source the common comparison function
source("R/compare_to_disser_common.R")

# Dissertation file path
disser_file <- "data/disser/lykhenko_Supplement1_sorted_include_present.xlsx"

# ===================================================================
# COMPARISON 1: 1_2 (First vs Second Trimester) - Sheet 1
# ===================================================================

cat("\n")
cat("######################################################################\n")
cat("#                    1_2 COMPARISON (Sheet 1)                        #\n")
cat("######################################################################\n")

summary_1_2 <- run_dissertation_comparison(
  results_dir = "output/test_first_datasets/BASELINE",
  disser_file = disser_file,
  disser_sheet = 1,
  output_dir = "output/comparison_1_2_vs_disser",
  comparison_name = "Trim 1_2"
)

# ===================================================================
# COMPARISON 2: 2_3 (Second Trimester vs Term) - Sheet 2
# ===================================================================

cat("\n")
cat("######################################################################\n")
cat("#                    2_3 COMPARISON (Sheet 2)                        #\n")
cat("######################################################################\n")

summary_2_3 <- run_dissertation_comparison(
  results_dir = "output/test_term_datasets/archive/2026-01-10_214735/BASELINE",
  disser_file = disser_file,
  disser_sheet = 2,
  output_dir = "output/comparison_2_3_vs_disser",
  comparison_name = "Trim 2_3"
)

# ===================================================================
# COMBINED SUMMARY
# ===================================================================

cat("\n")
cat("######################################################################\n")
cat("#                      COMBINED SUMMARY                              #\n")
cat("######################################################################\n\n")

cat("1_2 Comparison:\n")
cat("  Dissertation overlap:", summary_1_2$Value[summary_1_2$Metric == "Overlap (% of dissertation)"], "%\n")
cat("  Pearson correlation:", summary_1_2$Value[summary_1_2$Metric == "Pearson correlation"], "\n")
cat("  Direction concordance:", summary_1_2$Value[summary_1_2$Metric == "Direction concordance (%)"], "%\n\n")

cat("2_3 Comparison:\n")
cat("  Dissertation overlap:", summary_2_3$Value[summary_2_3$Metric == "Overlap (% of dissertation)"], "%\n")
cat("  Pearson correlation:", summary_2_3$Value[summary_2_3$Metric == "Pearson correlation"], "\n")
cat("  Direction concordance:", summary_2_3$Value[summary_2_3$Metric == "Direction concordance (%)"], "%\n\n")

cat("All comparisons complete!\n")
cat("Results saved to:\n")
cat("  - output/comparison_1_2_vs_disser/\n")
cat("  - output/comparison_2_3_vs_disser/\n\n")
