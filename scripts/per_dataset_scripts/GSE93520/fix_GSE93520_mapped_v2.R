#!/usr/bin/env Rscript

#' Fix GSE93520 mapped file TSV format (version 2 - correct approach)
#'
#' The file has 36 sample columns but no empty first header field.
#' We need to add an empty first field to the header line only.
#'
#' @author Expression Integration Pipeline
#' @date 2025-10-26

cat("\n=== Fixing GSE93520 Mapped File (v2) ===\n\n")

# Define paths
original_file <- "/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis/data/mapped/GSE93520_matrix_no_filtering.tsv"
backup_file <- "/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis/data/mapped/GSE93520_matrix_no_filtering.tsv.backup"
temp_file <- "/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis/data/mapped/GSE93520_matrix_no_filtering.tsv.temp"

# ===================================================================
# APPROACH: Add tab to beginning of header line only
# ===================================================================

cat("Reading file and adding tab to header...\n")

# Read all lines
all_lines <- readLines(original_file)

# Modify only the first line (header)
all_lines[1] <- paste0("\t", all_lines[1])

# Write to temp file
writeLines(all_lines, temp_file)

cat("  ✓ Added empty first field to header\n\n")

# ===================================================================
# VERIFY
# ===================================================================

cat("Verifying the fix...\n")

# Load with row.names=1
exprs_verify <- read.table(temp_file, header = TRUE, sep = "\t", row.names = 1,
                           stringsAsFactors = FALSE, check.names = FALSE)

cat("  Dimensions:", nrow(exprs_verify), "genes ×", ncol(exprs_verify), "samples\n")
cat("  First 3 samples:", paste(head(colnames(exprs_verify), 3), collapse = ", "), "\n")
cat("  Last 3 samples:", paste(tail(colnames(exprs_verify), 3), collapse = ", "), "\n")
cat("  Has GSM2453819_CV06.P10:", "GSM2453819_CV06.P10" %in% colnames(exprs_verify), "\n\n")

# Check some expression values
cat("Sample expression values for first gene:\n")
cat("  GSM2453819_CV06.P10:", exprs_verify["1", "GSM2453819_CV06.P10"], "\n")
cat("  GSM2453820_CV06.P18:", exprs_verify["1", "GSM2453820_CV06.P18"], "\n")
cat("  GSM2453821_CV06.P41:", exprs_verify["1", "GSM2453821_CV06.P41"], "\n\n")

# ===================================================================
# FINALIZE
# ===================================================================

if (ncol(exprs_verify) == 36 && "GSM2453819_CV06.P10" %in% colnames(exprs_verify)) {
  cat("✓ Fix successful! Replacing original file...\n")

  # Create backup if it doesn't exist
  if (!file.exists(backup_file)) {
    file.copy(original_file, backup_file)
    cat("  ✓ Backup created:", backup_file, "\n")
  }

  # Replace original with fixed version
  file.copy(temp_file, original_file, overwrite = TRUE)
  file.remove(temp_file)

  cat("  ✓ Original file updated\n\n")

  cat("=================================================================\n")
  cat("                    FIX COMPLETE                                 \n")
  cat("=================================================================\n\n")

  cat("✓ All 36 samples now accessible\n")
  cat("✓ GSM2453819_CV06.P10 recovered\n")
  cat("✓ Backup saved to:", backup_file, "\n\n")

} else {
  cat("✗ Fix failed\n")
  cat("  Expected 36 samples, got:", ncol(exprs_verify), "\n")
  cat("  Has GSM2453819:", "GSM2453819_CV06.P10" %in% colnames(exprs_verify), "\n\n")
  file.remove(temp_file)
}
