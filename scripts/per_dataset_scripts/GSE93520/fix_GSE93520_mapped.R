#!/usr/bin/env Rscript

#' Fix GSE93520 mapped file to have proper TSV format
#'
#' @author Expression Integration Pipeline
#' @date 2025-10-26

cat("\n=== Fixing GSE93520 Mapped File ===\n\n")

# Define paths
original_file <- "/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis/data/mapped/GSE93520_matrix_no_filtering.tsv"
backup_file <- "/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis/data/mapped/GSE93520_matrix_no_filtering.tsv.backup"
fixed_file <- "/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis/data/mapped/GSE93520_matrix_no_filtering.tsv"

# ===================================================================
# STEP 1: BACKUP ORIGINAL FILE
# ===================================================================

cat("Step 1: Backing up original file...\n")
file.copy(original_file, backup_file, overwrite = TRUE)
cat("  ✓ Backup created:", backup_file, "\n\n")

# ===================================================================
# STEP 2: LOAD THE DATA PROPERLY
# ===================================================================

cat("Step 2: Loading data from original file...\n")

# Load WITHOUT row.names to get all 36 samples
exprs_data <- read.table(original_file, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)

cat("  Original dimensions:", nrow(exprs_data), "genes ×", ncol(exprs_data), "samples\n")
cat("  First 3 sample names:", paste(head(colnames(exprs_data), 3), collapse = ", "), "\n")
cat("  Last 3 sample names:", paste(tail(colnames(exprs_data), 3), collapse = ", "), "\n")
cat("  Has GSM2453819_CV06.P10:", "GSM2453819_CV06.P10" %in% colnames(exprs_data), "\n\n")

# Check that we have all 36 samples
if (ncol(exprs_data) != 36) {
  stop("Expected 36 samples but got ", ncol(exprs_data))
}

# ===================================================================
# STEP 3: SET ROW NAMES
# ===================================================================

cat("Step 3: Setting row names (gene IDs)...\n")

# The row names should be the gene/probe IDs from the GEO data
# For now, we'll use sequential numbers as they appear in the data
# (The original file has genes as 1, 2, 3, etc.)

# Check if we need to extract row names from somewhere
# Read first column of original data
first_col <- read.table(original_file, header = FALSE, sep = "\t",
                        stringsAsFactors = FALSE, colClasses = "character")[, 1]

# Remove header
gene_ids <- first_col[-1]

cat("  Number of gene IDs:", length(gene_ids), "\n")
cat("  First 3 gene IDs:", paste(head(gene_ids, 3), collapse = ", "), "\n\n")

# Set row names
rownames(exprs_data) <- gene_ids

# ===================================================================
# STEP 4: SAVE WITH PROPER FORMAT
# ===================================================================

cat("Step 4: Saving with proper TSV format...\n")

# Write with col.names=NA to add empty first field for row names column
write.table(exprs_data, fixed_file, sep = "\t", quote = FALSE,
            row.names = TRUE, col.names = NA)

cat("  ✓ Saved fixed file:", fixed_file, "\n\n")

# ===================================================================
# STEP 5: VERIFY THE FIX
# ===================================================================

cat("Step 5: Verifying the fixed file...\n")

# Load with row.names=1 (the way it will be used in the pipeline)
exprs_verify <- read.table(fixed_file, header = TRUE, sep = "\t", row.names = 1,
                           stringsAsFactors = FALSE, check.names = FALSE)

cat("  Verified dimensions:", nrow(exprs_verify), "genes ×", ncol(exprs_verify), "samples\n")
cat("  First 3 sample names:", paste(head(colnames(exprs_verify), 3), collapse = ", "), "\n")
cat("  Last 3 sample names:", paste(tail(colnames(exprs_verify), 3), collapse = ", "), "\n")
cat("  Has GSM2453819_CV06.P10:", "GSM2453819_CV06.P10" %in% colnames(exprs_verify), "\n\n")

# Check header format
cat("Step 6: Checking header format...\n")
header_line <- readLines(fixed_file, n = 1)
header_fields <- strsplit(header_line, "\t")[[1]]
cat("  Number of fields in header:", length(header_fields), "\n")
cat("  First field (should be empty or column name):",
    if (header_fields[1] == "") "EMPTY (correct!)" else paste0("'", header_fields[1], "'"), "\n")
cat("  Second field (first sample):", header_fields[2], "\n")
cat("  Third field (second sample):", header_fields[3], "\n\n")

# ===================================================================
# FINAL VERIFICATION
# ===================================================================

cat("=================================================================\n")
cat("                    FIX COMPLETE                                 \n")
cat("=================================================================\n\n")

if (ncol(exprs_verify) == 36 && "GSM2453819_CV06.P10" %in% colnames(exprs_verify)) {
  cat("✓ SUCCESS! All 36 samples are now accessible, including GSM2453819_CV06.P10\n\n")
  cat("Summary:\n")
  cat("  Original file (broken): 35 samples accessible\n")
  cat("  Fixed file (now):       36 samples accessible\n")
  cat("  GSM2453819_CV06.P10:    ✓ RECOVERED\n\n")
  cat("Backup of original file:", backup_file, "\n")
  cat("Fixed file:", fixed_file, "\n\n")
} else {
  cat("✗ ERROR: Fix did not work as expected\n")
  cat("  Expected 36 samples, got:", ncol(exprs_verify), "\n")
  cat("  Has GSM2453819:", "GSM2453819_CV06.P10" %in% colnames(exprs_verify), "\n")
}
