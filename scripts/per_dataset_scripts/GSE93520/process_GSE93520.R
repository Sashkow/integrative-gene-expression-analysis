#!/usr/bin/env Rscript

#' Process downloaded GSE93520 series matrix
#'
#' @author Expression Integration Pipeline
#' @date 2025-10-26

cat("\n=== Processing GSE93520 Series Matrix ===\n\n")

# Load required libraries
suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

# Set working directory
setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis/GSE93520_redownload")

# ===================================================================
# STEP 1: LOAD DOWNLOADED SERIES MATRIX
# ===================================================================

cat("Step 1: Loading series matrix file...\n")

matrix_file <- "GSE93520_series_matrix.txt.gz"

if (!file.exists(matrix_file)) {
  stop("Series matrix file not found: ", matrix_file)
}

# Parse the series matrix file
gse <- getGEO(filename = matrix_file, getGPL = FALSE)

cat("✓ Successfully loaded GSE93520 series matrix\n\n")

# ===================================================================
# STEP 2: EXTRACT PHENODATA
# ===================================================================

cat("Step 2: Extracting phenodata...\n")

# Get phenotype data
pdata <- pData(gse)

cat("  Phenodata dimensions:", nrow(pdata), "samples ×", ncol(pdata), "variables\n")
cat("  Sample IDs column:", "geo_accession\n")
cat("  First few sample IDs:\n")
print(head(pdata$geo_accession, 10))
cat("  Last few sample IDs:\n")
print(tail(pdata$geo_accession, 10))

# Save phenodata
pdata_file <- "GSE93520_phenodata.csv"
write.csv(pdata, pdata_file, row.names = FALSE, quote = TRUE)
cat("\n  ✓ Saved phenodata to:", pdata_file, "\n\n")

# ===================================================================
# STEP 3: EXTRACT EXPRESSION DATA
# ===================================================================

cat("Step 3: Extracting expression data...\n")

# Get expression data
exprs_data <- exprs(gse)

cat("  Expression data dimensions:", nrow(exprs_data), "probes/genes ×",
    ncol(exprs_data), "samples\n")
cat("  First few sample IDs:\n")
print(head(colnames(exprs_data), 10))
cat("  Last few sample IDs:\n")
print(tail(colnames(exprs_data), 10))

# Save expression data
exprs_file <- "GSE93520_expression.tsv"
write.table(exprs_data, exprs_file, sep = "\t", quote = FALSE,
            row.names = TRUE, col.names = NA)
cat("\n  ✓ Saved expression data to:", exprs_file, "\n\n")

# ===================================================================
# STEP 4: CHECK SAMPLE CONSISTENCY
# ===================================================================

cat("Step 4: Checking sample consistency...\n")

# Get sample IDs from phenodata
pheno_samples <- pdata$geo_accession

# Get sample IDs from expression data
exprs_samples <- colnames(exprs_data)

cat("  Samples in phenodata:", length(pheno_samples), "\n")
cat("  Samples in expression data:", length(exprs_samples), "\n")

# Check if they match
if (length(pheno_samples) == length(exprs_samples)) {
  cat("  ✓ Sample counts MATCH!\n\n")
} else {
  cat("  ✗ WARNING: Sample counts DO NOT match!\n\n")
}

# Check if all samples are present in both
samples_in_both <- intersect(pheno_samples, exprs_samples)
samples_pheno_only <- setdiff(pheno_samples, exprs_samples)
samples_exprs_only <- setdiff(exprs_samples, pheno_samples)

cat("  Samples in both phenodata and expression:", length(samples_in_both), "\n")

if (length(samples_pheno_only) > 0) {
  cat("  ✗ Samples in phenodata but NOT in expression (",
      length(samples_pheno_only), "):\n", sep = "")
  print(samples_pheno_only)
} else {
  cat("  ✓ All phenodata samples are in expression data\n")
}

if (length(samples_exprs_only) > 0) {
  cat("  ✗ Samples in expression but NOT in phenodata (",
      length(samples_exprs_only), "):\n", sep = "")
  print(samples_exprs_only)
} else {
  cat("  ✓ All expression samples are in phenodata\n")
}

cat("\n")

# ===================================================================
# STEP 5: CHECK FOR GSM2453819
# ===================================================================

cat("Step 5: Checking for GSM2453819 specifically...\n")

# Check if GSM2453819 is present
gsm2453819_in_pheno <- "GSM2453819" %in% pheno_samples
gsm2453819_in_exprs <- "GSM2453819" %in% exprs_samples

cat("  GSM2453819 in phenodata:", gsm2453819_in_pheno, "\n")
cat("  GSM2453819 in expression data:", gsm2453819_in_exprs, "\n\n")

if (gsm2453819_in_pheno && gsm2453819_in_exprs) {
  cat("  ✓ GSM2453819 is present in BOTH datasets\n")
} else if (gsm2453819_in_pheno && !gsm2453819_in_exprs) {
  cat("  ✗ GSM2453819 is in phenodata but NOT in expression data\n")
} else if (!gsm2453819_in_pheno && gsm2453819_in_exprs) {
  cat("  ✗ GSM2453819 is in expression data but NOT in phenodata\n")
} else {
  cat("  ✗ GSM2453819 is NOT in either dataset\n")
}

# If GSM2453819 is in phenodata, show its details
if (gsm2453819_in_pheno) {
  cat("\n  GSM2453819 phenodata details:\n")
  gsm_idx <- which(pdata$geo_accession == "GSM2453819")
  if (length(gsm_idx) > 0) {
    cat("    Title:", as.character(pdata$title[gsm_idx]), "\n")
    if ("source_name_ch1" %in% colnames(pdata)) {
      cat("    Source:", as.character(pdata$source_name_ch1[gsm_idx]), "\n")
    }
  }
}

cat("\n")

# ===================================================================
# STEP 6: DISPLAY ALL SAMPLE IDS
# ===================================================================

cat("Step 6: Complete list of samples in GEO data:\n")
cat("  All", length(pheno_samples), "samples:\n")
for (i in 1:length(pheno_samples)) {
  cat("    ", i, ":", pheno_samples[i], "\n", sep = "")
}
cat("\n")

# ===================================================================
# STEP 7: SAVE SUMMARY REPORT
# ===================================================================

cat("Step 7: Saving summary report...\n")

summary_file <- "GSE93520_download_summary.txt"
sink(summary_file)

cat("GSE93520 Download Summary\n")
cat(rep("=", 60), "\n\n", sep = "")
cat("Processing date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("GEO accession: GSE93520\n\n")

cat("Phenodata:\n")
cat("  File:", pdata_file, "\n")
cat("  Samples:", length(pheno_samples), "\n")
cat("  Variables:", ncol(pdata), "\n\n")

cat("Expression data:\n")
cat("  File:", exprs_file, "\n")
cat("  Probes/genes:", nrow(exprs_data), "\n")
cat("  Samples:", length(exprs_samples), "\n\n")

cat("Sample consistency:\n")
cat("  Samples in both:", length(samples_in_both), "\n")
cat("  Samples in phenodata only:", length(samples_pheno_only), "\n")
cat("  Samples in expression only:", length(samples_exprs_only), "\n\n")

if (length(samples_pheno_only) > 0) {
  cat("Samples in phenodata but NOT in expression:\n")
  for (s in samples_pheno_only) {
    cat("  -", s, "\n")
  }
  cat("\n")
}

if (length(samples_exprs_only) > 0) {
  cat("Samples in expression but NOT in phenodata:\n")
  for (s in samples_exprs_only) {
    cat("  -", s, "\n")
  }
  cat("\n")
}

cat("GSM2453819 check:\n")
cat("  In phenodata:", gsm2453819_in_pheno, "\n")
cat("  In expression:", gsm2453819_in_exprs, "\n\n")

cat("All sample IDs (", length(pheno_samples), "):\n", sep = "")
for (i in 1:length(pheno_samples)) {
  cat("  ", i, ". ", pheno_samples[i], "\n", sep = "")
}
cat("\n")

if (length(pheno_samples) == length(exprs_samples) &&
    length(samples_in_both) == length(pheno_samples)) {
  cat("CONCLUSION: All samples have both phenodata and expression data ✓\n")
} else {
  cat("CONCLUSION: Some samples are missing from one or both datasets ✗\n")
}

sink()

cat("✓ Summary saved to:", summary_file, "\n\n")

# ===================================================================
# FINAL SUMMARY
# ===================================================================

cat("=================================================================\n")
cat("                  PROCESSING COMPLETE                            \n")
cat("=================================================================\n\n")

cat("Dataset: GSE93520\n")
cat("Phenodata samples:", length(pheno_samples), "\n")
cat("Expression samples:", length(exprs_samples), "\n")
cat("Samples with both:", length(samples_in_both), "\n\n")

if (length(pheno_samples) == length(exprs_samples) &&
    length(samples_in_both) == length(pheno_samples)) {
  cat("✓ All samples have COMPLETE data!\n")
} else {
  cat("✗ Some samples are missing data\n")
  cat("  Missing from expression:", length(samples_pheno_only), "\n")
  cat("  Missing from phenodata:", length(samples_exprs_only), "\n")
}

cat("\nFiles saved in:", getwd(), "\n\n")
