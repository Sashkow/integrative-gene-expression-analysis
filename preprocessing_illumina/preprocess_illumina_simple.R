#!/usr/bin/env Rscript

#' Simple Illumina Microarray Preprocessing Pipeline
#'
#' Downloads raw Illumina BeadArray data from ArrayExpress/GEO,
#' performs standard preprocessing (background correction, normalization),
#' maps probes to genes using max-mean strategy,
#' and outputs expression matrix with ENTREZID identifiers
#'
#' @author Integrative Gene Expression Analysis Pipeline
#' @date 2025-10-07

cat("\n=== Illumina Microarray Preprocessing ===\n\n")

# =====================================================================
# 1. LOAD LIBRARIES
# =====================================================================

library(ArrayExpress)            # Download data from ArrayExpress
library(GEOquery)                # Download data from GEO
library(AnnotationDbi)           # Gene annotation
library(illuminaHumanv3.db)      # Illumina HumanRef-8 v3 annotation
library(illuminaHumanv4.db)      # Illumina HumanHT-12 v4 annotation
library(limma)                   # General microarray tools
# Note: lumi not required for this simplified pipeline

cat("✓ Libraries loaded\n\n")

# =====================================================================
# 2. CONFIGURATION
# =====================================================================

# Dataset to download
# Note: This script requires raw Illumina data (non-normalized text files)
# E-GEOD-55439 does not have raw files available, using a dataset with raw data
ACCESSION <- "GSE100051"      # ArrayExpress/GEO accession with raw data

# Paths
raw_dir <- "data/raw/raw_illumina"
output_dir <- "data/preprocessed_illumina"

# Create directories
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Annotation database (E-GEOD-27272 uses HumanRef-8 v3)
ANNOTATION_DB <- illuminaHumanv3.db  # Change based on your platform
# Options: illuminaHumanv1.db, illuminaHumanv2.db, illuminaHumanv3.db, illuminaHumanv4.db

cat("Configuration:\n")
cat("  Accession:", ACCESSION, "\n")
cat("  Raw data dir:", raw_dir, "\n")
cat("  Output dir:", output_dir, "\n\n")

# =====================================================================
# 3. DOWNLOAD RAW DATA
# =====================================================================

# Detect if accession is from GEO or ArrayExpress
is_geo <- grepl("^GSE\\d+$", ACCESSION)
is_arrayexpress <- grepl("^E-", ACCESSION)

if (is_geo) {
  cat("Step 1: Downloading raw data from GEO...\n")
} else if (is_arrayexpress) {
  cat("Step 1: Downloading raw data from ArrayExpress...\n")
} else {
  stop("Unknown accession format: ", ACCESSION,
       "\nSupported formats: GSE##### (GEO) or E-**** (ArrayExpress)")
}

dataset_dir <- file.path(raw_dir, ACCESSION)
dir.create(dataset_dir, showWarnings = FALSE, recursive = TRUE)

# Check if already downloaded
existing_files <- list.files(dataset_dir, pattern = "\\.txt(\\.gz)?$", full.names = TRUE)

if (length(existing_files) > 0) {
  cat("  Data already downloaded\n")
} else {
  if (is_geo) {
    # Download from GEO
    cat("  Downloading from GEO...\n")

    # Get supplementary files (contains raw data)
    tryCatch({
      getGEOSuppFiles(ACCESSION, baseDir = raw_dir, makeDirectory = TRUE)
      cat("  ✓ Downloaded supplementary files to:", dataset_dir, "\n")

      # Untar/unzip if needed
      supp_files <- list.files(dataset_dir, full.names = TRUE)
      for (file in supp_files) {
        if (grepl("\\.tar$", file)) {
          cat("  Extracting tar archive...\n")
          untar(file, exdir = dataset_dir)
        } else if (grepl("\\.gz$", file) && !grepl("\\.tar\\.gz$", file)) {
          cat("  Decompressing gzip file...\n")
          system(paste("gunzip -f", shQuote(file)))
        } else if (grepl("\\.tar\\.gz$", file)) {
          cat("  Extracting tar.gz archive...\n")
          untar(file, exdir = dataset_dir)
        }
      }
      cat("  ✓ Extraction complete\n")
    }, error = function(e) {
      cat("  Warning: Could not download supplementary files\n")
      cat("  Error:", conditionMessage(e), "\n")
      cat("  Trying to download series matrix instead...\n")

      # Fallback: try to get series matrix (normalized data)
      gse <- getGEO(ACCESSION, destdir = dataset_dir, GSEMatrix = TRUE)
      if (length(gse) > 0) {
        # Extract expression data from series matrix
        exprs_data <- exprs(gse[[1]])
        output_file <- file.path(dataset_dir, paste0(ACCESSION, "_series_matrix.txt"))
        write.table(exprs_data, output_file, sep = "\t", quote = FALSE,
                    row.names = TRUE, col.names = NA)
        cat("  ✓ Downloaded series matrix to:", output_file, "\n")
      }
    })
  } else if (is_arrayexpress) {
    # Download from ArrayExpress
    cat("  Downloading from ArrayExpress...\n")
    aeData <- getAE(
      ACCESSION,
      path = dataset_dir,
      type = "raw"
    )
    cat("  ✓ Downloaded to:", dataset_dir, "\n")
  }
}

# =====================================================================
# 4. LOAD RAW DATA
# =====================================================================

cat("\nStep 2: Loading raw Illumina data...\n")

# Find the raw data file (usually non-normalized or raw)
raw_files <- list.files(dataset_dir, pattern = "non.normalized|raw|_RAW\\.txt|series_matrix\\.txt",
                        ignore.case = TRUE, full.names = TRUE, recursive = TRUE)

if (length(raw_files) == 0) {
  # If no raw file found, look for any .txt file
  raw_files <- list.files(dataset_dir, pattern = "\\.txt$", full.names = TRUE, recursive = TRUE)
}

if (length(raw_files) == 0) {
  stop("No raw data files found in ", dataset_dir)
}

cat("  Found raw file:", basename(raw_files[1]), "\n")

# Read expression data
# Try different read methods for different file formats
exprs_raw <- tryCatch({
  read.table(raw_files[1], header = TRUE, sep = '\t',
             comment.char = "!", row.names = 1, check.names = FALSE)
}, error = function(e) {
  cat("  Trying alternative read method...\n")
  # Try without comment.char
  tryCatch({
    read.table(raw_files[1], header = TRUE, sep = '\t',
               row.names = 1, check.names = FALSE)
  }, error = function(e2) {
    # Try with header in first row
    read.delim(raw_files[1], header = TRUE, row.names = 1, check.names = FALSE)
  })
})

cat("  Raw data dimensions:", nrow(exprs_raw), "probes x", ncol(exprs_raw), "columns\n")

# =====================================================================
# 4.5. EXTRACT EXPRESSION VALUES
# =====================================================================

# For Illumina data, we may have multiple columns per sample (e.g., AVG_Signal, Detection)
# Extract only the expression columns
cat("  Extracting expression columns...\n")

# Check column names
col_names <- colnames(exprs_raw)

# Find expression columns (typically AVG_Signal, Signal, or similar)
expr_cols <- grep("AVG_Signal|Signal|^SAMPLE|^GSM", col_names, ignore.case = TRUE)

if (length(expr_cols) > 0) {
  # Filter to expression columns only
  exprs_subset <- exprs_raw[, expr_cols, drop = FALSE]

  # Check if all columns are numeric
  numeric_cols <- sapply(exprs_subset, is.numeric)

  if (!all(numeric_cols)) {
    cat("  Warning: Some expression columns are not numeric. Attempting to convert...\n")
    for (i in which(!numeric_cols)) {
      exprs_subset[, i] <- as.numeric(as.character(exprs_subset[, i]))
    }
  }

  exprs_raw <- exprs_subset
  cat("  Extracted", ncol(exprs_raw), "expression columns\n")
} else {
  # If no specific columns found, assume all are numeric
  # Convert all columns to numeric
  for (i in seq_len(ncol(exprs_raw))) {
    if (!is.numeric(exprs_raw[, i])) {
      exprs_raw[, i] <- as.numeric(as.character(exprs_raw[, i]))
    }
  }
}

cat("  Final data dimensions:", nrow(exprs_raw), "probes x", ncol(exprs_raw), "samples\n")


# =====================================================================
# 5. PREPROCESSING (Log transformation + Normalization)
# =====================================================================

cat("\nStep 3: Preprocessing (log transformation + normalization)...\n")

# Convert to matrix if needed
exprs_raw <- as.matrix(exprs_raw)

# Remove probes with NA/NaN values
cat("  Removing probes with NA/NaN values...\n")
valid_probes <- rowSums(is.na(exprs_raw) | is.nan(exprs_raw)) == 0
exprs_clean <- exprs_raw[valid_probes, ]
cat("  Removed", sum(!valid_probes), "probes with NA/NaN values\n")

# Log2 transformation
# Add small offset to avoid log(0)
cat("  Applying log2 transformation...\n")
min_positive <- min(exprs_clean[exprs_clean > 0], na.rm = TRUE)
offset <- min_positive / 10  # Use 10% of minimum positive value as offset

exprs_log <- log2(exprs_clean + offset)
cat("  ✓ Log2 transformation complete (offset =", round(offset, 4), ")\n")

# Quantile normalization
cat("  Applying quantile normalization...\n")
exprs_norm <- normalizeBetweenArrays(exprs_log, method = "quantile")

cat("  ✓ Quantile normalization complete\n")
cat("  Normalized data dimensions:", nrow(exprs_norm), "probes x", ncol(exprs_norm), "samples\n")

# =====================================================================
# 6. MAP PROBES TO GENES (ENTREZID)
# =====================================================================

cat("\nStep 4: Mapping probes to ENTREZID...\n")

# Get probe IDs
probe_ids <- rownames(exprs_norm)

# Map to ENTREZID using annotation database
annotation <- tryCatch({
  AnnotationDbi::select(ANNOTATION_DB,
                        keys = probe_ids,
                        columns = "ENTREZID",
                        keytype = "PROBEID")
}, error = function(e) {
  # If PROBEID doesn't work, try other keytypes
  AnnotationDbi::select(ANNOTATION_DB,
                        keys = probe_ids,
                        columns = "ENTREZID")
})

cat("  Initial mapping:", nrow(annotation), "probe-gene pairs\n")

# Remove probes without ENTREZID
annotation <- annotation[!is.na(annotation$ENTREZID), ]
cat("  After removing NA:", nrow(annotation), "pairs\n")

# Filter expression matrix to only mapped probes
exprs_mapped <- exprs_norm[annotation$PROBEID, ]

cat("  Expression matrix:", nrow(exprs_mapped), "probes\n")

# =====================================================================
# 7. COLLAPSE MULTIPLE PROBES PER GENE (MAX MEAN STRATEGY)
# =====================================================================

cat("\nStep 5: Collapsing multiple probes per gene (max-mean strategy)...\n")

# Calculate mean expression per probe
probe_means <- rowMeans(exprs_mapped, na.rm = TRUE)

# Create data frame with probe, gene, mean expression
probe_gene_df <- data.frame(
  PROBEID = annotation$PROBEID,
  ENTREZID = annotation$ENTREZID,
  mean_expr = probe_means[annotation$PROBEID],
  stringsAsFactors = FALSE
)

# For each gene, keep the probe with maximum mean expression
probe_gene_df <- probe_gene_df[order(probe_gene_df$ENTREZID,
                                      -probe_gene_df$mean_expr), ]
probe_gene_df <- probe_gene_df[!duplicated(probe_gene_df$ENTREZID), ]

cat("  Probes before collapsing:", nrow(exprs_mapped), "\n")
cat("  Unique genes after collapsing:", nrow(probe_gene_df), "\n")

# Select probes
exprs_final <- exprs_mapped[probe_gene_df$PROBEID, ]

# Replace probe IDs with ENTREZID
rownames(exprs_final) <- probe_gene_df$ENTREZID

cat("  ✓ Final expression matrix:", nrow(exprs_final), "genes x",
    ncol(exprs_final), "samples\n")

# =====================================================================
# 8. SAVE RESULTS
# =====================================================================

cat("\nStep 6: Saving results...\n")

output_file <- file.path(output_dir, paste0(ACCESSION, "_mapped_illumina.tsv"))

write.table(exprs_final, output_file, sep = "\t", quote = FALSE,
            row.names = TRUE, col.names = NA)

cat("  ✓ Saved expression matrix to:", output_file, "\n")

# Save probe-to-gene mapping
mapping_file <- file.path(output_dir, paste0(ACCESSION, "_probe_mapping.tsv"))
write.table(probe_gene_df, mapping_file, sep = "\t", quote = FALSE,
            row.names = FALSE)

cat("  ✓ Saved probe mapping to:", mapping_file, "\n")

# =====================================================================
# 10. SUMMARY
# =====================================================================

cat("\n=== Preprocessing Summary ===\n")
cat("Dataset:", ACCESSION, "\n")
cat("Raw probes:", nrow(exprs_raw), "\n")
cat("After log2 transformation:", nrow(exprs_log), "\n")
cat("After normalization:", nrow(exprs_norm), "\n")
cat("After mapping:", nrow(exprs_mapped), "\n")
cat("Final unique genes:", nrow(exprs_final), "\n")
cat("Samples:", ncol(exprs_final), "\n")
cat("\nOutput file:", output_file, "\n")

cat("\n✓ Preprocessing complete!\n\n")

# =====================================================================
# NOTES
# =====================================================================

# Annotation databases for different Illumina platforms:
# - illuminaHumanv1.db  -> HumanRef-8 v1
# - illuminaHumanv2.db  -> HumanRef-8 v2
# - illuminaHumanv3.db  -> HumanHT-12 v3
# - illuminaHumanv4.db  -> HumanHT-12 v4
# - illuminaHumanWGDASLv3.db -> Whole-genome DASL
# - illuminaHumanWGDASLv4.db -> Whole-genome DASL v4

# Max-mean strategy:
# When multiple probes map to the same gene, keep the probe with
# the highest mean expression across samples. This probe is likely
# most reliable and specific.

# =====================================================================
# CHECK ANNOTATION DATABASE VERSIONS
# =====================================================================

cat("\n=== Annotation Database Versions ===\n")
cat("illuminaHumanv3.db:\n")
print(packageDescription('illuminaHumanv3.db', fields = c('Version', 'Date', 'Built')))
cat("\nilluminaHumanv4.db:\n")
print(packageDescription('illuminaHumanv4.db', fields = c('Version', 'Date', 'Built')))
