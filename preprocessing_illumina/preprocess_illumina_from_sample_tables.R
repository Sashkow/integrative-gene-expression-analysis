#!/usr/bin/env Rscript

#' Illumina Microarray Preprocessing from Sample Tables
#'
#' Processes Illumina sample table files (GSM*_sample_table.txt),
#' performs normalization, maps probes to genes using max-mean strategy,
#' and outputs expression matrix with ENTREZID identifiers
#'
#' @author Integrative Gene Expression Analysis Pipeline
#' @date 2025-10-09

cat("\n=== Illumina Microarray Preprocessing ===\n\n")

# =====================================================================
# 1. LOAD LIBRARIES
# =====================================================================

library(AnnotationDbi)           # Gene annotation
library(illuminaHumanv3.db)      # Illumina HumanRef-8 v3 annotation
library(limma)                   # General microarray tools

cat("✓ Libraries loaded\n\n")

# =====================================================================
# 2. CONFIGURATION
# =====================================================================

# Dataset to process
ACCESSION <- "E-GEOD-35574"

# Paths
raw_dir <- "/home/shivers/a/r/article-microarrays/raws/illumina/E-GEOD-35574"
output_dir <- "data/preprocessed_illumina"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Annotation database
ANNOTATION_DB <- illuminaHumanv3.db

cat("Configuration:\n")
cat("  Accession:", ACCESSION, "\n")
cat("  Raw data dir:", raw_dir, "\n")
cat("  Output dir:", output_dir, "\n\n")

# =====================================================================
# 3. LOAD RAW DATA FROM SAMPLE TABLES
# =====================================================================

cat("Step 1: Loading Illumina sample tables...\n")

# Find all sample table files
sample_files <- list.files(raw_dir, pattern = "GSM.*_sample_table\\.txt$",
                          full.names = TRUE)

if (length(sample_files) == 0) {
  stop("No sample table files found in ", raw_dir)
}

cat("  Found", length(sample_files), "sample files\n")

# Read first file to get probe IDs
first_file <- read.table(sample_files[1], header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)
probe_ids <- first_file[, 1]  # First column is probe ID

# Create expression matrix
exprs_raw <- matrix(NA, nrow = length(probe_ids), ncol = length(sample_files))
rownames(exprs_raw) <- probe_ids
colnames(exprs_raw) <- gsub("_sample_table\\.txt", "", basename(sample_files))

# Read expression values from all files
for (i in seq_along(sample_files)) {
  sample_data <- read.table(sample_files[i], header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE)
  # Second column is usually VALUE or expression
  exprs_raw[, i] <- sample_data[, 2]

  if (i %% 10 == 0) {
    cat("    Loaded", i, "/", length(sample_files), "samples\n")
  }
}

cat("  ✓ Raw data dimensions:", nrow(exprs_raw), "probes x",
    ncol(exprs_raw), "samples\n")

# =====================================================================
# 4. PREPROCESSING (Normalization)
# =====================================================================

cat("\nStep 2: Preprocessing (quantile normalization)...\n")

# Quantile normalization
exprs_norm <- normalizeBetweenArrays(exprs_raw, method = "quantile")

cat("  ✓ Quantile normalization complete\n")
cat("  Normalized data dimensions:", nrow(exprs_norm), "probes x",
    ncol(exprs_norm), "samples\n")

# =====================================================================
# 5. FILTER LOW-EXPRESSED PROBES (OPTIONAL)
# =====================================================================

# Optional: filter low-expressed probes
# Skipping filtering to match original preprocessing
exprs_filtered <- exprs_norm

cat("\nStep 3: No filtering applied (using all normalized probes)\n")
cat("  Probes:", nrow(exprs_filtered), "\n")

# =====================================================================
# 6. MAP PROBES TO GENES (ENTREZID)
# =====================================================================

cat("\nStep 4: Mapping probes to ENTREZID...\n")

# Get probe IDs
probe_ids <- rownames(exprs_filtered)

# Map to ENTREZID using annotation database
annotation <- tryCatch({
  AnnotationDbi::select(ANNOTATION_DB,
                        keys = probe_ids,
                        columns = "ENTREZID",
                        keytype = "PROBEID")
}, error = function(e) {
  # Try with ACCNUM if PROBEID doesn't work
  AnnotationDbi::select(ANNOTATION_DB,
                        keys = probe_ids,
                        columns = "ENTREZID",
                        keytype = "ACCNUM")
})

cat("  Initial mapping:", nrow(annotation), "probe-gene pairs\n")

# Remove probes without ENTREZID
annotation <- annotation[!is.na(annotation$ENTREZID), ]
cat("  After removing NA:", nrow(annotation), "pairs\n")

# Get the correct keytype column name
key_col <- colnames(annotation)[1]  # First column is the probe ID column

# Filter expression matrix to only mapped probes
exprs_mapped <- exprs_filtered[annotation[[key_col]], ]

cat("  Expression matrix:", nrow(exprs_mapped), "probes\n")

# =====================================================================
# 7. COLLAPSE MULTIPLE PROBES PER GENE (MAX MEAN STRATEGY)
# =====================================================================

cat("\nStep 5: Collapsing multiple probes per gene (max-mean strategy)...\n")

# Calculate mean expression per probe
probe_means <- rowMeans(exprs_mapped, na.rm = TRUE)

# Create data frame with probe, gene, mean expression
probe_gene_df <- data.frame(
  PROBEID = annotation[[key_col]],
  ENTREZID = annotation$ENTREZID,
  mean_expr = probe_means[annotation[[key_col]]],
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
# 9. SUMMARY
# =====================================================================

cat("\n=== Preprocessing Summary ===\n")
cat("Dataset:", ACCESSION, "\n")
cat("Raw probes:", nrow(exprs_raw), "\n")
cat("After normalization:", nrow(exprs_norm), "\n")
cat("After filtering:", nrow(exprs_filtered), "\n")
cat("After mapping:", nrow(exprs_mapped), "\n")
cat("Final unique genes:", nrow(exprs_final), "\n")
cat("Samples:", ncol(exprs_final), "\n")
cat("\nOutput file:", output_file, "\n")

cat("\n✓ Preprocessing complete!\n\n")
