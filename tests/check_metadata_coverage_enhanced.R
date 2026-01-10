#!/usr/bin/env Rscript

# Check if mapped datasets have metadata in samples.csv
# AND verify that column names match metadata rows

cat("\n=== Checking Metadata Coverage for Mapped Datasets (Enhanced) ===\n\n")

# Read samples.csv
samples <- read.csv("data/phenodata/samples.csv", stringsAsFactors = FALSE, na.strings = c("NA", ""))

# Get mapped files
mapped_files <- list.files("data/mapped", pattern = "\\.tsv$", full.names = TRUE)

cat("Found", length(mapped_files), "mapped files\n\n")

# Extract dataset IDs from filenames
extract_dataset_id <- function(filename) {
  basename_file <- basename(filename)
  if (grepl("^GSE", basename_file)) {
    return(sub("_.*", "", basename_file))
  } else if (grepl("^E-GEOD-", basename_file)) {
    egeod <- sub("_.*", "", basename_file)
    gse <- sub("E-GEOD-", "GSE", egeod)
    return(list(egeod = egeod, gse = gse))
  }
  return(NULL)
}

results <- data.frame(
  mapped_file = character(),
  dataset_id = character(),
  samples_in_metadata = integer(),
  samples_in_file = integer(),
  has_healthy = logical(),
  has_diagnosis = logical(),
  has_gestational_category = logical(),
  has_biological_specimen = logical(),
  all_columns_match = logical(),
  pct_columns_match = numeric(),
  stringsAsFactors = FALSE
)

for (file in mapped_files) {
  basename_file <- basename(file)
  ids <- extract_dataset_id(file)
  
  if (is.list(ids)) {
    dataset_matches <- samples[samples$accession == ids$egeod | 
                                samples$secondaryaccession == ids$gse, ]
    dataset_id <- ids$gse
  } else if (!is.null(ids)) {
    dataset_matches <- samples[samples$secondaryaccession == ids, ]
    dataset_id <- ids
  } else {
    next
  }
  
  n_samples <- nrow(dataset_matches)
  has_healthy <- any(dataset_matches$Diagnosis == "Healthy", na.rm = TRUE)
  has_diagnosis <- any(!is.na(dataset_matches$Diagnosis) & 
                       dataset_matches$Diagnosis != "" & 
                       dataset_matches$Diagnosis != "NA", na.rm = TRUE)
  has_gest_cat <- any(!is.na(dataset_matches$Gestational.Age.Category) & 
                      dataset_matches$Gestational.Age.Category != "" & 
                      dataset_matches$Gestational.Age.Category != "Unknown Gestational Category" &
                      dataset_matches$Gestational.Age.Category != "NA", na.rm = TRUE)
  has_bio_spec <- any(!is.na(dataset_matches$Biological.Specimen) & 
                      dataset_matches$Biological.Specimen != "" & 
                      dataset_matches$Biological.Specimen != "NA", na.rm = TRUE)
  
  # NEW: Check if column names in mapped file match metadata
  mapped_data <- tryCatch({
    read.table(file, header = TRUE, sep = "\t", nrows = 1, check.names = FALSE)
  }, error = function(e) {
    NULL
  })
  
  all_columns_match <- NA
  pct_columns_match <- NA
  n_samples_in_file <- NA
  
  if (!is.null(mapped_data)) {
    # Get column names from mapped file (excluding first column which is gene names)
    file_columns <- colnames(mapped_data)[-1]
    n_samples_in_file <- length(file_columns)
    
    # Get expected column names from metadata
    # The arraydatafile_exprscolumnnames field should match the column names
    metadata_columns <- dataset_matches$arraydatafile_exprscolumnnames
    
    # Check how many file columns have matching metadata
    if (length(metadata_columns) > 0 && length(file_columns) > 0) {
      matches <- file_columns %in% metadata_columns
      pct_columns_match <- sum(matches) / length(file_columns) * 100
      all_columns_match <- all(matches)
    }
  }
  
  results <- rbind(results, data.frame(
    mapped_file = basename_file,
    dataset_id = dataset_id,
    samples_in_metadata = n_samples,
    samples_in_file = n_samples_in_file,
    has_healthy = has_healthy,
    has_diagnosis = has_diagnosis,
    has_gestational_category = has_gest_cat,
    has_biological_specimen = has_bio_spec,
    all_columns_match = all_columns_match,
    pct_columns_match = pct_columns_match,
    stringsAsFactors = FALSE
  ))
  
  # Print details
  cat("File:", basename_file, "\n")
  cat("  Dataset ID:", dataset_id, "\n")
  cat("  Samples in metadata:", n_samples, "\n")
  cat("  Samples in file:", n_samples_in_file, "\n")
  cat("  Has Diagnosis:", has_diagnosis, "\n")
  cat("  Has Healthy samples:", has_healthy, "\n")
  cat("  Has Gestational Category:", has_gest_cat, "\n")
  cat("  Has Biological Specimen:", has_bio_spec, "\n")
  
  if (!is.na(all_columns_match)) {
    cat("  Column names match:", ifelse(all_columns_match, "YES", "NO"), 
        sprintf("(%.1f%%)", pct_columns_match), "\n")
  }
  
  if (n_samples == 0) {
    cat("  ⚠️  WARNING: No metadata found!\n")
  } else if (!has_diagnosis || !has_gest_cat || !has_bio_spec) {
    cat("  ⚠️  WARNING: Missing key metadata fields!\n")
  } else if (!is.na(all_columns_match) && !all_columns_match) {
    cat("  ⚠️  WARNING: Column names don't all match metadata!\n")
  } else {
    cat("  ✓ Complete metadata with matching columns\n")
  }
  cat("\n")
}

cat("\n=== Summary ===\n")
cat("Total mapped files:", nrow(results), "\n")
cat("Files with metadata:", sum(results$samples_in_metadata > 0, na.rm = TRUE), "\n")
cat("Files with complete metadata:", sum(results$has_diagnosis & 
                                          results$has_gestational_category & 
                                          results$has_biological_specimen, na.rm = TRUE), "\n")
cat("Files with all columns matching:", sum(results$all_columns_match, na.rm = TRUE), "\n")
cat("Files missing metadata:", sum(results$samples_in_metadata == 0, na.rm = TRUE), "\n")

# Show files with issues
incomplete <- results[results$samples_in_metadata == 0 | 
                     !results$has_diagnosis | 
                     !results$has_gestational_category | 
                     !results$has_biological_specimen |
                     (!is.na(results$all_columns_match) & !results$all_columns_match), ]

if (nrow(incomplete) > 0) {
  cat("\n⚠️  Files with missing/incomplete metadata or column mismatches:\n")
  for (i in 1:nrow(incomplete)) {
    issues <- c()
    if (incomplete$samples_in_metadata[i] == 0) issues <- c(issues, "no metadata")
    if (!incomplete$has_diagnosis[i]) issues <- c(issues, "missing diagnosis")
    if (!incomplete$has_gestational_category[i]) issues <- c(issues, "missing gest. category")
    if (!incomplete$has_biological_specimen[i]) issues <- c(issues, "missing bio. specimen")
    if (!is.na(incomplete$all_columns_match[i]) && !incomplete$all_columns_match[i]) {
      issues <- c(issues, sprintf("column mismatch (%.1f%%)", incomplete$pct_columns_match[i]))
    }
    cat(sprintf("  - %s: %s\n", incomplete$dataset_id[i], paste(issues, collapse = ", ")))
  }
}

# Save detailed results
write.csv(results, "metadata_coverage_detailed.csv", row.names = FALSE)
cat("\n✓ Detailed results saved to: metadata_coverage_detailed.csv\n")

cat("\n✓ Check complete!\n\n")
