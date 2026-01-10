#!/usr/bin/env Rscript

# Check if mapped datasets have metadata in samples.csv

cat("\n=== Checking Metadata Coverage for Mapped Datasets ===\n\n")

# Read samples.csv
samples <- read.csv("data/phenodata/samples.csv", stringsAsFactors = FALSE)

# Get mapped files
mapped_files <- list.files("data/mapped", pattern = "\\.tsv$", full.names = FALSE)

cat("Found", length(mapped_files), "mapped files\n\n")

# Extract dataset IDs from filenames
extract_dataset_id <- function(filename) {
  # Try to extract GSE or E-GEOD accession
  if (grepl("^GSE", filename)) {
    return(sub("_.*", "", filename))
  } else if (grepl("^E-GEOD-", filename)) {
    # Extract E-GEOD-##### and corresponding GSE#####
    egeod <- sub("_.*", "", filename)
    gse <- sub("E-GEOD-", "GSE", egeod)
    return(list(egeod = egeod, gse = gse))
  }
  return(NULL)
}

results <- data.frame(
  mapped_file = character(),
  dataset_id = character(),
  samples_in_metadata = integer(),
  has_healthy = logical(),
  has_diagnosis = logical(),
  has_gestational_category = logical(),
  has_biological_specimen = logical(),
  stringsAsFactors = FALSE
)

for (file in mapped_files) {
  ids <- extract_dataset_id(file)
  
  if (is.list(ids)) {
    # Check both E-GEOD and GSE versions
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
  
  results <- rbind(results, data.frame(
    mapped_file = file,
    dataset_id = dataset_id,
    samples_in_metadata = n_samples,
    has_healthy = has_healthy,
    has_diagnosis = has_diagnosis,
    has_gestational_category = has_gest_cat,
    has_biological_specimen = has_bio_spec,
    stringsAsFactors = FALSE
  ))
  
  # Print details
  cat("File:", file, "\n")
  cat("  Dataset ID:", dataset_id, "\n")
  cat("  Samples in metadata:", n_samples, "\n")
  cat("  Has Diagnosis:", has_diagnosis, "\n")
  cat("  Has Healthy samples:", has_healthy, "\n")
  cat("  Has Gestational Category:", has_gest_cat, "\n")
  cat("  Has Biological Specimen:", has_bio_spec, "\n")
  
  if (n_samples == 0) {
    cat("  ⚠️  WARNING: No metadata found!\n")
  } else if (!has_diagnosis || !has_gest_cat || !has_bio_spec) {
    cat("  ⚠️  WARNING: Missing key metadata fields!\n")
  } else {
    cat("  ✓ Complete metadata\n")
  }
  cat("\n")
}

cat("\n=== Summary ===\n")
cat("Total mapped files:", nrow(results), "\n")
cat("Files with metadata:", sum(results$samples_in_metadata > 0), "\n")
cat("Files with complete metadata:", sum(results$has_diagnosis & 
                                          results$has_gestational_category & 
                                          results$has_biological_specimen), "\n")
cat("Files missing metadata:", sum(results$samples_in_metadata == 0), "\n")

# Show files with missing or incomplete metadata
incomplete <- results[results$samples_in_metadata == 0 | 
                     !results$has_diagnosis | 
                     !results$has_gestational_category | 
                     !results$has_biological_specimen, ]

if (nrow(incomplete) > 0) {
  cat("\n⚠️  Files with missing/incomplete metadata:\n")
  for (i in 1:nrow(incomplete)) {
    cat("  -", incomplete$dataset_id[i], "\n")
  }
}

cat("\n✓ Check complete!\n\n")
