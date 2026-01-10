#!/usr/bin/env Rscript

cat("=== Checking Sex Data Availability for First Trimester Datasets ===\n\n")

# Read phenodata
phenodata <- read.csv("data/phenodata/samples.csv", stringsAsFactors = FALSE)

# Define datasets
baseline_datasets <- c("GSE122214", "GSE22490", "GSE37901", "GSE9984")
first_trim_datasets <- c("GSE55439", "GSE93520", "GSE28551", "GSE100051")
all_datasets <- c(baseline_datasets, first_trim_datasets)

cat("Total samples in phenodata:", nrow(phenodata), "\n\n")

# Function to check sex data for a dataset
check_sex_data <- function(dataset_id, pdata) {

  # Filter samples for this dataset
  dataset_samples <- pdata[pdata$secondaryaccession == dataset_id, ]

  if (nrow(dataset_samples) == 0) {
    return(list(
      dataset = dataset_id,
      total_samples = 0,
      has_fetus_sex = 0,
      has_sex = 0,
      has_estimated = 0,
      has_combined = 0,
      has_any_sex = 0,
      pct_any_sex = 0
    ))
  }

  # Count samples with sex data in different columns
  has_fetus_sex <- sum(!is.na(dataset_samples$Fetus.Sex) &
                       dataset_samples$Fetus.Sex != "" &
                       dataset_samples$Fetus.Sex != "NA", na.rm = TRUE)

  has_sex <- sum(!is.na(dataset_samples$Sex) &
                 dataset_samples$Sex != "" &
                 dataset_samples$Sex != "NA", na.rm = TRUE)

  has_estimated <- sum(!is.na(dataset_samples$Estimated.Fetus.Sex) &
                       dataset_samples$Estimated.Fetus.Sex != "" &
                       dataset_samples$Estimated.Fetus.Sex != "NA", na.rm = TRUE)

  has_estimated_sex <- sum(!is.na(dataset_samples$estimated_sex) &
                           dataset_samples$estimated_sex != "" &
                           dataset_samples$estimated_sex != "NA", na.rm = TRUE)

  has_combined <- sum(!is.na(dataset_samples$Combined.Fetus.Sex) &
                      dataset_samples$Combined.Fetus.Sex != "" &
                      dataset_samples$Combined.Fetus.Sex != "NA", na.rm = TRUE)

  # Any sex data at all (from any column)
  has_any <- (!is.na(dataset_samples$Fetus.Sex) & dataset_samples$Fetus.Sex != "" & dataset_samples$Fetus.Sex != "NA") |
             (!is.na(dataset_samples$Sex) & dataset_samples$Sex != "" & dataset_samples$Sex != "NA") |
             (!is.na(dataset_samples$Estimated.Fetus.Sex) & dataset_samples$Estimated.Fetus.Sex != "" & dataset_samples$Estimated.Fetus.Sex != "NA") |
             (!is.na(dataset_samples$estimated_sex) & dataset_samples$estimated_sex != "" & dataset_samples$estimated_sex != "NA") |
             (!is.na(dataset_samples$Combined.Fetus.Sex) & dataset_samples$Combined.Fetus.Sex != "" & dataset_samples$Combined.Fetus.Sex != "NA")

  has_any_sex <- sum(has_any, na.rm = TRUE)
  pct_any_sex <- round(100 * has_any_sex / nrow(dataset_samples), 1)

  return(list(
    dataset = dataset_id,
    total_samples = nrow(dataset_samples),
    has_fetus_sex = has_fetus_sex,
    has_sex = has_sex,
    has_estimated_fetus_sex = has_estimated,
    has_estimated_sex = has_estimated_sex,
    has_combined = has_combined,
    has_any_sex = has_any_sex,
    pct_any_sex = pct_any_sex
  ))
}

# Check each dataset
results <- list()
for (dataset in all_datasets) {
  results[[dataset]] <- check_sex_data(dataset, phenodata)
}

# Convert to data frame
results_df <- do.call(rbind, lapply(results, as.data.frame))

# Add category
results_df$category <- ifelse(results_df$dataset %in% baseline_datasets,
                               "Baseline",
                               "First Trimester")

# Reorder columns
results_df <- results_df[, c("category", "dataset", "total_samples",
                             "has_fetus_sex", "has_sex",
                             "has_estimated_fetus_sex", "has_estimated_sex",
                             "has_combined", "has_any_sex", "pct_any_sex")]

# Print results
cat("========================================================================\n")
cat("SEX DATA AVAILABILITY BY DATASET\n")
cat("========================================================================\n\n")

print(results_df, row.names = FALSE)

cat("\n========================================================================\n")
cat("SUMMARY\n")
cat("========================================================================\n\n")

# Summary by category
cat("BASELINE DATASETS:\n")
baseline_results <- results_df[results_df$category == "Baseline", ]
cat("  Total datasets:", nrow(baseline_results), "\n")
cat("  Total samples:", sum(baseline_results$total_samples), "\n")
cat("  Samples with any sex data:", sum(baseline_results$has_any_sex),
    "(", round(100 * sum(baseline_results$has_any_sex) / sum(baseline_results$total_samples), 1), "%)\n")
cat("  Datasets with 100% sex data:", sum(baseline_results$pct_any_sex == 100), "\n")
cat("  Datasets with some sex data (>0%):", sum(baseline_results$pct_any_sex > 0), "\n")
cat("  Datasets with NO sex data:", sum(baseline_results$pct_any_sex == 0), "\n\n")

cat("FIRST TRIMESTER DATASETS:\n")
first_results <- results_df[results_df$category == "First Trimester", ]
cat("  Total datasets:", nrow(first_results), "\n")
cat("  Total samples:", sum(first_results$total_samples), "\n")
cat("  Samples with any sex data:", sum(first_results$has_any_sex),
    "(", round(100 * sum(first_results$has_any_sex) / sum(first_results$total_samples), 1), "%)\n")
cat("  Datasets with 100% sex data:", sum(first_results$pct_any_sex == 100), "\n")
cat("  Datasets with some sex data (>0%):", sum(first_results$pct_any_sex > 0), "\n")
cat("  Datasets with NO sex data:", sum(first_results$pct_any_sex == 0), "\n\n")

cat("ALL DATASETS:\n")
cat("  Total datasets:", nrow(results_df), "\n")
cat("  Total samples:", sum(results_df$total_samples), "\n")
cat("  Samples with any sex data:", sum(results_df$has_any_sex),
    "(", round(100 * sum(results_df$has_any_sex) / sum(results_df$total_samples), 1), "%)\n")
cat("  Datasets with 100% sex data:", sum(results_df$pct_any_sex == 100), "\n")
cat("  Datasets with some sex data (>0%):", sum(results_df$pct_any_sex > 0), "\n")
cat("  Datasets with NO sex data:", sum(results_df$pct_any_sex == 0), "\n\n")

# Save results
write.csv(results_df, "output/sex_data_availability.csv", row.names = FALSE)
cat("✓ Results saved to: output/sex_data_availability.csv\n\n")
