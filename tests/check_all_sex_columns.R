#!/usr/bin/env Rscript

samples <- read.csv("data/phenodata/samples.csv", stringsAsFactors = FALSE)
baseline_datasets <- c("GSE122214", "GSE22490", "GSE37901", "GSE9984")
baseline_samples <- samples[samples$secondaryaccession %in% baseline_datasets, ]

cat("Total baseline samples:", nrow(baseline_samples), "\n\n")

# Check all sex-related columns
sex_cols <- c("Fetus.Sex", "Sex", "Fetus", "Estimated.Fetus.Sex",
              "estimated_sex", "Combined.Fetus.Sex")

for (col in sex_cols) {
  if (col %in% names(baseline_samples)) {
    cat("===", col, "===\n")
    cat("Unique values:\n")
    print(unique(baseline_samples[[col]]))
    cat("\nDistribution:\n")
    print(table(baseline_samples[[col]], useNA="always"))

    # Count non-blank, non-underscore values
    valid <- sum(baseline_samples[[col]] != "" &
                 baseline_samples[[col]] != "_" &
                 !is.na(baseline_samples[[col]]))
    cat("Valid (not blank/underscore/NA):", valid, "\n")

    if (valid > 0) {
      cat("Valid values breakdown:\n")
      valid_values <- baseline_samples[[col]][baseline_samples[[col]] != "" &
                                               baseline_samples[[col]] != "_" &
                                               !is.na(baseline_samples[[col]])]
      print(table(valid_values))
    }
    cat("\n")
  }
}
