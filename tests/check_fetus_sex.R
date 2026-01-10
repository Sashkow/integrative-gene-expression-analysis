#!/usr/bin/env Rscript

samples <- read.csv("data/phenodata/samples.csv", stringsAsFactors = FALSE)
baseline_datasets <- c("GSE122214", "GSE22490", "GSE37901", "GSE9984")
baseline_samples <- samples[samples$secondaryaccession %in% baseline_datasets, ]

cat("Total baseline samples:", nrow(baseline_samples), "\n\n")
cat("Fetus.Sex distribution:\n")
print(table(baseline_samples$Fetus.Sex, useNA="always"))

cat("\nValid Fetus.Sex (not blank, not underscore):\n")
valid_count <- sum(baseline_samples$Fetus.Sex != "" & baseline_samples$Fetus.Sex != "_" & !is.na(baseline_samples$Fetus.Sex))
cat("Valid:", valid_count, "\n")
cat("Blank/underscore/NA:", nrow(baseline_samples) - valid_count, "\n")

if (valid_count > 0) {
  cat("\nValid values:\n")
  print(table(baseline_samples$Fetus.Sex[baseline_samples$Fetus.Sex != "" & baseline_samples$Fetus.Sex != "_"]))
}
