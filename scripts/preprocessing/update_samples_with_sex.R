#!/usr/bin/env Rscript

#' Update samples.csv with GSE93520 sex estimation results
#'
#' @author Expression Integration Pipeline
#' @date 2025-10-26

cat("\n=== Updating samples.csv with GSE93520 Sex Estimation ===\n\n")

# Load sex estimation results
sex_results <- read.csv("output/sex_estimation/GSE93520_sex_estimation.csv",
                        stringsAsFactors = FALSE)

cat("Loaded sex estimation results:", nrow(sex_results), "samples\n")
cat("  Female:", sum(sex_results$Estimated.Fetus.Sex == "Female"), "\n")
cat("  Male:", sum(sex_results$Estimated.Fetus.Sex == "Male"), "\n\n")

# Load current phenodata
phenodata <- read.csv("data/phenodata/samples.csv", stringsAsFactors = FALSE)

cat("Loaded phenodata:", nrow(phenodata), "total samples\n")
cat("  GSE93520 samples:", sum(phenodata$secondaryaccession == "GSE93520"), "\n\n")

# Backup phenodata
backup_file <- "data/phenodata/samples.csv.backup_before_sex_update"
if (!file.exists(backup_file)) {
  file.copy("data/phenodata/samples.csv", backup_file)
  cat("Created backup:", backup_file, "\n\n")
}

# Update phenodata with sex estimates
cat("Updating phenodata...\n")
updated_count <- 0

for (i in 1:nrow(sex_results)) {
  sample_id <- sex_results$arraydatafile_exprscolumnnames[i]

  # Find matching row in phenodata
  idx <- which(phenodata$arraydatafile_exprscolumnnames == sample_id)

  if (length(idx) > 0) {
    # Update sex fields
    phenodata$Estimated.Fetus.Sex[idx] <- sex_results$Estimated.Fetus.Sex[i]
    phenodata$estimated_sex[idx] <- sex_results$estimated_sex[i]
    phenodata$etimated_sex_unimodality[idx] <- sex_results$etimated_sex_unimodality[i]
    phenodata$etimated_sex_unimodality_p[idx] <- sex_results$etimated_sex_unimodality_p[i]

    updated_count <- updated_count + 1
  } else {
    cat("  Warning: Sample not found in phenodata:", sample_id, "\n")
  }
}

cat("  Updated", updated_count, "samples\n\n")

# Save updated phenodata
write.csv(phenodata, "data/phenodata/samples.csv", row.names = FALSE, quote = TRUE)
cat("✓ Saved updated phenodata to: data/phenodata/samples.csv\n\n")

# Verify update
cat("Verification:\n")
gse93520_samples <- phenodata[phenodata$secondaryaccession == "GSE93520", ]
cat("  GSE93520 samples with Estimated.Fetus.Sex:",
    sum(!is.na(gse93520_samples$Estimated.Fetus.Sex) &
        gse93520_samples$Estimated.Fetus.Sex != ""), "\n")
cat("  GSE93520 Female:",
    sum(gse93520_samples$Estimated.Fetus.Sex == "Female", na.rm = TRUE), "\n")
cat("  GSE93520 Male:",
    sum(gse93520_samples$Estimated.Fetus.Sex == "Male", na.rm = TRUE), "\n\n")

cat("=================================================================\n")
cat("                    UPDATE COMPLETE                              \n")
cat("=================================================================\n\n")
