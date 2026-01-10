#!/usr/bin/env Rscript

#' Add Combined.Fetus.Sex column to samples.csv
#'
#' This script adds a Combined.Fetus.Sex column that uses:
#' 1. Fetus.Sex if not empty or "_"
#' 2. Otherwise Estimated.Fetus.Sex
#' 3. Otherwise blank
#'
#' @author Expression Integration Pipeline
#' @date 2025-10-26

cat("\n=== Adding Combined.Fetus.Sex Column ===\n\n")

# Load data
samples_file <- "data/phenodata/samples.csv"
backup_file <- "data/phenodata/samples_backup_before_combined_sex.csv"

if (!file.exists(samples_file)) {
  stop("File not found: ", samples_file)
}

cat("Reading", samples_file, "...\n")
samples <- read.csv(samples_file, stringsAsFactors = FALSE)

cat("Original data:\n")
cat("  Rows:", nrow(samples), "\n")
cat("  Columns:", ncol(samples), "\n")

# Check if Combined.Fetus.Sex already exists
if ("Combined.Fetus.Sex" %in% colnames(samples)) {
  cat("\n⚠ Combined.Fetus.Sex column already exists. Recreating it...\n")
  samples$Combined.Fetus.Sex <- NULL
}

# Create backup
cat("\nCreating backup:", backup_file, "\n")
write.csv(samples, backup_file, row.names = FALSE)

# Check which columns exist
has_fetus_sex <- "Fetus.Sex" %in% colnames(samples)
has_estimated_sex <- "Estimated.Fetus.Sex" %in% colnames(samples)

cat("\nAvailable columns:\n")
cat("  Fetus.Sex:", has_fetus_sex, "\n")
cat("  Estimated.Fetus.Sex:", has_estimated_sex, "\n\n")

if (!has_fetus_sex && !has_estimated_sex) {
  stop("Neither Fetus.Sex nor Estimated.Fetus.Sex columns found!")
}

# Initialize Combined.Fetus.Sex with empty strings
samples$Combined.Fetus.Sex <- ""

# Logic: Use Fetus.Sex if valid, else Estimated.Fetus.Sex, else blank
for (i in 1:nrow(samples)) {

  # Check Fetus.Sex first
  if (has_fetus_sex) {
    fetus_sex <- samples$Fetus.Sex[i]
    # Use it if not NA, not empty, and not "_"
    if (!is.na(fetus_sex) && fetus_sex != "" && fetus_sex != "_") {
      samples$Combined.Fetus.Sex[i] <- fetus_sex
      next
    }
  }

  # If Fetus.Sex not valid, try Estimated.Fetus.Sex
  if (has_estimated_sex) {
    estimated_sex <- samples$Estimated.Fetus.Sex[i]
    # Use it if not NA and not empty
    if (!is.na(estimated_sex) && estimated_sex != "") {
      samples$Combined.Fetus.Sex[i] <- estimated_sex
      next
    }
  }

  # Otherwise leave it blank (already initialized to "")
}

# Report statistics
cat("Combined.Fetus.Sex statistics:\n")
cat("\nValue counts:\n")
print(table(samples$Combined.Fetus.Sex, useNA = "always"))

# Count sources
from_fetus_sex <- 0
from_estimated_sex <- 0
blank <- 0

for (i in 1:nrow(samples)) {
  combined <- samples$Combined.Fetus.Sex[i]

  if (combined == "") {
    blank <- blank + 1
  } else {
    # Check if it came from Fetus.Sex
    if (has_fetus_sex) {
      fetus_sex <- samples$Fetus.Sex[i]
      if (!is.na(fetus_sex) && fetus_sex != "" && fetus_sex != "_" && fetus_sex == combined) {
        from_fetus_sex <- from_fetus_sex + 1
        next
      }
    }
    # Otherwise it came from Estimated.Fetus.Sex
    from_estimated_sex <- from_estimated_sex + 1
  }
}

cat("\nSource breakdown:\n")
cat("  From Fetus.Sex:", from_fetus_sex, "\n")
cat("  From Estimated.Fetus.Sex:", from_estimated_sex, "\n")
cat("  Blank:", blank, "\n\n")

# Save updated file
cat("Saving updated file to:", samples_file, "\n")
write.csv(samples, samples_file, row.names = FALSE)

cat("\n✓ Combined.Fetus.Sex column added successfully!\n")
cat("\nBackup saved to:", backup_file, "\n")
cat("Original file updated:", samples_file, "\n\n")
