#!/usr/bin/env Rscript

#' Check for multiple expression files matching same dataset
#'
#' Detects potential merge conflicts when multiple files exist for the same
#' dataset ID (e.g., GSE100051.tsv and GSE100051_full.tsv would both be loaded
#' when requesting GSE100051, causing column name conflicts).
#'
#' Can be run standalone or sourced and called with validate_dataset_files()

library(yaml)

#' Extract dataset ID from filename
#'
#' @param filename Expression file name
#' @return Dataset ID (e.g., "GSE100051") or NA
extract_dataset_id <- function(filename) {

  base <- sub("\\.tsv$", "", basename(filename))
  if (grepl("^GSE[0-9]+", base)) {
    return(regmatches(base, regexpr("^GSE[0-9]+", base)))
  } else if (grepl("^E-GEOD-[0-9]+", base)) {
    egeod <- regmatches(base, regexpr("^E-GEOD-[0-9]+", base))
    return(sub("E-GEOD-", "GSE", egeod))
  }
  return(NA)
}

#' Validate dataset files and metadata matching
#'
#' @param mapped_path Path to mapped expression files
#' @param phenodata_path Path to phenodata CSV file
#' @param datasets Optional vector of dataset IDs to check (NULL = all)
#' @param pattern File pattern for expression files
#' @param sample_col Column name for sample IDs in phenodata
#' @param batch_col Column name for dataset/batch IDs in phenodata
#' @param check_metadata_fields Check for required metadata fields
#' @param verbose Print detailed output
#' @return List with validation results
validate_dataset_files <- function(mapped_path = "data/mapped",
                                    phenodata_path = "data/phenodata/samples.csv",
                                    datasets = NULL,
                                    pattern = "\\.tsv$",
                                    sample_col = "arraydatafile_exprscolumnnames",
                                    batch_col = "secondaryaccession",
                                    check_metadata_fields = TRUE,
                                    verbose = TRUE) {

  results <- list(
    valid = TRUE,
    duplicate_files = list(),
    sample_mismatches = data.frame(),
    metadata_coverage = data.frame(),
    warnings = character()
  )

  # Get all expression files
  all_files <- list.files(mapped_path, pattern = pattern, full.names = FALSE)

  if (verbose) {
    cat("\n=== Checking for Duplicate Dataset Files ===\n\n")
    cat("Found", length(all_files), "expression files in", mapped_path, "\n\n")
  }

  # Map files to dataset IDs
  file_dataset_map <- data.frame(
    file = all_files,
    dataset_id = sapply(all_files, extract_dataset_id),
    stringsAsFactors = FALSE
  )

  # Filter to requested datasets if specified
  if (!is.null(datasets)) {
    check_datasets <- datasets
  } else {
    check_datasets <- unique(file_dataset_map$dataset_id)
    check_datasets <- check_datasets[!is.na(check_datasets)]
  }

  # Check 1: Duplicate files for same dataset
  if (verbose) cat("Checking for duplicate dataset files...\n")

  for (ds in check_datasets) {
    matching_files <- file_dataset_map$file[file_dataset_map$dataset_id == ds]

    if (length(matching_files) > 1) {
      results$valid <- FALSE
      results$duplicate_files[[ds]] <- matching_files
      warning_msg <- paste0(ds, ": Multiple files - ",
                            paste(matching_files, collapse = ", "))
      results$warnings <- c(results$warnings, warning_msg)

      if (verbose) {
        cat("  ⚠️  ", ds, ": Multiple files found!\n", sep = "")
        for (f in matching_files) cat("      - ", f, "\n", sep = "")

        # Check column overlap
        cols_list <- lapply(matching_files, function(f) {
          tryCatch({
            h <- read.table(file.path(mapped_path, f), header = TRUE,
                            sep = "\t", nrows = 1, check.names = FALSE,
                            row.names = 1)
            colnames(h)
          }, error = function(e) character(0))
        })

        if (length(cols_list) >= 2) {
          overlap <- Reduce(intersect, cols_list)
          if (length(overlap) > 0) {
            cat("      → ", length(overlap),
                " overlapping samples will get .x/.y suffixes!\n", sep = "")
          }
        }
      }
    } else if (length(matching_files) == 1) {
      if (verbose) cat("  ✓ ", ds, ": ", matching_files, "\n", sep = "")
    } else {
      results$valid <- FALSE
      warning_msg <- paste0(ds, ": No expression file found")
      results$warnings <- c(results$warnings, warning_msg)
      if (verbose) cat("  ⚠️  ", ds, ": No expression file found!\n", sep = "")
    }
  }

  # Check 2: Metadata-expression sample matching
  if (verbose) cat("\n=== Checking Metadata-Expression Sample Matching ===\n\n")

  if (file.exists(phenodata_path)) {
    pdata <- read.csv(phenodata_path, stringsAsFactors = FALSE)

    mismatch_summary <- data.frame(
      dataset = character(),
      metadata_samples = integer(),
      expr_columns = integer(),
      matched = integer(),
      in_metadata_only = integer(),
      in_expr_only = integer(),
      missing_samples = character(),
      stringsAsFactors = FALSE
    )

    for (ds in check_datasets) {
      # Get metadata samples (exclude samples marked as excluded)
      ds_mask <- pdata[[batch_col]] == ds
      if ("X.excluded." %in% colnames(pdata)) {
        ds_mask <- ds_mask & (is.na(pdata$X.excluded.) | pdata$X.excluded. == FALSE)
      }
      meta_samples <- pdata[[sample_col]][ds_mask]
      meta_samples <- meta_samples[!is.na(meta_samples) & meta_samples != ""]

      # Find expression files
      dataset_files <- file_dataset_map$file[file_dataset_map$dataset_id == ds]
      if (length(dataset_files) == 0) next

      # Get expression columns
      expr_columns <- character()
      for (f in dataset_files) {
        tryCatch({
          h <- read.table(file.path(mapped_path, f), header = TRUE,
                          sep = "\t", nrows = 1, check.names = FALSE,
                          row.names = 1)
          expr_columns <- union(expr_columns, colnames(h))
        }, error = function(e) NULL)
      }

      # Compare
      matched <- sum(meta_samples %in% expr_columns)
      missing_from_expr <- meta_samples[!meta_samples %in% expr_columns]
      in_meta_only <- length(missing_from_expr)
      in_expr_only <- sum(!expr_columns %in% meta_samples)

      mismatch_summary <- rbind(mismatch_summary, data.frame(
        dataset = ds,
        metadata_samples = length(meta_samples),
        expr_columns = length(expr_columns),
        matched = matched,
        in_metadata_only = in_meta_only,
        in_expr_only = in_expr_only,
        missing_samples = paste(missing_from_expr, collapse = ";"),
        stringsAsFactors = FALSE
      ))

      if (in_meta_only > 0) {
        results$valid <- FALSE
        warning_msg <- paste0(ds, ": ", in_meta_only,
                              " samples in metadata but not in expression")
        results$warnings <- c(results$warnings, warning_msg)

        if (verbose) {
          cat("  ⚠️  ", ds, ": ", in_meta_only,
              " samples in metadata but NOT in expression\n", sep = "")
          if (in_meta_only <= 5) {
            cat("      Missing: ", paste(missing_from_expr, collapse = ", "),
                "\n", sep = "")
          }
        }
      } else if (matched > 0) {
        if (verbose) cat("  ✓ ", ds, ": ", matched, " samples matched\n", sep = "")
      }
    }

    results$sample_mismatches <- mismatch_summary

    # Check 3: Metadata field coverage (Diagnosis, Gestational Category, etc.)
    if (check_metadata_fields) {
      if (verbose) {
        cat("\n=== Checking Metadata Field Coverage ===\n\n")
      }

      coverage_summary <- data.frame(
        dataset = character(),
        n_samples = integer(),
        has_diagnosis = logical(),
        has_healthy = logical(),
        has_gestational_category = logical(),
        has_biological_specimen = logical(),
        stringsAsFactors = FALSE
      )

      for (ds in check_datasets) {
        ds_mask <- pdata[[batch_col]] == ds
        if ("X.excluded." %in% colnames(pdata)) {
          ds_mask <- ds_mask & (is.na(pdata$X.excluded.) | pdata$X.excluded. == FALSE)
        }
        ds_samples <- pdata[ds_mask, ]
        if (nrow(ds_samples) == 0) next

        has_diag <- any(
          !is.na(ds_samples$Diagnosis) &
          ds_samples$Diagnosis != "" &
          ds_samples$Diagnosis != "NA",
          na.rm = TRUE
        )

        has_healthy <- any(
          ds_samples$Diagnosis == "Healthy",
          na.rm = TRUE
        )

        has_gest <- any(
          !is.na(ds_samples$Gestational.Age.Category) &
          ds_samples$Gestational.Age.Category != "" &
          ds_samples$Gestational.Age.Category != "NA" &
          ds_samples$Gestational.Age.Category != "Unknown Gestational Category",
          na.rm = TRUE
        )

        has_bio <- any(
          !is.na(ds_samples$Biological.Specimen) &
          ds_samples$Biological.Specimen != "" &
          ds_samples$Biological.Specimen != "NA",
          na.rm = TRUE
        )

        coverage_summary <- rbind(coverage_summary, data.frame(
          dataset = ds,
          n_samples = nrow(ds_samples),
          has_diagnosis = has_diag,
          has_healthy = has_healthy,
          has_gestational_category = has_gest,
          has_biological_specimen = has_bio,
          stringsAsFactors = FALSE
        ))

        missing_fields <- c()
        if (!has_diag) missing_fields <- c(missing_fields, "Diagnosis")
        if (!has_healthy) missing_fields <- c(missing_fields, "Healthy samples")
        if (!has_gest) missing_fields <- c(missing_fields, "Gestational.Age.Category")
        if (!has_bio) missing_fields <- c(missing_fields, "Biological.Specimen")

        if (length(missing_fields) > 0) {
          results$valid <- FALSE
          warning_msg <- paste0(ds, ": Missing ",
                                paste(missing_fields, collapse = ", "))
          results$warnings <- c(results$warnings, warning_msg)

          if (verbose) {
            cat("  ⚠️  ", ds, ": Missing ",
                paste(missing_fields, collapse = ", "), "\n", sep = "")
          }
        } else {
          if (verbose) {
            cat("  ✓ ", ds, ": All required fields present\n", sep = "")
          }
        }
      }

      results$metadata_coverage <- coverage_summary
    }
  }

  if (verbose) {
    cat("\n")
    if (!results$valid) {
      cat("───────────────────────────────────────────────────────────\n")
      cat("⚠️  VALIDATION WARNINGS FOUND\n")
      cat("   Some samples may be lost during analysis.\n")
      cat("───────────────────────────────────────────────────────────\n")
    } else {
      cat("✓ All datasets validated successfully\n")
    }
    cat("\n")
  }

  return(results)
}

#' Load config and run validation
#'
#' @param config_file Path to YAML config file
#' @param verbose Print detailed output
#' @return Validation results
validate_from_config <- function(config_file, verbose = TRUE) {
  config <- yaml::read_yaml(config_file)

  # Extract paths from config

  mapped_path <- config$paths$mapped_data
  phenodata_path <- config$paths$phenodata
  pattern <- config$merging$pattern
  sample_col <- config$merging$sample_column
  batch_col <- config$batch_correction$batch_column

  # Get all datasets to check
  datasets <- c(
    config$baseline$datasets,
    if (!is.null(config$term_datasets)) config$term_datasets$datasets,
    if (!is.null(config$first_trimester_datasets))
      config$first_trimester_datasets$datasets
  )
  datasets <- unique(datasets)

  if (verbose) {
    cat("\n=== Validating datasets from config: ", config_file, " ===\n", sep = "")
    cat("Datasets to check: ", paste(datasets, collapse = ", "), "\n", sep = "")
  }

  validate_dataset_files(
    mapped_path = mapped_path,
    phenodata_path = phenodata_path,
    datasets = datasets,
    pattern = pattern,
    sample_col = sample_col,
    batch_col = batch_col,
    verbose = verbose
  )
}


# Run standalone if executed directly (not sourced)
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) > 0 && file.exists(args[1])) {
    # Config file provided
    results <- validate_from_config(args[1])
  } else {
    # Use defaults
    results <- validate_dataset_files()
  }

  # Save results
  if (nrow(results$sample_mismatches) > 0) {
    write.csv(results$sample_mismatches, "dataset_sample_matching.csv",
              row.names = FALSE)
    cat("✓ Results saved to: dataset_sample_matching.csv\n")
  }

  cat("\n✓ Check complete!\n\n")

  # Exit with error code if validation failed
  if (!results$valid) {
    quit(status = 1)
  }
}
