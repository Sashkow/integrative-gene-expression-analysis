#' Data Merging Module
#'
#' Functions for merging expression data from multiple studies/datasets
#'
#' @author Expression Integration Pipeline
#' @date 2025

#' Merge multiple expression matrices by common genes
#'
#' @param mapped_path Path to directory containing mapped expression files
#' @param pattern File pattern to match (default: "*.tsv")
#' @return List with merged expression matrix and metadata
#' @export
merge_expression_data <- function(mapped_path, pattern = "\\.tsv$") {

  # Get list of expression files
  exprs_files <- list.files(mapped_path, pattern = pattern, full.names = FALSE)

  if (length(exprs_files) == 0) {
    stop("No expression files found in ", mapped_path)
  }

  message("Found ", length(exprs_files), " expression files to merge")

  # Read first file
  mrgd <- read.table(file.path(mapped_path, exprs_files[1]),
                     header = TRUE, sep = '\t', row.names = 1)

  message("Starting with ", nrow(mrgd), " genes from ", exprs_files[1])

  # Merge remaining files
  for (exprs_file in exprs_files[-1]) {
    current_exprs <- read.table(file.path(mapped_path, exprs_file),
                                header = TRUE, sep = '\t', row.names = 1)

    mrgd <- merge(mrgd, current_exprs, by = "row.names", all = FALSE)
    rownames(mrgd) <- mrgd$Row.names
    mrgd <- mrgd[, !(colnames(mrgd) == "Row.names")]

    message("  Merged ", exprs_file, " - ", nrow(mrgd), " common genes")
  }

  message("Final merged data: ", nrow(mrgd), " genes x ", ncol(mrgd), " samples")

  return(mrgd)
}


#' Align phenodata with expression data
#'
#' @param mrgd Merged expression matrix
#' @param pdata Phenotype/metadata dataframe
#' @param sample_col Column name in pdata containing sample identifiers
#' @return Aligned phenodata dataframe
#' @export
align_phenodata <- function(mrgd, pdata, sample_col = "arraydatafile_exprscolumnnames") {

  # Get samples in expression data
  pdata_aligned <- pdata[make.names(pdata[[sample_col]]) %in% colnames(mrgd), ]

  message("Aligned ", nrow(pdata_aligned), " samples between expression and phenodata")

  # Check alignment
  missing_exprs <- setdiff(colnames(mrgd), make.names(pdata_aligned[[sample_col]]))
  if (length(missing_exprs) > 0) {
    warning("Missing metadata for ", length(missing_exprs), " expression samples")
  }

  return(pdata_aligned)
}


#' Filter samples by criteria
#'
#' @param mrgd Expression matrix
#' @param pdata Phenodata
#' @param filters Named list of filter criteria
#' @return List with filtered expression matrix and phenodata
#' @export
filter_samples <- function(mrgd, pdata, filters = list()) {

  # Start with all samples
  keep_idx <- rep(TRUE, nrow(pdata))

  # Apply filters
  for (col in names(filters)) {
    if (col %in% colnames(pdata)) {
      values <- filters[[col]]
      keep_idx <- keep_idx & (pdata[[col]] %in% values)
      message("Filter '", col, "' in [", paste(values, collapse = ", "),
              "]: ", sum(keep_idx), " samples")
    } else {
      warning("Column '", col, "' not found in phenodata")
    }
  }

  # Filter data
  pdata_filtered <- pdata[keep_idx, ]
  mrgd_filtered <- mrgd[, make.names(pdata_filtered$arraydatafile_exprscolumnnames)]

  message("Final filtered data: ", ncol(mrgd_filtered), " samples")

  return(list(
    exprs = mrgd_filtered,
    pdata = pdata_filtered
  ))
}


#' Filter phenodata by criteria (metadata-only filtering)
#'
#' @param pdata Phenodata dataframe
#' @param filters Named list of filter criteria
#' @return Filtered phenodata dataframe
#' @export
filter_phenodata <- function(pdata, filters = list()) {

  # Start with all samples
  keep_idx <- rep(TRUE, nrow(pdata))

  # Apply filters
  for (col in names(filters)) {
    if (col %in% colnames(pdata)) {
      values <- filters[[col]]
      keep_idx <- keep_idx & (pdata[[col]] %in% values)
      message("Filter '", col, "' in [", paste(values, collapse = ", "),
              "]: ", sum(keep_idx), " samples")
    } else {
      warning("Column '", col, "' not found in phenodata")
    }
  }

  pdata_filtered <- pdata[keep_idx, ]

  message("Filtered phenodata: ", nrow(pdata_filtered), " samples")

  return(pdata_filtered)
}


#' Identify datasets needed based on filtered phenodata
#'
#' @param pdata_filtered Filtered phenodata
#' @param batch_col Column name containing dataset/batch identifier
#' @return Vector of unique dataset identifiers
#' @export
identify_datasets <- function(pdata_filtered, batch_col = "secondaryaccession") {

  datasets <- unique(pdata_filtered[[batch_col]])
  datasets <- datasets[!is.na(datasets) & datasets != ""]

  message("Identified ", length(datasets), " datasets needed: ",
          paste(datasets, collapse = ", "))

  return(datasets)
}


#' Print metadata summary table for selected samples
#'
#' @param pdata_filtered Filtered phenodata
#' @param batch_col Column name containing dataset/batch identifier
#' @param group_col Column name for grouping (e.g., trimester)
#' @export
print_metadata_summary <- function(pdata_filtered,
                                   batch_col = "secondaryaccession",
                                   group_col = "Gestational.Age.Category") {

  cat("\n=== Selected Samples Metadata Summary ===\n\n")

  # Overall summary
  cat("Total samples selected:", nrow(pdata_filtered), "\n")
  cat("Number of datasets:", length(unique(pdata_filtered[[batch_col]])), "\n\n")

  # Create summary table by dataset and group
  if (group_col %in% colnames(pdata_filtered)) {
    # Cross-tabulation
    summary_table <- table(pdata_filtered[[batch_col]],
                           pdata_filtered[[group_col]])

    # Convert to data frame for better formatting
    summary_df <- as.data.frame.matrix(summary_table)
    summary_df$Dataset <- rownames(summary_df)
    summary_df$Total <- rowSums(summary_table)

    # Reorder columns
    summary_df <- summary_df[, c("Dataset", colnames(summary_df)[colnames(summary_df) != "Dataset"])]

    cat("Samples by Dataset and", group_col, ":\n")
    cat(paste(rep("-", 80), collapse = ""), "\n")

    # Print header
    header <- sprintf("%-15s", "Dataset")
    for (col in colnames(summary_df)[-1]) {
      if (col == "Total") {
        header <- paste0(header, sprintf("%10s", col))
      } else {
        # Truncate long column names
        col_short <- substr(col, 1, 10)
        header <- paste0(header, sprintf("%12s", col_short))
      }
    }
    cat(header, "\n")
    cat(paste(rep("-", 80), collapse = ""), "\n")

    # Print rows
    for (i in 1:nrow(summary_df)) {
      row_str <- sprintf("%-15s", substr(summary_df$Dataset[i], 1, 15))
      for (j in 2:ncol(summary_df)) {
        if (colnames(summary_df)[j] == "Total") {
          row_str <- paste0(row_str, sprintf("%10d", summary_df[i, j]))
        } else {
          row_str <- paste0(row_str, sprintf("%12d", summary_df[i, j]))
        }
      }
      cat(row_str, "\n")
    }
    cat(paste(rep("-", 80), collapse = ""), "\n")

    # Print totals
    total_str <- sprintf("%-15s", "TOTAL")
    for (j in 2:ncol(summary_df)) {
      if (colnames(summary_df)[j] == "Total") {
        total_str <- paste0(total_str, sprintf("%10d", sum(summary_df[, j])))
      } else {
        total_str <- paste0(total_str, sprintf("%12d", sum(summary_df[, j])))
      }
    }
    cat(total_str, "\n")
    cat(paste(rep("-", 80), collapse = ""), "\n\n")

  } else {
    # Simple count by dataset
    dataset_counts <- table(pdata_filtered[[batch_col]])
    cat("Samples by Dataset:\n")
    print(dataset_counts)
    cat("\n")
  }

  # Additional metadata if available
  if ("Diagnosis" %in% colnames(pdata_filtered)) {
    cat("Diagnosis distribution:\n")
    print(table(pdata_filtered$Diagnosis))
    cat("\n")
  }

  if ("Biological.Specimen" %in% colnames(pdata_filtered)) {
    cat("Biological Specimen distribution:\n")
    print(table(pdata_filtered$Biological.Specimen))
    cat("\n")
  }

  if ("estimated_sex" %in% colnames(pdata_filtered)) {
    sex_counts <- table(pdata_filtered$estimated_sex, useNA = "ifany")
    cat("Estimated Sex distribution:\n")
    print(sex_counts)
    cat("\n")
  }

  cat("=========================================\n\n")
}


#' Merge expression data from specific datasets only
#'
#' @param mapped_path Path to directory containing mapped expression files
#' @param datasets Vector of dataset identifiers to include
#' @param pattern File pattern to match
#' @return Merged expression matrix
#' @export
merge_selected_datasets <- function(mapped_path, datasets, pattern = "\\.tsv$") {

  # Get all expression files
  all_files <- list.files(mapped_path, pattern = pattern, full.names = FALSE)

  # Filter files to only include selected datasets
  selected_files <- character(0)
  for (dataset in datasets) {
    # Match files containing the dataset ID
    # Try multiple patterns: GSE12345, E-GEOD-12345, E-GEOD-GSE12345
    patterns_to_try <- c(
      dataset,                           # GSE12345
      gsub("^GSE", "", dataset),         # 12345 (if dataset starts with GSE)
      paste0("E-GEOD-", gsub("^GSE", "", dataset)),  # E-GEOD-12345
      paste0("E-GEOD-", dataset)         # E-GEOD-GSE12345
    )

    matching <- character(0)
    for (pat in patterns_to_try) {
      found <- grep(pat, all_files, value = TRUE, fixed = TRUE)
      matching <- c(matching, found)
    }

    matching <- unique(matching)

    if (length(matching) > 0) {
      selected_files <- c(selected_files, matching)
      message("  Found file(s) for ", dataset, ": ", paste(matching, collapse = ", "))
    } else {
      warning("No file found for dataset: ", dataset)
    }
  }

  selected_files <- unique(selected_files)

  if (length(selected_files) == 0) {
    stop("No expression files found for datasets: ", paste(datasets, collapse = ", "))
  }

  message("Merging ", length(selected_files), " selected dataset files:")
  for (f in selected_files) {
    message("  - ", f)
  }

  # Read first file
  mrgd <- read.table(file.path(mapped_path, selected_files[1]),
                     header = TRUE, sep = '\t', row.names = 1)

  message("Starting with ", nrow(mrgd), " genes from ", selected_files[1])

  # Merge remaining files
  if (length(selected_files) > 1) {
    for (exprs_file in selected_files[-1]) {
      current_exprs <- read.table(file.path(mapped_path, exprs_file),
                                  header = TRUE, sep = '\t', row.names = 1)

      mrgd <- merge(mrgd, current_exprs, by = "row.names", all = FALSE)
      rownames(mrgd) <- mrgd$Row.names
      mrgd <- mrgd[, !(colnames(mrgd) == "Row.names")]

      message("  Merged ", exprs_file, " - ", nrow(mrgd), " common genes")
    }
  }

  message("Final merged data: ", nrow(mrgd), " genes x ", ncol(mrgd), " samples")

  return(mrgd)
}


#' Filter expression matrix to only include samples from filtered phenodata
#'
#' @param mrgd Merged expression matrix
#' @param pdata_filtered Filtered phenodata
#' @param sample_col Column name in pdata containing sample identifiers
#' @return Filtered expression matrix
#' @export
subset_expression_by_phenodata <- function(mrgd, pdata_filtered, sample_col = "arraydatafile_exprscolumnnames") {

  # Get sample names from phenodata
  target_samples <- make.names(pdata_filtered[[sample_col]])

  # Get samples present in expression data
  available_samples <- intersect(target_samples, colnames(mrgd))

  if (length(available_samples) == 0) {
    stop("No samples from filtered phenodata found in expression matrix")
  }

  # Subset expression matrix
  mrgd_subset <- mrgd[, available_samples]

  message("Subset expression matrix: ", ncol(mrgd_subset), " samples (from ",
          nrow(pdata_filtered), " requested)")

  if (length(available_samples) < length(target_samples)) {
    missing <- setdiff(target_samples, available_samples)
    warning("Missing ", length(missing), " samples from expression data")
  }

  return(mrgd_subset)
}
