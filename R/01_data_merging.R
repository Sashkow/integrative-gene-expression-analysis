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
