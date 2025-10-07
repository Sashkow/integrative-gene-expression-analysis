#' Utility Functions
#'
#' Helper functions for the expression integration pipeline
#'
#' @author Expression Integration Pipeline
#' @date 2025

#' Save data to file
#'
#' @param data Data to save
#' @param filepath File path
#' @param format File format ("csv", "tsv", "rds")
#' @export
save_data <- function(data, filepath, format = "csv") {

  dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)

  if (format == "csv") {
    write.csv(data, filepath, row.names = FALSE)
  } else if (format == "tsv") {
    write.table(data, filepath, sep = "\t", row.names = FALSE, quote = FALSE)
  } else if (format == "rds") {
    saveRDS(data, filepath)
  } else {
    stop("Unknown format: ", format)
  }

  message("Saved: ", filepath)
  invisible(filepath)
}


#' Load data from file
#'
#' @param filepath File path
#' @param format File format ("csv", "tsv", "rds")
#' @return Loaded data
#' @export
load_data <- function(filepath, format = "csv") {

  if (!file.exists(filepath)) {
    stop("File not found: ", filepath)
  }

  if (format == "csv") {
    data <- read.csv(filepath, stringsAsFactors = FALSE)
  } else if (format == "tsv") {
    data <- read.table(filepath, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else if (format == "rds") {
    data <- readRDS(filepath)
  } else {
    stop("Unknown format: ", format)
  }

  message("Loaded: ", filepath, " (", nrow(data), " rows)")
  return(data)
}


#' Create timestamp for output files
#'
#' @param format Date format (default: "%Y%m%d_%H%M%S")
#' @return Timestamp string
#' @export
timestamp <- function(format = "%Y%m%d_%H%M%S") {
  format(Sys.time(), format)
}


#' Print session information
#'
#' @export
print_session_info <- function() {
  message("\n=== Session Information ===")
  print(sessionInfo())
  message("===========================\n")
}


#' Check required packages
#'
#' @param packages Vector of package names
#' @export
check_packages <- function(packages) {

  missing <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]

  if (length(missing) > 0) {
    stop("Missing required packages: ", paste(missing, collapse = ", "),
         "\nInstall with: BiocManager::install(c('", paste(missing, collapse = "', '"), "'))")
  }

  message("All required packages available")
  invisible(TRUE)
}


#' Create output directory structure
#'
#' @param base_dir Base output directory
#' @return List of created directories
#' @export
create_output_dirs <- function(base_dir) {

  dirs <- list(
    base = base_dir,
    difexp = file.path(base_dir, "difexp"),
    clusters = file.path(base_dir, "clusters"),
    enrichment = file.path(base_dir, "enrichment"),
    plots = file.path(base_dir, "plots"),
    reports = file.path(base_dir, "reports")
  )

  for (dir in dirs) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }

  message("Created output directory structure at: ", base_dir)
  return(dirs)
}


#' Write summary statistics
#'
#' @param stats Named list of statistics
#' @param output_file Output file path
#' @export
write_summary <- function(stats, output_file) {

  summary_df <- data.frame(
    Metric = names(stats),
    Value = unlist(stats),
    stringsAsFactors = FALSE
  )

  write.csv(summary_df, output_file, row.names = FALSE)
  message("Saved summary: ", output_file)

  invisible(summary_df)
}


#' Log message with timestamp
#'
#' @param ... Messages to log
#' @export
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
  message(msg)
}


#' Validate input data structure
#'
#' @param exprs Expression matrix
#' @param pdata Phenodata
#' @export
validate_input <- function(exprs, pdata) {

  errors <- c()

  # Check dimensions
  if (ncol(exprs) == 0) {
    errors <- c(errors, "Expression matrix has no samples")
  }
  if (nrow(exprs) == 0) {
    errors <- c(errors, "Expression matrix has no genes")
  }
  if (nrow(pdata) == 0) {
    errors <- c(errors, "Phenodata has no samples")
  }

  # Check for NAs
  if (any(is.na(colnames(exprs)))) {
    errors <- c(errors, "Expression matrix has NA column names")
  }
  if (any(is.na(rownames(exprs)))) {
    errors <- c(errors, "Expression matrix has NA row names")
  }

  if (length(errors) > 0) {
    stop("Input validation failed:\n  ", paste(errors, collapse = "\n  "))
  }

  message("Input validation passed")
  invisible(TRUE)
}
