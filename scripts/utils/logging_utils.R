#' Logging and Archiving Utilities
#'
#' Common functions for logging output and archiving previous results.
#' Reusable across all analysis scripts.
#'
#' @author Expression Integration Pipeline
#' @date 2025

# ============================================================================
# Archiving Functions
# ============================================================================

#' Archive previous results in output directory
#'
#' Moves existing files (except archive folder and logs) to a timestamped
#' subdirectory within archive/.
#'
#' @param output_dir Path to output directory
#' @param exclude_patterns Regex patterns to exclude from archiving
#' @return Path to archive subdirectory (or NULL if nothing archived)
#' @export
archive_previous_results <- function(output_dir,
                                     exclude_patterns = "^archive$|^run_log_") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    return(NULL)
  }

  existing_files <- list.files(output_dir, full.names = FALSE)
  existing_files <- existing_files[!grepl(exclude_patterns, existing_files)]

  if (length(existing_files) == 0) {
    return(NULL)
  }

  archive_dir <- file.path(output_dir, "archive")
  dir.create(archive_dir, showWarnings = FALSE, recursive = TRUE)

  timestamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  archive_subdir <- file.path(archive_dir, timestamp)
  dir.create(archive_subdir, showWarnings = FALSE)

  for (f in existing_files) {
    file.rename(file.path(output_dir, f), file.path(archive_subdir, f))
  }

  cat("Archived", length(existing_files), "files to:", archive_subdir, "\n")
  return(archive_subdir)
}


# ============================================================================
# Logging Functions
# ============================================================================

#' Setup logging to file
#'
#' Redirects console output to both file and console using sink().
#' Call close_logging() when done.
#'
#' @param output_dir Directory to save log file
#' @param prefix Log filename prefix (default: "run_log")
#' @return Path to log file
#' @export
setup_logging <- function(output_dir, prefix = "run_log") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  timestamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  log_file <- file.path(output_dir, paste0(prefix, "_", timestamp, ".txt"))

  # Sink output to file (split = TRUE keeps console output)
  sink(log_file, split = TRUE)

  cat("Log started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Log file:", log_file, "\n\n")

  return(log_file)
}


#' Close logging
#'
#' Closes the sink connection opened by setup_logging().
#'
#' @param log_file Optional path to log file for final message
#' @export
close_logging <- function(log_file = NULL) {
  cat("\nLog ended:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

  sink()

  if (!is.null(log_file)) {
    message("Log saved to: ", log_file)
  }
}


#' Log a section header
#'
#' Prints a formatted section header for readability.
#'
#' @param title Section title
#' @param width Total width of header line (default: 60)
#' @export
log_section <- function(title, width = 60) {
  cat("\n")
  cat(paste(rep("=", width), collapse = ""), "\n")
  cat(" ", title, "\n")
  cat(paste(rep("=", width), collapse = ""), "\n\n")
}


#' Log a subsection header
#'
#' @param title Subsection title
#' @export
log_subsection <- function(title) {
  cat("\n--- ", title, " ---\n\n")
}


#' Log key-value pair
#'
#' @param key Parameter name
#' @param value Parameter value
#' @param indent Number of spaces to indent (default: 2)
#' @export
log_param <- function(key, value, indent = 2) {
  spaces <- paste(rep(" ", indent), collapse = "")
  cat(sprintf("%s%-20s %s\n", spaces, paste0(key, ":"), value))
}


#' Log elapsed time
#'
#' @param start_time Start time from Sys.time()
#' @param task_name Name of task for message
#' @export
log_elapsed <- function(start_time, task_name = "Task") {
  elapsed <- difftime(Sys.time(), start_time, units = "secs")
  if (elapsed < 60) {
    cat(sprintf("%s completed in %.1f seconds\n", task_name, elapsed))
  } else {
    cat(sprintf("%s completed in %.1f minutes\n", task_name, elapsed / 60))
  }
}
