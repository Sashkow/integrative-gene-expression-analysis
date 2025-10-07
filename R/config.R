#' Configuration Management
#'
#' Functions for loading and managing pipeline configuration
#'
#' @author Expression Integration Pipeline
#' @date 2025

#' Load configuration from YAML file
#'
#' @param config_file Path to config YAML file
#' @return List with configuration parameters
#' @export
load_config <- function(config_file = "config/config.yaml") {

  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' required. Install with: install.packages('yaml')")
  }

  if (!file.exists(config_file)) {
    stop("Config file not found: ", config_file)
  }

  config <- yaml::read_yaml(config_file)

  message("Loaded configuration from: ", config_file)
  message("  Project: ", config$project$name)

  return(config)
}


#' Get configuration parameter
#'
#' @param config Configuration list
#' @param ... Path to parameter (e.g., "paths", "output")
#' @param default Default value if not found
#' @return Parameter value
#' @export
get_config <- function(config, ..., default = NULL) {

  keys <- list(...)
  value <- config

  for (key in keys) {
    if (key %in% names(value)) {
      value <- value[[key]]
    } else {
      if (!is.null(default)) {
        return(default)
      } else {
        stop("Configuration parameter not found: ", paste(keys, collapse = " -> "))
      }
    }
  }

  return(value)
}


#' Validate configuration
#'
#' @param config Configuration list
#' @export
validate_config <- function(config) {

  errors <- c()

  # Check required sections
  required_sections <- c("paths", "batch_correction", "differential_expression", "stringdb")

  for (section in required_sections) {
    if (!section %in% names(config)) {
      errors <- c(errors, paste("Missing required section:", section))
    }
  }

  # Check required paths
  if ("paths" %in% names(config)) {
    if (!"mapped_data" %in% names(config$paths)) {
      errors <- c(errors, "Missing required path: mapped_data")
    }
    if (!"phenodata" %in% names(config$paths)) {
      errors <- c(errors, "Missing required path: phenodata")
    }
  }

  # Check batch correction
  if ("batch_correction" %in% names(config)) {
    if (!"batch_column" %in% names(config$batch_correction)) {
      errors <- c(errors, "Missing batch_column in batch_correction")
    }
  }

  if (length(errors) > 0) {
    stop("Configuration validation failed:\n  ", paste(errors, collapse = "\n  "))
  }

  message("Configuration validation passed")
  invisible(TRUE)
}


#' Print configuration summary
#'
#' @param config Configuration list
#' @export
print_config <- function(config) {

  cat("\n=== Pipeline Configuration ===\n\n")

  cat("Project Information:\n")
  cat("  Name:", config$project$name, "\n")
  cat("  Author:", config$project$author, "\n")
  cat("  Date:", config$project$date, "\n\n")

  cat("Input Paths:\n")
  cat("  Mapped data:", config$paths$mapped_data, "\n")
  cat("  Phenodata:", config$paths$phenodata, "\n")
  cat("  Output:", config$paths$output, "\n\n")

  cat("Batch Correction:\n")
  cat("  Batch column:", config$batch_correction$batch_column, "\n")
  cat("  Model formula:", config$batch_correction$model_formula, "\n\n")

  cat("Differential Expression:\n")
  cat("  Design formula:", config$differential_expression$design_formula, "\n")
  cat("  Contrasts:", paste(config$differential_expression$contrasts, collapse = "; "), "\n")
  cat("  LogFC threshold:", config$differential_expression$logfc_threshold, "\n\n")

  cat("STRING Database:\n")
  cat("  Version:", config$stringdb$version, "\n")
  cat("  Species:", config$stringdb$species, "\n")
  cat("  Score threshold:", config$stringdb$score_threshold, "\n\n")

  cat("==============================\n\n")

  invisible(NULL)
}


#' Update configuration parameter
#'
#' @param config Configuration list
#' @param value New value
#' @param ... Path to parameter
#' @return Updated configuration list
#' @export
set_config <- function(config, value, ...) {

  keys <- list(...)
  temp <- config

  # Navigate to parent
  for (i in seq_len(length(keys) - 1)) {
    key <- keys[[i]]
    if (!key %in% names(temp)) {
      temp[[key]] <- list()
    }
    temp <- temp[[key]]
  }

  # Set value
  last_key <- keys[[length(keys)]]
  temp[[last_key]] <- value

  message("Updated config: ", paste(keys, collapse = " -> "), " = ", value)

  return(config)
}
