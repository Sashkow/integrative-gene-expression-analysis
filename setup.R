#!/usr/bin/env Rscript

#' Setup Script for Expression Integration Pipeline
#'
#' This script installs required packages and sets up the pipeline environment
#'
#' @author Expression Integration Pipeline
#' @date 2025-10-04

cat("\n=== Expression Integration Pipeline Setup ===\n\n")

# Set up user library path
user_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(user_lib)) {
  cat("Creating user library directory:", user_lib, "\n")
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}

# Ensure user library is first in search path
.libPaths(c(user_lib, .libPaths()))
cat("Using library:", .libPaths()[1], "\n\n")

# Check R version
r_version <- getRversion()
cat("R version:", as.character(r_version), "\n")

if (r_version < "4.0.0") {
  stop("R version 4.0.0 or higher is required")
}

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager to user library...\n")
  install.packages("BiocManager", lib = .libPaths()[1], repos = "https://cran.r-project.org")
}

# Required Bioconductor packages
bioc_packages <- c(
  "limma",
  "sva",
  "org.Hs.eg.db",
  "STRINGdb",
  "AnnotationDbi"
)

# Required CRAN packages
cran_packages <- c(
  "yaml",
  "ggplot2",
  "factoextra",
  "igraph",
  "pheatmap",
  "ggrepel",
  "enrichR",
  "writexl",
  "rmarkdown",
  "knitr"
)

# Function to check and install packages
install_if_missing <- function(pkg, repo = "bioc") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    if (repo == "bioc") {
      BiocManager::install(pkg, update = FALSE, ask = FALSE, lib = .libPaths()[1])
    } else {
      install.packages(pkg, lib = .libPaths()[1], repos = "https://cran.r-project.org", quiet = TRUE)
    }
    return(TRUE)
  }
  return(FALSE)
}

# Install Bioconductor packages
cat("\nChecking Bioconductor packages...\n")
for (pkg in bioc_packages) {
  if (install_if_missing(pkg, "bioc")) {
    cat("  ✓", pkg, "installed\n")
  } else {
    cat("  ✓", pkg, "already installed\n")
  }
}

# Install CRAN packages
cat("\nChecking CRAN packages...\n")
for (pkg in cran_packages) {
  if (install_if_missing(pkg, "cran")) {
    cat("  ✓", pkg, "installed\n")
  } else {
    cat("  ✓", pkg, "already installed\n")
  }
}

# Create .gitkeep files for empty directories
dirs_to_keep <- c(
  "data/raw",
  "data/mapped",
  "data/merged",
  "data/temp",
  "output"
)

for (dir in dirs_to_keep) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  gitkeep_file <- file.path(dir, ".gitkeep")
  if (!file.exists(gitkeep_file)) {
    file.create(gitkeep_file)
  }
}

# Test loading modules
cat("\nTesting pipeline modules...\n")

test_module <- function(module_path) {
  tryCatch({
    source(module_path)
    cat("  ✓", basename(module_path), "\n")
    return(TRUE)
  }, error = function(e) {
    cat("  ✗", basename(module_path), "-", e$message, "\n")
    return(FALSE)
  })
}

modules <- list.files("R", pattern = "\\.R$", full.names = TRUE)
all_ok <- all(sapply(modules, test_module))

# Verify configuration
cat("\nVerifying configuration...\n")
if (file.exists("config/config.yaml")) {
  cat("  ✓ config.yaml found\n")
} else {
  cat("  ✗ config.yaml not found\n")
  all_ok <- FALSE
}

# Print summary
cat("\n=== Setup Summary ===\n")
cat("R version:", as.character(r_version), "\n")
cat("Packages installed:", length(c(bioc_packages, cran_packages)), "\n")
cat("Modules loaded:", length(modules), "\n")

if (all_ok) {
  cat("\n✓ Setup completed successfully!\n")
  cat("\nNext steps:\n")
  cat("1. Place your expression files in data/mapped/\n")
  cat("2. Place your phenodata in data/phenodata/samples.csv\n")
  cat("3. Edit config/config.yaml to customize pipeline parameters\n")
  cat("4. Run: rmarkdown::render('notebooks/pipeline.Rmd')\n\n")
} else {
  cat("\n✗ Setup completed with errors. Please check the messages above.\n\n")
}

# Print session info
cat("\n=== Session Information ===\n")
sessionInfo()
