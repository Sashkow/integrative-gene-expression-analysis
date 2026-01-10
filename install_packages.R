#!/usr/bin/env Rscript

# Install missing packages for comparison script
packages <- c("VennDiagram", "readxl")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cran.r-project.org")
  } else {
    cat(pkg, "is already installed\n")
  }
}

cat("\nAll packages installed successfully!\n")
