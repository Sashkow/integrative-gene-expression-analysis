#!/usr/bin/env Rscript

# Install missing packages for integrative analysis pipeline
#
# SYSTEM DEPENDENCIES:
# Some packages require system libraries. On Ubuntu/Debian, install:
#   sudo apt-get install libgmp-dev libmpfr-dev
#
# These are needed for: RankProd (via gmp, Rmpfr dependencies)

# ============================================================================
# CRAN packages
# ============================================================================
cran_packages <- c(
  "VennDiagram",
  "readxl",
  "yaml",
  "metafor",            # General meta-analysis (Phase 2 alternative)
  "softImpute",         # Matrix completion imputation (Phase 2B)
  "cluster",            # Silhouette scores for normalization comparison
  "ruv"                 # RUVinv batch correction (auto k selection)
)

cat("=== Checking CRAN packages ===\n")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cran.r-project.org")
  } else {
    cat(pkg, "is already installed\n")
  }
}

# ============================================================================
# Bioconductor packages
# ============================================================================
bioc_packages <- c(
  # Core pipeline
  "limma",
  "sva",              # ComBat batch correction
  "org.Hs.eg.db",     # Gene annotation


  # Phase 2: Meta-analysis
  "DExMA",            # Effect-size meta-analysis with imputation
  "RankProd",         # Rank-based meta-analysis

  # Batch correction
  "HarmonizR"         # Matrix dissection + ComBat/limma (Voß et al. 2022)
)

cat("\n=== Checking Bioconductor packages ===\n")

# Ensure BiocManager is available
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager", repos = "https://cran.r-project.org")
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "from Bioconductor...\n")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat(pkg, "is already installed\n")
  }
}

cat("\nAll packages installed successfully!\n")
