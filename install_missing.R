# Install missing packages
user_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE)
}
.libPaths(c(user_lib, .libPaths()))

cat("Installing to:", .libPaths()[1], "\n\n")

missing_cran <- c("factoextra", "pheatmap", "ggrepel", "writexl", "enrichR")
missing_bioc <- c("org.Hs.eg.db", "AnnotationDbi", "STRINGdb")

# Install CRAN packages
for (pkg in missing_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, lib = .libPaths()[1], repos = "https://cran.r-project.org")
  } else {
    cat(pkg, "already installed\n")
  }
}

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = .libPaths()[1], repos = "https://cran.r-project.org")
}

library(BiocManager)
for (pkg in missing_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    BiocManager::install(pkg, update = FALSE, ask = FALSE, lib = .libPaths()[1])
  } else {
    cat(pkg, "already installed\n")
  }
}

cat("\n=== Installation Complete ===\n")
cat("Installed packages location:", .libPaths()[1], "\n")
