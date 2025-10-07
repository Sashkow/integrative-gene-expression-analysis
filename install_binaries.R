# Force binary installation instead of source compilation
options(install.packages.compile.from.source = "never")
options(repos = c(CRAN = "https://cloud.r-project.org"))

user_lib <- Sys.getenv("R_LIBS_USER")
.libPaths(c(user_lib, .libPaths()))

cat("Installing to:", .libPaths()[1], "\n\n")

# Try binary first, source if needed
install_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, lib = .libPaths()[1], type = "binary")
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("  Binary failed, trying source...\n")
      install.packages(pkg, lib = .libPaths()[1], type = "source")
    }
  } else {
    cat(pkg, "already installed\n")
  }
}

# Install critical missing packages
pkgs <- c("factoextra", "ggrepel")
for (pkg in pkgs) {
  install_pkg(pkg)
}

# Check STRINGdb dependencies
install_pkg("KernSmooth")
install_pkg("gplots")
install_pkg("plotrix")
install_pkg("RColorBrewer")
install_pkg("sqldf")
install_pkg("igraph")

cat("\n=== Checking installation ===\n")
for (pkg in c("factoextra", "ggrepel", "igraph", "ggplot2", "yaml")) {
  status <- if(requireNamespace(pkg, quietly=TRUE)) "OK" else "MISSING"
  cat(pkg, ":", status, "\n")
}
