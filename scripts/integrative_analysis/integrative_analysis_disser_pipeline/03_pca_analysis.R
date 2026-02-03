#' PCA Analysis Module
#'
#' Functions for principal component analysis and visualization
#'
#' @author Expression Integration Pipeline
#' @date 2025

#' Perform PCA on expression data
#'
#' @param exprs Expression matrix (genes x samples)
#' @param center Center data (default: TRUE)
#' @param scale Scale data (default: FALSE)
#' @return PCA object from prcomp
#' @export
perform_pca <- function(exprs, center = TRUE, scale = FALSE) {

  library(stats)

  # Remove genes with missing values
  exprs_clean <- na.omit(exprs)

  n_removed <- nrow(exprs) - nrow(exprs_clean)
  if (n_removed > 0) {
    message("Removed ", n_removed, " genes with missing values")
  }

  # Transpose for PCA (samples x genes)
  pca <- prcomp(t(exprs_clean), center = center, scale. = scale)

  # Calculate variance explained
  n_pcs <- min(10, ncol(pca$x))
  var_explained <- summary(pca)$importance[2, 1:n_pcs] * 100
  message("Variance explained by first ", n_pcs, " PCs: ",
          paste(round(var_explained, 1), collapse = "%, "), "%")

  return(pca)
}


#' Create PCA plots colored by metadata variables
#'
#' @param pca PCA object
#' @param pdata Phenodata
#' @param variables Vector of variable names to color by
#' @param output_dir Directory to save plots
#' @param prefix Plot filename prefix
#' @param axes Principal components to plot (default: c(1,2))
#' @export
plot_pca_metadata <- function(pca, pdata, variables,
                               output_dir, prefix = "pca",
                               axes = c(1, 2)) {

  library(factoextra)
  library(ggplot2)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  for (var in variables) {
    if (!var %in% colnames(pdata)) {
      warning("Variable '", var, "' not found in phenodata")
      next
    }

    # Get grouping variable and convert to factor
    groups <- as.factor(pdata[[var]])

    # Skip if all NA
    if (all(is.na(groups))) {
      warning("Variable '", var, "' contains only NA values")
      next
    }

    # Verify sample alignment
    if (length(groups) != nrow(pca$x)) {
      warning("Sample count mismatch for '", var, "': ", length(groups),
              " groups vs ", nrow(pca$x), " PCA samples")
      next
    }

    tryCatch({
      # Create plot
      p <- fviz_pca_ind(
        pca,
        axes = axes,
        label = "none",
        habillage = groups,
        repel = TRUE,
        palette = "jco",
        addEllipses = TRUE,
        ellipse.type = "confidence",
        title = paste("PCA colored by", var)
      )

      # Save plot
      filename <- file.path(output_dir,
                            paste0(prefix, "_", var, "_PC", axes[1], "_PC", axes[2], ".svg"))
      svg(filename, width = 10, height = 10)
      print(p)
      dev.off()

      message("Saved: ", filename)
    }, error = function(e) {
      warning("Failed to create PCA plot for '", var, "': ", e$message)
    })
  }

  invisible(NULL)
}


#' Create scree plot
#'
#' @param pca PCA object
#' @param output_dir Output directory (for SVG)
#' @param output_dir_png Output directory for PNG (optional)
#' @param prefix Filename prefix
#' @param n_pcs Number of PCs to show (default: 20)
#' @export
plot_pca_scree <- function(pca, output_dir, output_dir_png = NULL, prefix = "pca", n_pcs = 20) {

  library(factoextra)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  p <- fviz_eig(pca, addlabels = TRUE, ncp = n_pcs,
                title = "Scree Plot - Variance Explained")

  filename_svg <- file.path(output_dir, paste0(prefix, "_scree.svg"))
  svg(filename_svg, width = 12, height = 8)
  print(p)
  dev.off()

  if (!is.null(output_dir_png)) {
    filename_png <- file.path(output_dir_png, paste0(prefix, "_scree.png"))
  } else {
    filename_png <- file.path(output_dir, paste0(prefix, "_scree.png"))
  }
  png(filename_png, width = 12, height = 8, units = "in", res = 300)
  print(p)
  dev.off()

  message("Saved: ", filename_svg, " and ", filename_png)

  invisible(NULL)
}


#' Create biplot with gene loadings
#'
#' @param pca PCA object
#' @param pdata Phenodata
#' @param color_by Variable to color samples
#' @param n_genes Number of top genes to show (default: 10)
#' @param output_dir Output directory
#' @param prefix Filename prefix
#' @param axes Principal components to plot (default: c(1,2))
#' @export
plot_pca_biplot <- function(pca, pdata, color_by,
                             n_genes = 10, output_dir, prefix = "pca",
                             axes = c(1, 2)) {

  library(factoextra)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  p <- fviz_pca_biplot(
    pca,
    axes = axes,
    label = "var",
    habillage = pdata[[color_by]],
    repel = TRUE,
    palette = "Set2",
    select.var = list(contrib = n_genes),
    addEllipses = TRUE,
    title = paste("PCA Biplot - Top", n_genes, "genes by contribution")
  )

  filename <- file.path(output_dir,
                        paste0(prefix, "_biplot_", color_by, "_PC", axes[1], "_PC", axes[2], ".svg"))
  svg(filename, width = 12, height = 10)
  print(p)
  dev.off()

  message("Saved: ", filename)

  invisible(NULL)
}
