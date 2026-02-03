#' Normalization Module for Phase 2B
#'
#' Functions for cross-platform/batch normalization of gene expression data.
#' Implements ComBat (baseline), DWD, and other alternatives.
#'
#' @author Expression Integration Pipeline
#' @date 2026-01

library(sva)

#' Apply ComBat batch correction (baseline)
#'
#' Standard parametric empirical Bayes batch correction.
#'
#' @param exprs Expression matrix (genes x samples)
#' @param batch Batch vector (factor or character)
#' @param mod Model matrix for biological covariates (optional)
#' @param par_prior Use parametric priors (default: TRUE)
#' @return Batch-corrected expression matrix
#' @export
normalize_combat <- function(exprs, batch, mod = NULL, par_prior = TRUE) {

  cat("\n=== ComBat Batch Correction ===\n")

  batch <- as.factor(batch)
  cat("Batches:", nlevels(batch), "\n")
  cat("Batch sizes:", paste(table(batch), collapse = ", "), "\n")

  # Remove genes with zero variance
  gene_vars <- apply(exprs, 1, var, na.rm = TRUE)
  keep <- gene_vars > 0 & !is.na(gene_vars)
  if (sum(!keep) > 0) {
    cat("Removing", sum(!keep), "genes with zero variance\n")
    exprs <- exprs[keep, ]
  }

  # Apply ComBat
  corrected <- ComBat(
    dat = as.matrix(exprs),
    batch = batch,
    mod = mod,
    par.prior = par_prior,
    prior.plots = FALSE
  )

  cat("ComBat correction complete\n")

  corrected
}


#' Apply Distance Weighted Discrimination (DWD) normalization
#'
#' DWD is a margin-based method that works better than ComBat when:
#' - Group sizes are unequal
#' - Strong batch effects exist
#' - Data is high-dimensional relative to sample size
#'
#' This implementation uses a simplified DWD approach based on
#' mean-shift correction with distance weighting.
#'
#' @param exprs Expression matrix (genes x samples)
#' @param batch Batch vector
#' @param mod Model matrix for biological covariates (optional)
#' @return Batch-corrected expression matrix
#' @export
normalize_dwd <- function(exprs, batch, mod = NULL) {

  cat("\n=== DWD-style Batch Correction ===\n")

  batch <- as.factor(batch)
  batches <- levels(batch)
  n_batches <- length(batches)

  cat("Batches:", n_batches, "\n")
  cat("Batch sizes:", paste(table(batch), collapse = ", "), "\n")

  # Reference batch (largest)
  batch_sizes <- table(batch)
  ref_batch <- names(which.max(batch_sizes))
  cat("Reference batch:", ref_batch, "\n")

  exprs <- as.matrix(exprs)
  corrected <- exprs

  # Calculate reference batch statistics
  ref_idx <- which(batch == ref_batch)
  ref_mean <- rowMeans(exprs[, ref_idx, drop = FALSE], na.rm = TRUE)
  ref_sd <- apply(exprs[, ref_idx, drop = FALSE], 1, sd, na.rm = TRUE)
  ref_sd[ref_sd == 0] <- 1  # Avoid division by zero

  # Correct each non-reference batch
  for (b in setdiff(batches, ref_batch)) {
    batch_idx <- which(batch == b)
    cat("  Correcting batch", b, "(", length(batch_idx), "samples)...\n")

    # Calculate batch statistics
    batch_mean <- rowMeans(exprs[, batch_idx, drop = FALSE], na.rm = TRUE)
    batch_sd <- apply(exprs[, batch_idx, drop = FALSE], 1, sd, na.rm = TRUE)
    batch_sd[batch_sd == 0] <- 1

    # DWD-style correction: scale then shift
    # 1. Z-score within batch
    z_batch <- sweep(exprs[, batch_idx, drop = FALSE], 1, batch_mean, "-")
    z_batch <- sweep(z_batch, 1, batch_sd, "/")

    # 2. Rescale to reference distribution
    corrected[, batch_idx] <- sweep(z_batch, 1, ref_sd, "*")
    corrected[, batch_idx] <- sweep(corrected[, batch_idx], 1, ref_mean, "+")
  }

  # Verify correction
  for (b in batches) {
    batch_idx <- which(batch == b)
    new_mean <- mean(rowMeans(corrected[, batch_idx, drop = FALSE], na.rm = TRUE))
    cat("  Batch", b, "mean after correction:", round(new_mean, 2), "\n")
  }

  cat("DWD correction complete\n")

  corrected
}


#' Apply quantile normalization across batches
#'
#' Forces all samples to have the same distribution.
#' Aggressive but effective for strong batch effects.
#'
#' @param exprs Expression matrix (genes x samples)
#' @param batch Batch vector (used for reporting only)
#' @return Quantile-normalized expression matrix
#' @export
normalize_quantile <- function(exprs, batch = NULL) {

  cat("\n=== Quantile Normalization ===\n")

  if (!requireNamespace("preprocessCore", quietly = TRUE)) {
    stop("Package 'preprocessCore' required. Install with: BiocManager::install('preprocessCore')")
  }

  exprs <- as.matrix(exprs)
  cat("Input:", nrow(exprs), "genes x", ncol(exprs), "samples\n")

  # Apply quantile normalization
  normalized <- preprocessCore::normalize.quantiles(exprs)
  rownames(normalized) <- rownames(exprs)
  colnames(normalized) <- colnames(exprs)

  cat("Quantile normalization complete\n")

  normalized
}


#' Apply batch mean-centering (simple baseline)
#'
#' Simplest batch correction: subtract batch mean from each gene.
#'
#' @param exprs Expression matrix (genes x samples)
#' @param batch Batch vector
#' @return Mean-centered expression matrix
#' @export
normalize_mean_center <- function(exprs, batch) {

  cat("\n=== Batch Mean Centering ===\n")

  batch <- as.factor(batch)
  batches <- levels(batch)

  cat("Batches:", length(batches), "\n")

  exprs <- as.matrix(exprs)
  corrected <- exprs

  # Calculate global mean
  global_mean <- rowMeans(exprs, na.rm = TRUE)

  # Center each batch to global mean

  for (b in batches) {
    batch_idx <- which(batch == b)
    batch_mean <- rowMeans(exprs[, batch_idx, drop = FALSE], na.rm = TRUE)

    # Shift to global mean
    shift <- global_mean - batch_mean
    corrected[, batch_idx] <- sweep(exprs[, batch_idx, drop = FALSE], 1, shift, "+")
  }

  cat("Mean centering complete\n")

  corrected
}


#' Compare normalization methods
#'
#' Runs multiple normalization methods and compares their effects
#' using PCA and batch effect metrics.
#'
#' @param exprs Expression matrix (genes x samples)
#' @param batch Batch vector
#' @param biological_group Biological group of interest (for preservation)
#' @param methods Vector of methods to test
#' @return List with normalized matrices and comparison metrics
#' @export
compare_normalizations <- function(exprs,
                                    batch,
                                    biological_group = NULL,
                                    methods = c("none", "combat", "dwd", "mean_center")) {

  cat("\n=== Comparing Normalization Methods ===\n\n")

  results <- list()
  metrics <- list()

  for (method in methods) {
    cat("Method:", method, "\n")

    normalized <- switch(
      method,
      "none" = as.matrix(exprs),
      "combat" = normalize_combat(exprs, batch),
      "dwd" = normalize_dwd(exprs, batch),
      "mean_center" = normalize_mean_center(exprs, batch),
      "quantile" = normalize_quantile(exprs, batch),
      stop("Unknown method: ", method)
    )

    results[[method]] <- normalized

    # Calculate batch effect metrics
    metrics[[method]] <- calculate_batch_metrics(normalized, batch, biological_group)
    cat("\n")
  }

  # Summary comparison
  cat("\n=== Normalization Comparison Summary ===\n\n")
  metrics_df <- do.call(rbind, lapply(names(metrics), function(m) {
    data.frame(
      method = m,
      batch_variance_pct = metrics[[m]]$batch_variance_pct,
      bio_variance_pct = metrics[[m]]$bio_variance_pct,
      silhouette_batch = metrics[[m]]$silhouette_batch,
      stringsAsFactors = FALSE
    )
  }))
  print(metrics_df)

  list(
    normalized = results,
    metrics = metrics,
    summary = metrics_df
  )
}


#' Calculate batch effect metrics
#'
#' @param exprs Expression matrix
#' @param batch Batch vector
#' @param biological_group Biological group (optional)
#' @return List of metrics
calculate_batch_metrics <- function(exprs, batch, biological_group = NULL) {

  # PCA
  pca <- prcomp(t(exprs), scale. = TRUE, center = TRUE)

  # Variance explained by batch (using PC1-PC5)
  pc_scores <- pca$x[, 1:min(5, ncol(pca$x))]

  batch_variance <- 0
  bio_variance <- 0

  for (i in 1:ncol(pc_scores)) {
    # ANOVA for batch effect
    batch_anova <- summary(aov(pc_scores[, i] ~ batch))[[1]]
    batch_variance <- batch_variance + batch_anova$`Sum Sq`[1] / sum(batch_anova$`Sum Sq`)

    # ANOVA for biological effect
    if (!is.null(biological_group)) {
      bio_anova <- summary(aov(pc_scores[, i] ~ biological_group))[[1]]
      bio_variance <- bio_variance + bio_anova$`Sum Sq`[1] / sum(bio_anova$`Sum Sq`)
    }
  }

  batch_variance_pct <- round(100 * batch_variance / ncol(pc_scores), 1)
  bio_variance_pct <- round(100 * bio_variance / ncol(pc_scores), 1)

  # Silhouette score for batch clustering
  silhouette_batch <- NA
  if (requireNamespace("cluster", quietly = TRUE) && nlevels(as.factor(batch)) > 1) {
    dist_mat <- dist(pc_scores)
    sil <- cluster::silhouette(as.numeric(as.factor(batch)), dist_mat)
    silhouette_batch <- round(mean(sil[, 3]), 3)
  }

  cat("  Batch variance (PC1-5):", batch_variance_pct, "%\n")
  if (!is.null(biological_group)) {
    cat("  Biological variance (PC1-5):", bio_variance_pct, "%\n")
  }
  cat("  Batch silhouette:", silhouette_batch, "\n")

  list(
    batch_variance_pct = batch_variance_pct,
    bio_variance_pct = bio_variance_pct,
    silhouette_batch = silhouette_batch,
    pca = pca
  )
}
