#' Imputation Module for Phase 2B
#'
#' Functions for imputing missing gene expression values across datasets
#' before merging. Implements softImpute (nuclear norm matrix completion)
#' and KNN imputation.
#'
#' @author Expression Integration Pipeline
#' @date 2026-01

library(softImpute)

#' Create combined expression matrix with NAs for missing genes
#'
#' @param exprs_list Named list of expression matrices (genes x samples)
#' @param min_coverage Minimum proportion of datasets a gene must appear in (0-1)
#' @return List with: matrix (combined with NAs), gene_info, sample_info
#' @export
create_incomplete_matrix <- function(exprs_list, min_coverage = 0.5) {

  n_datasets <- length(exprs_list)
  cat("Creating incomplete matrix from", n_datasets, "datasets\n")

  # Get all genes across datasets
  all_genes <- unique(unlist(lapply(exprs_list, rownames)))
  cat("  Total unique genes:", length(all_genes), "\n")

  # Calculate gene coverage (vectorized for speed)
  # Create a presence matrix: genes x datasets
  gene_sets <- lapply(exprs_list, rownames)
  presence_matrix <- sapply(gene_sets, function(genes) all_genes %in% genes)
  rownames(presence_matrix) <- all_genes

  gene_presence <- rowSums(presence_matrix)
  gene_coverage <- gene_presence / n_datasets

  # Filter genes by coverage
  genes_to_keep <- names(gene_coverage)[gene_coverage >= min_coverage]
  cat("  Genes with >=", round(min_coverage * 100), "% coverage:",
      length(genes_to_keep), "\n")

  # Get all samples
  all_samples <- unlist(lapply(exprs_list, colnames))
  cat("  Total samples:", length(all_samples), "\n")

  # Create sample-to-dataset mapping
  sample_dataset <- rep(names(exprs_list), sapply(exprs_list, ncol))
  names(sample_dataset) <- all_samples

  # Create incomplete matrix
  combined <- matrix(
    NA,
    nrow = length(genes_to_keep),
    ncol = length(all_samples),
    dimnames = list(genes_to_keep, all_samples)
  )

  # Fill in values
  for (ds_name in names(exprs_list)) {
    exprs <- exprs_list[[ds_name]]
    common_genes <- intersect(genes_to_keep, rownames(exprs))
    combined[common_genes, colnames(exprs)] <- as.matrix(exprs[common_genes, ])
  }

  # Calculate missingness
  n_missing <- sum(is.na(combined))
  pct_missing <- round(100 * n_missing / length(combined), 1)
  cat("  Missing values:", n_missing, "(", pct_missing, "%)\n")

  # Gene info
  gene_info <- data.frame(
    gene = genes_to_keep,
    n_datasets = gene_presence[genes_to_keep],
    coverage = gene_coverage[genes_to_keep],
    n_missing = rowSums(is.na(combined)),
    stringsAsFactors = FALSE
  )

  # Sample info
  sample_info <- data.frame(
    sample = all_samples,
    dataset = sample_dataset[all_samples],
    n_missing = colSums(is.na(combined)),
    stringsAsFactors = FALSE
  )

  list(
    matrix = combined,
    gene_info = gene_info,
    sample_info = sample_info,
    n_datasets = n_datasets,
    coverage_threshold = min_coverage
  )
}


#' Impute missing values using softImpute
#'
#' Uses nuclear norm regularization to complete the matrix,
#' exploiting low-rank structure typical of gene expression data.
#'
#' @param incomplete List from create_incomplete_matrix
#' @param rank_max Maximum rank for approximation
#' @param lambda Regularization parameter (0 = auto via warm starts)
#' @param thresh Convergence threshold
#' @param maxit Maximum iterations
#' @param type Algorithm type: "als" or "svd"
#' @return List with: matrix (imputed), fit (softImpute object), validation
#' @export
impute_softimpute <- function(incomplete,
                               rank_max = 50,
                               lambda = 0,
                               thresh = 1e-5,
                               maxit = 100,
                               type = "als") {

  cat("\n=== softImpute Matrix Completion ===\n")

  X <- incomplete$matrix
  cat("Input matrix:", nrow(X), "genes x", ncol(X), "samples\n")
  cat("Missing values:", sum(is.na(X)), "/", length(X),
      "(", round(100 * sum(is.na(X)) / length(X), 1), "%)\n")

  # Center and scale for better convergence (reduced iterations for speed)
  X_centered <- biScale(X, maxit = 20, thresh = 1e-3)

  # Auto-select lambda via warm starts if lambda = 0
  if (lambda == 0) {
    cat("Auto-selecting lambda via warm starts...\n")

    # Get lambda sequence (8 values for speed, was 20)
    lambda_max <- lambda0(X_centered)
    lambda_seq <- exp(seq(log(lambda_max), log(lambda_max * 0.01), length.out = 8))

    # Warm start path
    fits <- list()
    fit <- NULL
    for (i in seq_along(lambda_seq)) {
      lam <- lambda_seq[i]
      fit <- softImpute(
        X_centered,
        rank.max = rank_max,
        lambda = lam,
        type = type,
        thresh = thresh,
        maxit = maxit,
        warm.start = fit
      )
      fits[[i]] <- fit

      # Check rank
      if (!is.null(fit$d) && length(fit$d) >= rank_max * 0.8) {
        cat("  Lambda", round(lam, 4), ": rank =", length(fit$d),
            "(approaching rank_max)\n")
      }
    }

    # Use last fit (smallest lambda, highest rank)
    final_fit <- fit
    lambda_used <- lambda_seq[length(lambda_seq)]
    cat("Using lambda =", round(lambda_used, 4), "\n")

  } else {
    cat("Using specified lambda =", lambda, "\n")
    final_fit <- softImpute(
      X_centered,
      rank.max = rank_max,
      lambda = lambda,
      type = type,
      thresh = thresh,
      maxit = maxit
    )
    lambda_used <- lambda
  }

  # Complete the matrix
  X_complete <- complete(X_centered, final_fit)

  # Unscale
  X_imputed <- X_complete * attr(X_centered, "biScale:column")$scale
  X_imputed <- sweep(X_imputed, 2, attr(X_centered, "biScale:column")$center, "+")
  X_imputed <- sweep(X_imputed, 1, attr(X_centered, "biScale:row")$center, "+")

  # Preserve original values (only impute NAs)
  X_final <- incomplete$matrix
  X_final[is.na(X_final)] <- X_imputed[is.na(X_final)]

  # Report
  final_rank <- if (!is.null(final_fit$d)) length(final_fit$d) else NA
  cat("\nImputation complete:\n")
  cat("  Final rank:", final_rank, "\n")
  cat("  Imputed values:", sum(is.na(incomplete$matrix)), "\n")

  list(
    matrix = X_final,
    fit = final_fit,
    lambda = lambda_used,
    rank = final_rank,
    gene_info = incomplete$gene_info,
    sample_info = incomplete$sample_info
  )
}


#' Impute missing values using KNN
#'
#' Simple k-nearest neighbors imputation based on gene similarity.
#'
#' @param incomplete List from create_incomplete_matrix
#' @param k Number of neighbors
#' @return List with: matrix (imputed), validation
#' @export
impute_knn <- function(incomplete, k = 10) {

  cat("\n=== KNN Imputation ===\n")

  if (!requireNamespace("impute", quietly = TRUE)) {
    stop("Package 'impute' required. Install with: BiocManager::install('impute')")
  }

  X <- incomplete$matrix
  cat("Input matrix:", nrow(X), "genes x", ncol(X), "samples\n")
  cat("Missing values:", sum(is.na(X)), "\n")

  # impute.knn expects genes as rows
  X_imputed <- impute::impute.knn(X, k = k)$data

  cat("Imputation complete\n")

  list(
    matrix = X_imputed,
    k = k,
    gene_info = incomplete$gene_info,
    sample_info = incomplete$sample_info
  )
}


#' Validate imputation accuracy using leave-out cross-validation
#'
#' Randomly masks a fraction of observed values, imputes them,
#' and measures correlation with true values.
#'
#' @param exprs_list Named list of expression matrices
#' @param min_coverage Minimum gene coverage for merging
#' @param leave_out_fraction Fraction of values to mask (0-1)
#' @param n_repeats Number of CV repetitions
#' @param methods Character vector of methods to test: "softimpute", "knn"
#' @param ... Additional parameters passed to imputation functions
#' @return Data frame with validation metrics
#' @export
validate_imputation <- function(exprs_list,
                                 min_coverage = 0.5,
                                 leave_out_fraction = 0.1,
                                 n_repeats = 5,
                                 methods = c("softimpute"),
                                 ...) {

  cat("\n=== Imputation Validation ===\n")
  cat("Leave-out fraction:", leave_out_fraction, "\n")
  cat("Repeats:", n_repeats, "\n")
  cat("Methods:", paste(methods, collapse = ", "), "\n\n")

  # Create the base incomplete matrix
  incomplete <- create_incomplete_matrix(exprs_list, min_coverage)
  X <- incomplete$matrix

  # Get indices of observed (non-NA) values
  observed_idx <- which(!is.na(X))
  n_observed <- length(observed_idx)
  n_mask <- round(n_observed * leave_out_fraction)

  cat("Observed values:", n_observed, "\n")
  cat("Values to mask per repeat:", n_mask, "\n\n")

  results <- list()

  for (rep in seq_len(n_repeats)) {
    cat("Repeat", rep, "/", n_repeats, "...\n")

    # Randomly select values to mask
    set.seed(rep * 123)  # Reproducible
    mask_idx <- sample(observed_idx, n_mask)

    # Create masked matrix
    X_masked <- X
    true_values <- X[mask_idx]
    X_masked[mask_idx] <- NA

    # Create new incomplete object with masked data
    incomplete_masked <- incomplete
    incomplete_masked$matrix <- X_masked

    # Test each method
    for (method in methods) {
      imputed <- switch(
        method,
        "softimpute" = impute_softimpute(incomplete_masked, ...),
        "knn" = impute_knn(incomplete_masked, ...),
        stop("Unknown method: ", method)
      )

      # Get imputed values at masked positions
      imputed_values <- imputed$matrix[mask_idx]

      # Calculate metrics
      correlation <- cor(true_values, imputed_values, use = "complete.obs")
      rmse <- sqrt(mean((true_values - imputed_values)^2, na.rm = TRUE))
      mae <- mean(abs(true_values - imputed_values), na.rm = TRUE)

      results[[length(results) + 1]] <- data.frame(
        method = method,
        repeat_n = rep,
        correlation = correlation,
        rmse = rmse,
        mae = mae,
        n_masked = n_mask,
        stringsAsFactors = FALSE
      )
    }
  }

  results_df <- do.call(rbind, results)

  # Summary statistics
  cat("\n=== Validation Summary ===\n")
  summary_df <- aggregate(
    cbind(correlation, rmse, mae) ~ method,
    data = results_df,
    FUN = function(x) c(mean = mean(x), sd = sd(x))
  )

  for (method in unique(results_df$method)) {
    method_results <- results_df[results_df$method == method, ]
    cat("\n", method, ":\n", sep = "")
    cat("  Correlation: ", round(mean(method_results$correlation), 3),
        " (+/- ", round(sd(method_results$correlation), 3), ")\n", sep = "")
    cat("  RMSE: ", round(mean(method_results$rmse), 3),
        " (+/- ", round(sd(method_results$rmse), 3), ")\n", sep = "")
    cat("  MAE: ", round(mean(method_results$mae), 3),
        " (+/- ", round(sd(method_results$mae), 3), ")\n", sep = "")
  }

  results_df
}


#' Compare gene recovery across imputation methods and thresholds
#'
#' @param exprs_list Named list of expression matrices
#' @param thresholds Vector of coverage thresholds to compare
#' @param methods Vector of methods: "none", "softimpute", "knn"
#' @return Data frame with gene counts
#' @export
compare_gene_recovery <- function(exprs_list,
                                   thresholds = c(1.0, 0.75, 0.5, 0.25),
                                   methods = c("none", "softimpute")) {

  cat("\n=== Gene Recovery Comparison ===\n\n")

  results <- list()

  for (thresh in thresholds) {
    for (method in methods) {
      cat("Threshold:", thresh, "Method:", method, "\n")

      incomplete <- create_incomplete_matrix(exprs_list, min_coverage = thresh)
      n_genes <- nrow(incomplete$matrix)
      n_missing_total <- sum(is.na(incomplete$matrix))
      pct_missing <- round(100 * n_missing_total / length(incomplete$matrix), 1)

      results[[length(results) + 1]] <- data.frame(
        threshold = thresh,
        method = method,
        n_genes = n_genes,
        n_missing = n_missing_total,
        pct_missing = pct_missing,
        stringsAsFactors = FALSE
      )
    }
  }

  results_df <- do.call(rbind, results)

  cat("\nSummary:\n")
  print(results_df)

  results_df
}
