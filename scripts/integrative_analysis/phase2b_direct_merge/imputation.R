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
#' @param min_datasets Minimum number of datasets a gene must appear in (integer, 1..N)
#' @return List with: matrix (combined with NAs), gene_info, sample_info
#' @export
create_incomplete_matrix <- function(exprs_list, min_datasets = 1L) {

  n_datasets <- length(exprs_list)
  cat("Creating incomplete matrix from", n_datasets, "datasets\n")

  # Get all genes across datasets
  all_genes <- unique(unlist(lapply(exprs_list, rownames)))
  cat("  Total unique genes:", length(all_genes), "\n")

  # Calculate gene presence (vectorized for speed)
  # Create a presence matrix: genes x datasets
  gene_sets <- lapply(exprs_list, rownames)
  presence_matrix <- sapply(gene_sets, function(genes) all_genes %in% genes)
  rownames(presence_matrix) <- all_genes

  gene_presence <- rowSums(presence_matrix)
  gene_coverage <- gene_presence / n_datasets

  # Filter genes by minimum dataset count
  genes_to_keep <- names(gene_presence)[gene_presence >= min_datasets]
  cat("  Genes in >=", min_datasets, "of", n_datasets, "datasets:",
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
    min_datasets = min_datasets
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

  # Complete the matrix. softImpute::complete() reads the biScale
  # attributes on X_centered and returns values on the ORIGINAL scale,
  # so no manual un-scaling is needed (and doing so double-corrects
  # and produces values off by row+column centers).
  X_imputed <- complete(X_centered, final_fit)

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


#' Impute missing values using iterative regularized PCA (missMDA)
#'
#' Wraps \code{missMDA::imputePCA}, an EM-style algorithm that
#' alternates between a low-rank PCA fit and re-imputation of the
#' missing cells. Conceptually closest to \code{softImpute} --- both
#' exploit low-rank structure --- but uses a different regularisation
#' scheme and convergence criterion.
#'
#' missMDA expects observations in rows and variables in columns, so we
#' transpose internally (samples become rows, genes become columns).
#'
#' @param incomplete List from create_incomplete_matrix
#' @param ncp Number of principal components used for imputation
#' @param scale Whether to scale variables to unit variance before PCA
#' @param method Either "Regularized" or "EM"
#' @return List with: matrix (imputed), ncp, gene_info, sample_info
#' @export
impute_missmda <- function(incomplete, ncp = 5, scale = TRUE,
                           method = "Regularized") {

  cat("\n=== missMDA imputePCA ===\n")

  if (!requireNamespace("missMDA", quietly = TRUE)) {
    stop("Package 'missMDA' required. Install with: install.packages('missMDA')")
  }

  X <- incomplete$matrix
  cat("Input matrix:", nrow(X), "genes x", ncol(X), "samples\n")
  cat("Missing values:", sum(is.na(X)), "\n")
  cat("ncp =", ncp, " method =", method, "\n")
  if (sum(is.na(X)) == 0L) {
    cat("Nothing to impute (no NAs). Returning input.\n")
    return(list(matrix = X, ncp = ncp,
                gene_info = incomplete$gene_info,
                sample_info = incomplete$sample_info))
  }
  cat("[", format(Sys.time(), "%H:%M:%S"),
      "] imputePCA EM loop starting (silent until convergence)...\n",
      sep = "")

  # imputePCA: observations = rows, variables = columns
  X_samples_by_genes <- t(X)
  t0 <- Sys.time()
  res <- missMDA::imputePCA(
    X_samples_by_genes,
    ncp = ncp,
    scale = scale,
    method = method
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  X_imputed <- t(res$completeObs)
  # Preserve original cells exactly
  X_final <- incomplete$matrix
  X_final[is.na(X_final)] <- X_imputed[is.na(X_final)]

  cat(sprintf("[%s] Imputation complete (%.1fs)\n",
              format(Sys.time(), "%H:%M:%S"), elapsed))

  list(
    matrix = X_final,
    ncp = ncp,
    gene_info = incomplete$gene_info,
    sample_info = incomplete$sample_info
  )
}


#' Impute missing values using sample-space KNN
#'
#' A hand-rolled sample-space KNN imputer. Contrast with
#' \code{impute::impute.knn}, which searches in gene space. Here each
#' missing cell is filled by averaging k nearest SAMPLES (columns)
#' that have an observed value for that gene. Sample-sample distances
#' are computed once on the fully-observed gene subset, standardised
#' per gene so each contributes equally.
#'
#' @param incomplete List from create_incomplete_matrix
#' @param k Number of neighbours
#' @return List with: matrix (imputed), k, gene_info, sample_info
#' @export
impute_sample_knn <- function(incomplete, k = 10) {

  cat("\n=== Sample-space KNN Imputation ===\n")

  X <- incomplete$matrix
  n_genes <- nrow(X)
  n_samples <- ncol(X)
  cat("Input matrix:", n_genes, "genes x", n_samples, "samples\n")
  cat("Missing values:", sum(is.na(X)), "\n")
  cat("k =", k, "\n")

  t0 <- Sys.time()

  na_per_gene <- rowSums(is.na(X))
  basis_idx <- which(na_per_gene == 0)
  cat(sprintf("[%s] Distance basis: %d fully-observed genes\n",
              format(Sys.time(), "%H:%M:%S"), length(basis_idx)))
  if (length(basis_idx) < 10) {
    stop("Too few fully-observed genes (", length(basis_idx),
         ") to form a distance basis")
  }

  basis <- X[basis_idx, , drop = FALSE]
  row_mu <- rowMeans(basis)
  row_sd <- apply(basis, 1, sd)
  row_sd[row_sd == 0 | is.na(row_sd)] <- 1
  basis_std <- (basis - row_mu) / row_sd

  d_mat <- as.matrix(stats::dist(t(basis_std)))
  diag(d_mat) <- Inf
  cat(sprintf("[%s] Sample distance matrix: %dx%d\n",
              format(Sys.time(), "%H:%M:%S"),
              nrow(d_mat), ncol(d_mat)))

  na_genes <- which(na_per_gene > 0)
  n_to_impute <- length(na_genes)
  cat(sprintf("[%s] Genes needing imputation: %d / %d (%.1f%%)\n",
              format(Sys.time(), "%H:%M:%S"),
              n_to_impute, n_genes,
              100 * n_to_impute / n_genes))

  if (n_to_impute == 0) {
    cat("Nothing to impute (no NAs). Returning input.\n")
    return(list(
      matrix = X,
      k = k,
      gene_info = incomplete$gene_info,
      sample_info = incomplete$sample_info
    ))
  }

  X_imputed <- X
  progress_every <- max(1L, floor(n_to_impute / 20))
  for (idx in seq_along(na_genes)) {
    g <- na_genes[idx]
    row_vals <- X[g, ]
    missing_samples <- which(is.na(row_vals))
    observed_samples <- which(!is.na(row_vals))
    if (length(observed_samples) == 0L) next
    for (s in missing_samples) {
      dists <- d_mat[s, observed_samples]
      ord <- order(dists)
      nn <- observed_samples[ord[seq_len(min(k, length(ord)))]]
      X_imputed[g, s] <- mean(row_vals[nn])
    }
    if (idx %% progress_every == 0L || idx == n_to_impute) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      pct <- 100 * idx / n_to_impute
      rate <- idx / max(elapsed, 1e-6)
      eta <- (n_to_impute - idx) / max(rate, 1e-6)
      cat(sprintf("[%s] %d/%d genes (%.0f%%)  elapsed=%.1fs  ETA=%.1fs\n",
                  format(Sys.time(), "%H:%M:%S"),
                  idx, n_to_impute, pct, elapsed, eta))
    }
  }

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("[%s] Imputation complete (%.1fs)\n",
              format(Sys.time(), "%H:%M:%S"), elapsed))

  list(
    matrix = X_imputed,
    k = k,
    gene_info = incomplete$gene_info,
    sample_info = incomplete$sample_info
  )
}


#' Validate imputation accuracy using leave-out cross-validation
#'
#' Masks a fraction of observed values, imputes them, and measures
#' correlation with true values. Two masking strategies are supported:
#'
#' \itemize{
#'   \item \code{"random_cells"} — uniform random draw of individual
#'     matrix cells (the classical in-distribution leave-out).
#'   \item \code{"gene_dataset_block"} — mask all cells for a
#'     randomly chosen (gene, dataset) pair, stopping once the total
#'     number of masked cells is at least
#'     \code{leave_out_fraction * n_observed}. This simulates the
#'     real missingness pattern where a gene is absent from one
#'     entire dataset because the chip did not probe it.
#' }
#'
#' @param exprs_list Named list of expression matrices
#' @param min_datasets Minimum number of datasets a gene must appear in
#' @param leave_out_fraction Fraction of values to mask (0-1)
#' @param n_repeats Number of CV repetitions
#' @param methods Character vector of methods to test: "softimpute", "knn"
#' @param rank_max Passed to impute_softimpute
#' @param k Passed to impute_knn
#' @param mask_type Either "random_cells" or "gene_dataset_block"
#' @return Data frame with validation metrics
#' @export
validate_imputation <- function(exprs_list,
                                 min_datasets = 1L,
                                 leave_out_fraction = 0.1,
                                 n_repeats = 5,
                                 methods = c("softimpute"),
                                 rank_max = 50,
                                 k = 10,
                                 mask_type = "random_cells",
                                 min_obs_per_gene = 4L) {

  mask_type <- match.arg(mask_type, c("random_cells", "gene_dataset_block"))

  cat("\n=== Imputation Validation ===\n")
  cat("Mask type:", mask_type, "\n")
  cat("Leave-out fraction:", leave_out_fraction, "\n")
  cat("Repeats:", n_repeats, "\n")
  cat("Methods:", paste(methods, collapse = ", "), "\n\n")

  # Create the base incomplete matrix
  incomplete <- create_incomplete_matrix(exprs_list, min_datasets)
  X <- incomplete$matrix

  # Get indices of observed (non-NA) values
  observed_idx <- which(!is.na(X))
  n_observed <- length(observed_idx)
  n_mask_target <- round(n_observed * leave_out_fraction)

  cat("Observed values:", n_observed, "\n")
  cat("Target masked cells per repeat:", n_mask_target, "\n")

  # Pre-compute (gene, dataset) pair pool for block masking
  obs_pairs <- NULL
  ds_cols_cache <- NULL
  row_idx_map <- NULL
  if (mask_type == "gene_dataset_block") {
    sample_info <- incomplete$sample_info
    datasets <- unique(sample_info$dataset)
    ds_cols_cache <- split(seq_len(ncol(X)), sample_info$dataset)
    row_idx_map <- setNames(seq_len(nrow(X)), rownames(X))

    # Boolean presence matrix: genes x datasets (TRUE if gene has >=1 obs)
    presence <- sapply(datasets, function(ds) {
      rowSums(!is.na(X[, ds_cols_cache[[ds]], drop = FALSE])) > 0
    })
    colnames(presence) <- datasets
    n_ds_per_gene <- rowSums(presence)

    # Only include pairs where the gene is present in ALL datasets.
    # Masking then removes one dataset's contribution from an otherwise
    # fully-observed row, which keeps the row's non-masked density high
    # enough for softImpute's biScale to converge. Partially-observed
    # genes are excluded because block-masking them can produce rows
    # that are >50% empty, breaking row-scaling.
    n_datasets_total <- length(datasets)
    pair_list <- list()
    for (ds in datasets) {
      eligible <- presence[, ds] & (n_ds_per_gene == n_datasets_total)
      genes_in_ds <- rownames(X)[eligible]
      if (length(genes_in_ds) > 0) {
        pair_list[[ds]] <- data.frame(
          gene = genes_in_ds,
          dataset = ds,
          stringsAsFactors = FALSE
        )
      }
    }
    obs_pairs <- do.call(rbind, pair_list)
    rownames(obs_pairs) <- NULL
    cat("Observed (gene, dataset) pairs (gene in all",
        n_datasets_total, "datasets):", nrow(obs_pairs), "\n")
  }
  cat("\n")

  results <- list()

  for (rep in seq_len(n_repeats)) {
    cat(sprintf("[%s] Repeat %d / %d ...\n",
                format(Sys.time(), "%H:%M:%S"), rep, n_repeats))

    set.seed(rep * 123)  # Reproducible

    if (mask_type == "random_cells") {
      mask_idx <- sample(observed_idx, n_mask_target)
      n_blocks <- NA_integer_
    } else {
      shuffled <- sample.int(nrow(obs_pairs))
      mask_idx_chunks <- vector("list", length(shuffled))
      n_masked_so_far <- 0L
      n_blocks <- 0L
      n_blocks_skipped_floor <- 0L
      gene_obs_remaining <- rowSums(!is.na(X))
      for (p in shuffled) {
        g  <- obs_pairs$gene[p]
        ds <- obs_pairs$dataset[p]
        cols <- ds_cols_cache[[ds]]
        r <- row_idx_map[[g]]
        # linear indices into X
        block_idx <- (cols - 1L) * nrow(X) + r
        block_idx <- block_idx[!is.na(X[block_idx])]
        if (length(block_idx) == 0) next
        # Reject blocks that would push the gene below the per-row floor:
        # softImpute's biScale and ALS fit need at least min_obs_per_gene
        # non-NA cells per row to standardize and converge.
        if (gene_obs_remaining[r] - length(block_idx) < min_obs_per_gene) {
          n_blocks_skipped_floor <- n_blocks_skipped_floor + 1L
          next
        }
        n_blocks <- n_blocks + 1L
        mask_idx_chunks[[n_blocks]] <- block_idx
        gene_obs_remaining[r] <- gene_obs_remaining[r] - length(block_idx)
        n_masked_so_far <- n_masked_so_far + length(block_idx)
        if (n_masked_so_far >= n_mask_target) break
      }
      mask_idx <- unlist(mask_idx_chunks[seq_len(n_blocks)], use.names = FALSE)
      cat("  Blocks selected:", n_blocks,
          " (skipped to keep >=", min_obs_per_gene, "obs/gene:",
          n_blocks_skipped_floor, ")\n")
    }
    n_mask <- length(mask_idx)
    cat("  Cells masked:", n_mask, "\n")

    # Create masked matrix
    X_masked <- X
    true_values <- X[mask_idx]
    X_masked[mask_idx] <- NA

    # Create new incomplete object with masked data
    incomplete_masked <- incomplete
    incomplete_masked$matrix <- X_masked

    # Test each method (both score against the SAME mask_idx / true_values)
    for (method in methods) {
      method_cfg <- list(rank_max = rank_max, k = k)
      imputed <- tryCatch(
        run_imputer(method, incomplete_masked, method_cfg),
        error = function(e) {
          cat("  ", method, " failed: ", conditionMessage(e), "\n", sep = "")
          NULL
        }
      )

      if (is.null(imputed)) {
        results[[length(results) + 1]] <- data.frame(
          method = method,
          mask_type = mask_type,
          repeat_n = rep,
          correlation = NA_real_,
          rmse = NA_real_,
          mae = NA_real_,
          n_masked = n_mask,
          n_blocks = n_blocks,
          stringsAsFactors = FALSE
        )
        next
      }

      # Get imputed values at masked positions
      imputed_values <- imputed$matrix[mask_idx]

      # Calculate metrics
      correlation <- cor(true_values, imputed_values, use = "complete.obs")
      rmse <- sqrt(mean((true_values - imputed_values)^2, na.rm = TRUE))
      mae <- mean(abs(true_values - imputed_values), na.rm = TRUE)

      results[[length(results) + 1]] <- data.frame(
        method = method,
        mask_type = mask_type,
        repeat_n = rep,
        correlation = correlation,
        rmse = rmse,
        mae = mae,
        n_masked = n_mask,
        n_blocks = n_blocks,
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


#' Imputation method registry
#'
#' Maps an imputation method name (as used in configs and CLI args) to a
#' uniform adapter function. Each adapter takes \code{(incomplete, cfg)}
#' where \code{incomplete} is the list returned by
#' \code{create_incomplete_matrix()} and \code{cfg} is a named list of
#' method-specific options (usually \code{config$imputation[[method]]}).
#' Every adapter must return a list with at least a \code{matrix} element
#' containing the imputed gene x sample matrix. Adding a new method
#' means adding one entry here --- the dispatch sites in
#' \code{run_phase2b.R}, \code{validate_imputation()}, and
#' \code{compare_gene_recovery()} do not need to change.
#'
#' @export
IMPUTERS <- list(
  softimpute = function(incomplete, cfg = list()) {
    impute_softimpute(
      incomplete,
      rank_max = cfg$rank_max %||% 50,
      lambda   = cfg$lambda   %||% 0,
      thresh   = cfg$thresh   %||% 1e-5,
      maxit    = cfg$maxit    %||% 100,
      type     = cfg$type     %||% "als"
    )
  },
  knn = function(incomplete, cfg = list()) {
    impute_knn(incomplete, k = cfg$k %||% 10)
  },
  missmda = function(incomplete, cfg = list()) {
    impute_missmda(
      incomplete,
      ncp    = cfg$ncp    %||% 5,
      scale  = cfg$scale  %||% TRUE,
      method = cfg$method %||% "Regularized"
    )
  },
  sample_knn = function(incomplete, cfg = list()) {
    impute_sample_knn(incomplete, k = cfg$k %||% 10)
  }
)

# Local fallback so IMPUTERS adapters work even when run_phase2b.R's
# %||% has not been sourced yet (e.g. when imputation.R is used from
# tests or standalone scripts).
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}


#' Run a registered imputation method by name
#'
#' Thin wrapper around \code{IMPUTERS} that gives a clear error for
#' unknown methods and keeps dispatch sites one-liners.
#'
#' @param method Name of a registered imputation method (e.g. "softimpute").
#' @param incomplete List from \code{create_incomplete_matrix()}.
#' @param cfg Named list of method-specific options. Typically
#'   \code{config$imputation[[method]]}.
#' @return Whatever the adapter returns (at least a list with
#'   \code{matrix}).
#' @export
run_imputer <- function(method, incomplete, cfg = list()) {
  if (!method %in% names(IMPUTERS)) {
    stop("Unknown imputation method: '", method,
         "'. Registered methods: ",
         paste(names(IMPUTERS), collapse = ", "))
  }
  IMPUTERS[[method]](incomplete, cfg %||% list())
}


#' Compare gene recovery across imputation methods and dataset-count thresholds
#'
#' @param exprs_list Named list of expression matrices
#' @param min_datasets_range Integer vector of minimum dataset counts to compare
#'   (defaults to 1:N where N = length(exprs_list))
#' @param methods Vector of methods: "none", "softimpute", "knn"
#' @return Data frame with gene counts
#' @export
compare_gene_recovery <- function(exprs_list,
                                   min_datasets_range = NULL,
                                   methods = c("none", "softimpute")) {

  if (is.null(min_datasets_range))
    min_datasets_range <- seq_len(length(exprs_list))

  cat("\n=== Gene Recovery Comparison ===\n\n")

  results <- list()

  for (min_ds in min_datasets_range) {
    incomplete <- create_incomplete_matrix(exprs_list, min_datasets = min_ds)

    for (method in methods) {
      cat("min_datasets:", min_ds, "Method:", method, "\n")

      mat <- if (method == "none") {
        incomplete$matrix
      } else {
        tryCatch(
          run_imputer(method, incomplete, list())$matrix,
          error = function(e) {
            cat("  ", method, " failed: ", conditionMessage(e), "\n", sep = "")
            NULL
          }
        )
      }

      if (is.null(mat)) next

      n_genes <- nrow(mat)
      n_missing_total <- sum(is.na(mat))
      pct_missing <- round(100 * n_missing_total / length(mat), 1)

      results[[length(results) + 1]] <- data.frame(
        min_datasets = min_ds,
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
