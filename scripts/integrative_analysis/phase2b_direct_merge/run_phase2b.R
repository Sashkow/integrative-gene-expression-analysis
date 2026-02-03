#!/usr/bin/env Rscript

#' Run Phase 2B: Direct Merge Improvements
#'
#' This phase evaluates imputation and alternative normalization methods
#' to improve gene recovery when directly merging expression datasets.
#'
#' Usage:
#'   Rscript run_phase2b.R [options]
#'
#' Options:
#'   --config=PATH         Path to config YAML (default: config_phase2b.yaml)
#'   --imputation=METHOD   Run only specific imputation: softimpute, knn, none
#'   --normalization=METHOD Run only specific normalization: combat, dwd, mean_center
#'   --validate_only       Only run imputation validation (no DE analysis)
#'   --output_dir=PATH     Override output directory
#'   --no_archive          Don't archive previous results
#'   --help                Show this help message

library(yaml)
library(limma)

# Null coalescing operator
`%||%` <- function(a, b) if (is.null(a)) b else a

# Get script directory
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  return("scripts/integrative_analysis/phase2b_direct_merge")
}

script_dir <- get_script_dir()

# Source utilities
utils_path <- "scripts/utils/logging_utils.R"
if (!file.exists(utils_path)) {
  utils_path <- file.path("..", "..", "utils", "logging_utils.R")
}
if (file.exists(utils_path)) {
  source(utils_path)
} else {
  # Fallback logging functions
  setup_logging <- function(dir, prefix = "log") {
    log_file <- file.path(dir, paste0(prefix, "_", format(Sys.time(), "%Y-%m-%d_%H%M%S"), ".txt"))
    cat("Logging to:", log_file, "\n")
    log_file
  }
  close_logging <- function(log_file) invisible(NULL)
  archive_previous_results <- function(output_dir) {
    if (dir.exists(output_dir) && length(list.files(output_dir)) > 0) {
      archive_dir <- file.path(output_dir, "archive",
                               format(Sys.time(), "%Y-%m-%d_%H%M%S"))
      dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
      files <- list.files(output_dir, full.names = TRUE)
      files <- files[!grepl("archive", files)]
      if (length(files) > 0) {
        file.copy(files, archive_dir, recursive = TRUE)
        unlink(files, recursive = TRUE)
        cat("Archived previous results to:", archive_dir, "\n")
      }
    }
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
}

# Source module files
source(file.path(script_dir, "imputation.R"))
source(file.path(script_dir, "normalization.R"))

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Defaults
config_path <- file.path(script_dir, "config_phase2b.yaml")
override_imputation <- NULL
override_normalization <- NULL
validate_only <- FALSE
override_output <- NULL
no_archive <- FALSE

# Parse arguments
for (arg in args) {
  if (arg == "--help" || arg == "-h") {
    cat("
Phase 2B: Direct Merge Improvements

Usage:
  Rscript run_phase2b.R [options]

Options:
  --config=PATH           Path to config YAML (default: config_phase2b.yaml)
  --imputation=METHOD     Run only specific imputation: softimpute, knn, none
  --normalization=METHOD  Run only specific normalization: combat, dwd, mean_center
  --validate_only         Only run imputation validation (no DE analysis)
  --output_dir=PATH       Override output directory
  --no_archive            Don't archive previous results
  --help                  Show this help message

Examples:
  Rscript run_phase2b.R
  Rscript run_phase2b.R --imputation=softimpute --normalization=dwd
  Rscript run_phase2b.R --validate_only
")
    quit(status = 0)
  } else if (grepl("^--config=", arg)) {
    config_path <- sub("^--config=", "", arg)
  } else if (grepl("^--imputation=", arg)) {
    override_imputation <- sub("^--imputation=", "", arg)
  } else if (grepl("^--normalization=", arg)) {
    override_normalization <- sub("^--normalization=", "", arg)
  } else if (arg == "--validate_only") {
    validate_only <- TRUE
  } else if (grepl("^--output_dir=", arg)) {
    override_output <- sub("^--output_dir=", "", arg)
  } else if (arg == "--no_archive") {
    no_archive <- TRUE
  }
}

# Load config
if (!file.exists(config_path)) {
  stop("Config file not found: ", config_path)
}
config <- yaml::read_yaml(config_path)

# Apply overrides
if (!is.null(override_output)) config$paths$output <- override_output

output_dir <- config$paths$output

cat("\n")
cat("============================================================\n")
cat("  Phase 2B: Direct Merge Improvements\n")
cat("============================================================\n\n")

cat("Configuration:\n")
cat("  Config file:     ", config_path, "\n")
cat("  Mapped path:     ", config$paths$mapped_data, "\n")
cat("  Phenodata:       ", config$paths$phenodata, "\n")
cat("  Output dir:      ", output_dir, "\n")
cat("  Comparison:      ", config$phenotype$contrast, "vs", config$phenotype$baseline, "\n")
cat("  Coverage thresh: ", config$coverage$threshold, "\n")
cat("\n")

# Archive previous results
if (!no_archive) {
  archive_previous_results(output_dir)
} else {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}

# Setup logging
log_file <- setup_logging(output_dir, prefix = "phase2b_log")

# ============================================================
# Step 1: Load datasets
# ============================================================

cat("=== Step 1: Loading Datasets ===\n\n")

# Load phenodata
phenodata <- read.csv(config$paths$phenodata, stringsAsFactors = FALSE)
cat("Loaded phenodata:", nrow(phenodata), "samples\n")

# Load expression files
datasets <- config$files$datasets
exprs_list <- list()

for (ds in datasets) {
  file_path <- file.path(config$paths$mapped_data, paste0(ds, ".tsv"))
  if (file.exists(file_path)) {
    exprs <- read.delim(file_path, row.names = 1, check.names = FALSE)
    exprs_list[[ds]] <- exprs
    cat("  Loaded", ds, ":", nrow(exprs), "genes x", ncol(exprs), "samples\n")
  } else {
    cat("  WARNING: File not found:", file_path, "\n")
  }
}

cat("\nLoaded", length(exprs_list), "datasets\n")

# ============================================================
# Step 2: Gene Coverage Analysis
# ============================================================

cat("\n=== Step 2: Gene Coverage Analysis ===\n")

gene_recovery <- compare_gene_recovery(
  exprs_list,
  thresholds = config$coverage$compare_thresholds %||% c(1.0, 0.75, 0.5, 0.25)
)

# Save gene recovery report
write.csv(gene_recovery, file.path(output_dir, "gene_recovery_comparison.csv"),
          row.names = FALSE)

# ============================================================
# Step 3: Imputation Validation (if enabled)
# ============================================================

validation_results <- NULL

if (config$validation$leave_out_fraction %||% 0 > 0) {
  cat("\n=== Step 3: Imputation Validation ===\n")

  methods_to_validate <- c()
  if (config$imputation$softimpute$enabled %||% FALSE) {
    methods_to_validate <- c(methods_to_validate, "softimpute")
  }

  if (length(methods_to_validate) > 0) {
    validation_results <- validate_imputation(
      exprs_list,
      min_coverage = config$coverage$threshold,
      leave_out_fraction = config$validation$leave_out_fraction %||% 0.1,
      n_repeats = config$validation$n_repeats %||% 5,
      methods = methods_to_validate,
      rank_max = config$imputation$softimpute$rank_max %||% 50
    )

    # Save validation results
    write.csv(validation_results,
              file.path(output_dir, "imputation_validation.csv"),
              row.names = FALSE)
    cat("\nSaved imputation validation results\n")
  }
}

if (validate_only) {
  cat("\n--validate_only flag set, stopping here.\n")
  close_logging(log_file)
  quit(status = 0)
}

# ============================================================
# Step 4: Create Merged Matrices with Different Methods
# ============================================================

cat("\n=== Step 4: Creating Merged Expression Matrices ===\n")

# Methods to run
imputation_methods <- c()
if (config$imputation$none$enabled %||% TRUE) imputation_methods <- c(imputation_methods, "none")
if (config$imputation$softimpute$enabled %||% FALSE) imputation_methods <- c(imputation_methods, "softimpute")
if (config$imputation$knn$enabled %||% FALSE) imputation_methods <- c(imputation_methods, "knn")

if (!is.null(override_imputation)) {
  imputation_methods <- override_imputation
}

normalization_methods <- c()
if (config$normalization$combat$enabled %||% TRUE) normalization_methods <- c(normalization_methods, "combat")
if (config$normalization$dwd$enabled %||% FALSE) normalization_methods <- c(normalization_methods, "dwd")

if (!is.null(override_normalization)) {
  normalization_methods <- override_normalization
}

# Create incomplete matrix
incomplete <- create_incomplete_matrix(
  exprs_list,
  min_coverage = config$coverage$threshold
)

# Store all results for comparison
all_results <- list()

for (imp_method in imputation_methods) {
  cat("\n--- Imputation:", imp_method, "---\n")

  # Apply imputation
  if (imp_method == "none") {
    # Inner join (no imputation)
    common_genes <- Reduce(intersect, lapply(exprs_list, rownames))
    cat("Common genes across all datasets:", length(common_genes), "\n")

    # Use unname to prevent list names from being prepended to column names
    merged_exprs <- do.call(cbind, unname(lapply(exprs_list, function(e) e[common_genes, ])))
    imputed <- list(matrix = merged_exprs, method = "none")
  } else if (imp_method == "softimpute") {
    imputed <- impute_softimpute(
      incomplete,
      rank_max = config$imputation$softimpute$rank_max %||% 50,
      lambda = config$imputation$softimpute$lambda %||% 0,
      thresh = config$imputation$softimpute$thresh %||% 1e-5,
      maxit = config$imputation$softimpute$maxit %||% 100
    )
    merged_exprs <- imputed$matrix
  } else if (imp_method == "knn") {
    imputed <- impute_knn(incomplete, k = config$imputation$knn$k %||% 10)
    merged_exprs <- imputed$matrix
  }

  # Save imputed matrix
  if (config$output$save_imputed_matrices %||% FALSE) {
    imputed_file <- file.path(output_dir, paste0("exprs_imputed_", imp_method, ".tsv"))
    write.table(merged_exprs, imputed_file, sep = "\t", quote = FALSE)
    cat("Saved imputed matrix:", imputed_file, "\n")
  }

  # Create phenodata for merged samples
  merged_samples <- colnames(merged_exprs)
  merged_pdata <- phenodata[phenodata$arraydatafile_exprscolumnnames %in% merged_samples, ]

  # Match order
  rownames(merged_pdata) <- merged_pdata$arraydatafile_exprscolumnnames
  merged_pdata <- merged_pdata[merged_samples, ]

  # Create batch variable
  batch <- as.factor(merged_pdata$secondaryaccession)

  # Create biological group
  bio_group <- merged_pdata[[config$phenotype$group_column]]

  for (norm_method in normalization_methods) {
    cat("\n  Normalization:", norm_method, "\n")

    result_key <- paste0(imp_method, "_", norm_method)

    # Apply normalization
    normalized <- switch(
      norm_method,
      "combat" = normalize_combat(merged_exprs, batch,
                                   mod = model.matrix(~bio_group)),
      "dwd" = normalize_dwd(merged_exprs, batch),
      "mean_center" = normalize_mean_center(merged_exprs, batch),
      merged_exprs
    )

    # Save normalized matrix
    if (config$output$save_normalized_matrices %||% FALSE) {
      norm_file <- file.path(output_dir, paste0("exprs_", result_key, ".tsv"))
      write.table(normalized, norm_file, sep = "\t", quote = FALSE)
      cat("  Saved normalized matrix:", norm_file, "\n")
    }

    # ============================================================
    # Step 5: Differential Expression Analysis
    # ============================================================

    cat("\n  Running differential expression...\n")

    # Filter to comparison groups
    keep_samples <- merged_pdata[[config$phenotype$group_column]] %in%
      c(config$phenotype$baseline, config$phenotype$contrast)

    de_exprs <- normalized[, keep_samples]
    de_pdata <- merged_pdata[keep_samples, ]

    # Create design matrix
    group <- factor(de_pdata[[config$phenotype$group_column]],
                    levels = c(config$phenotype$baseline, config$phenotype$contrast))
    design <- model.matrix(~group)

    # Fit linear model
    fit <- lmFit(de_exprs, design)
    fit <- eBayes(fit)

    # Get results
    de_results <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
    de_results$gene <- rownames(de_results)
    de_results <- de_results[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

    # Apply thresholds
    fdr_thresh <- config$thresholds$fdr %||% 0.05
    logfc_thresh <- config$thresholds$logfc %||% 1.0

    de_significant <- de_results[
      de_results$adj.P.Val < fdr_thresh & abs(de_results$logFC) > logfc_thresh,
    ]

    cat("  Total genes:", nrow(de_results), "\n")
    cat("  Significant (FDR <", fdr_thresh, ", |logFC| >", logfc_thresh, "):",
        nrow(de_significant), "\n")
    cat("    Up-regulated:", sum(de_significant$logFC > 0), "\n")
    cat("    Down-regulated:", sum(de_significant$logFC < 0), "\n")

    # Save DE results
    de_file <- file.path(output_dir, paste0("difexp_", result_key, ".tsv"))
    write.table(de_results, de_file, sep = "\t", row.names = FALSE, quote = FALSE)

    sig_file <- file.path(output_dir, paste0("difexp_significant_", result_key, ".tsv"))
    write.table(de_significant, sig_file, sep = "\t", row.names = FALSE, quote = FALSE)

    # Store results for comparison
    all_results[[result_key]] <- list(
      imputation = imp_method,
      normalization = norm_method,
      n_genes = nrow(de_results),
      n_significant = nrow(de_significant),
      n_up = sum(de_significant$logFC > 0),
      n_down = sum(de_significant$logFC < 0),
      de_results = de_results,
      de_significant = de_significant
    )
  }
}

# ============================================================
# Step 6: Method Comparison
# ============================================================

cat("\n=== Step 6: Method Comparison ===\n\n")

# Summary table
comparison_df <- do.call(rbind, lapply(names(all_results), function(key) {
  r <- all_results[[key]]
  data.frame(
    method = key,
    imputation = r$imputation,
    normalization = r$normalization,
    n_genes = r$n_genes,
    n_significant = r$n_significant,
    n_up = r$n_up,
    n_down = r$n_down,
    stringsAsFactors = FALSE
  )
}))

print(comparison_df)

# Save comparison
write.csv(comparison_df, file.path(output_dir, "method_comparison.csv"),
          row.names = FALSE)

# Calculate overlaps between methods
if (length(all_results) > 1) {
  cat("\n=== DEG Overlap Between Methods ===\n\n")

  method_names <- names(all_results)
  overlap_matrix <- matrix(NA, nrow = length(method_names), ncol = length(method_names),
                           dimnames = list(method_names, method_names))

  for (i in seq_along(method_names)) {
    for (j in seq_along(method_names)) {
      genes_i <- all_results[[method_names[i]]]$de_significant$gene
      genes_j <- all_results[[method_names[j]]]$de_significant$gene
      overlap_matrix[i, j] <- length(intersect(genes_i, genes_j))
    }
  }

  cat("Overlap matrix (significant genes):\n")
  print(overlap_matrix)

  write.csv(overlap_matrix, file.path(output_dir, "deg_overlap_matrix.csv"))
}

# ============================================================
# Save Summary
# ============================================================

summary_file <- file.path(output_dir, "summary.txt")
summary_conn <- file(summary_file, "w")
writeLines(c(
  "Phase 2B: Direct Merge Improvements Summary",
  "==========================================",
  "",
  paste("Date:", Sys.time()),
  paste("Config:", config_path),
  "",
  paste("Comparison:", config$phenotype$contrast, "vs", config$phenotype$baseline),
  paste("Datasets:", paste(datasets, collapse = ", ")),
  paste("Coverage threshold:", config$coverage$threshold),
  "",
  "Gene Recovery:",
  capture.output(print(gene_recovery)),
  "",
  "Method Comparison:",
  capture.output(print(comparison_df)),
  "",
  if (!is.null(validation_results)) {
    c("Imputation Validation:",
      capture.output(print(aggregate(correlation ~ method, validation_results, mean))))
  } else {
    "Imputation validation: not run"
  }
), summary_conn)
close(summary_conn)

cat("\nSaved summary to:", summary_file, "\n")

# Close logging
close_logging(log_file)

cat("\n============================================================\n")
cat("  Phase 2B Complete!\n")
cat("============================================================\n\n")

cat("Output files in", output_dir, ":\n")
for (f in list.files(output_dir, pattern = "\\.(tsv|csv|txt)$")) {
  cat("  ", f, "\n")
}
cat("\n")
