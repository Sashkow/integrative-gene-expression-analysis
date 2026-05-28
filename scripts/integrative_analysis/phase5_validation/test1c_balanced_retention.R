#!/usr/bin/env Rscript
#
# Test 1c: Retention of balanced-reference DEGs across subsamples
# For each subsample size, measures what fraction of the balanced
# 2-dataset (GSE100051+GSE9984) DEGs are retained in the subsample.
# Retention = |balanced_sig ∩ subsample_sig| / |balanced_sig|
#
# Usage: Rscript test1c_balanced_retention.R [--config=config_validation.yaml]

library(parallel)

script_dir <- if (length(grep("--file=", commandArgs(FALSE), value = TRUE)) > 0) {
  dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
} else {
  "scripts/integrative_analysis/phase5_validation"
}
source(file.path(script_dir, "subsampling_helpers.R"))

default_config <- file.path(script_dir, "config_validation.yaml")
config <- parse_config_arg(default_config)
n_iter <- parse_int_arg("n_iter", config$test1$n_iter %||% 30)
n_cores <- parse_int_arg("n_cores", config$parallel$n_cores %||% 10)
if (n_cores == 0) n_cores <- max(1, detectCores() - 2)

output_dir <- config$paths$output
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

bal_de_path <- config$paths$balanced_reference_de
if (is.null(bal_de_path) || !file.exists(bal_de_path))
  stop("balanced_reference_de not set or file not found: ", bal_de_path)

cat("=== Test 1c: Balanced DEG retention ===\n")
cat("Iterations per size:", n_iter, " Cores:", n_cores, "\n")

data <- load_pipeline_data(config)

bal_de <- read.delim(bal_de_path, stringsAsFactors = FALSE)
fdr_thresh <- config$thresholds$fdr %||% 0.05
logfc_thresh <- config$thresholds$logfc %||% 1.0

bal_sig <- bal_de$gene[
  !is.na(bal_de$adj.P.Val) & !is.na(bal_de$logFC) &
  bal_de$adj.P.Val < fdr_thresh & abs(bal_de$logFC) > logfc_thresh
]
bal_up <- bal_de$gene[
  !is.na(bal_de$adj.P.Val) & !is.na(bal_de$logFC) &
  bal_de$adj.P.Val < fdr_thresh & bal_de$logFC > logfc_thresh
]
bal_down <- bal_de$gene[
  !is.na(bal_de$adj.P.Val) & !is.na(bal_de$logFC) &
  bal_de$adj.P.Val < fdr_thresh & bal_de$logFC < -logfc_thresh
]

cat("Balanced ref DEGs:", length(bal_sig),
    " (up:", length(bal_up), " down:", length(bal_down), ")\n")
cat("Full ref DEGs:", length(data$ref_sig_genes), "\n\n")

sizes <- config$test1$sizes_1st_trim
first_pool <- data$first_trim_samples
second_pool <- data$second_trim_samples

cat("1st-trim pool:", length(first_pool),
    " 2nd-trim pool:", length(second_pool), "\n")
cat("Sizes:", paste(sizes, collapse = ", "), "\n\n")

all_rows <- list()

for (N in sizes) {
  if (N > length(first_pool)) {
    cat("Skipping N=", N, "\n")
    next
  }
  cat("--- N_1st =", N, "---\n")
  t0 <- Sys.time()

  iter_results <- mclapply(seq_len(n_iter), function(iter) {
    set.seed(iter * 1000 + N)
    sampled_1st <- sample(first_pool, N, replace = FALSE)
    keep <- c(sampled_1st, second_pool)
    sub <- subsample_exprs_list(data$exprs_list,
                                data$phenodata, keep)

    sink(tempfile())
    on.exit(sink(), add = TRUE)
    res <- run_lean_pipeline(sub$exprs_list, sub$phenodata,
                              config, data$ref_imputed)
    sink()

    if (res$status != "ok") {
      return(data.frame(
        N_1st = N, N_total = N + length(second_pool),
        iter = iter,
        n_deg = NA, n_up = NA, n_down = NA,
        retention_all = NA,
        retention_up = NA,
        retention_down = NA,
        false_pos_rate = NA,
        same_direction_pct = NA,
        logfc_r_matched = NA,
        status = res$status,
        stringsAsFactors = FALSE
      ))
    }

    sub_sig <- res$sig_genes
    sub_de <- res$de_results
    sub_sig_up <- sub_de$gene[
      !is.na(sub_de$adj.P.Val) & !is.na(sub_de$logFC) &
      sub_de$adj.P.Val < fdr_thresh & sub_de$logFC > logfc_thresh
    ]
    sub_sig_down <- sub_de$gene[
      !is.na(sub_de$adj.P.Val) & !is.na(sub_de$logFC) &
      sub_de$adj.P.Val < fdr_thresh & sub_de$logFC < -logfc_thresh
    ]

    retained_all <- intersect(sub_sig, bal_sig)
    retained_up <- intersect(sub_sig_up, bal_up)
    retained_down <- intersect(sub_sig_down, bal_down)
    novel <- setdiff(sub_sig, bal_sig)

    retention_all <- length(retained_all) / max(length(bal_sig), 1)
    retention_up <- length(retained_up) / max(length(bal_up), 1)
    retention_down <- length(retained_down) / max(length(bal_down), 1)
    false_pos_rate <- length(novel) / max(length(sub_sig), 1)

    shared_genes <- intersect(sub_de$gene, bal_de$gene)
    matched_bal <- bal_sig[bal_sig %in% shared_genes]
    if (length(matched_bal) > 1) {
      sub_lfc <- sub_de$logFC[match(matched_bal, sub_de$gene)]
      bal_lfc <- bal_de$logFC[match(matched_bal, bal_de$gene)]
      same_dir <- sum(sign(sub_lfc) == sign(bal_lfc), na.rm = TRUE) /
        length(matched_bal) * 100
      logfc_r <- cor(sub_lfc, bal_lfc, use = "complete.obs")
    } else {
      same_dir <- NA_real_
      logfc_r <- NA_real_
    }

    data.frame(
      N_1st = N, N_total = N + length(second_pool),
      iter = iter,
      n_deg = res$n_sig,
      n_up = res$n_up,
      n_down = res$n_down,
      retention_all = retention_all,
      retention_up = retention_up,
      retention_down = retention_down,
      false_pos_rate = false_pos_rate,
      same_direction_pct = same_dir,
      logfc_r_matched = logfc_r,
      status = res$status,
      stringsAsFactors = FALSE
    )
  }, mc.cores = n_cores)

  rows <- do.call(rbind, iter_results)
  all_rows[[length(all_rows) + 1]] <- rows

  ok <- sum(rows$status == "ok", na.rm = TRUE)
  elapsed <- round(as.numeric(difftime(Sys.time(), t0,
                                        units = "secs")), 1)
  cat(sprintf(
    "  ok: %d/%d  ret_all: %.3f  ret_up: %.3f  ret_down: %.3f  (%.1fs)\n",
    ok, n_iter,
    median(rows$retention_all, na.rm = TRUE),
    median(rows$retention_up, na.rm = TRUE),
    median(rows$retention_down, na.rm = TRUE),
    elapsed
  ))
}

results_df <- do.call(rbind, all_rows)
out_file <- file.path(output_dir, "test1c_balanced_retention.tsv")
write.table(results_df, out_file, sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("\nSaved:", out_file, "\n")
cat("Test 1c complete.\n")
