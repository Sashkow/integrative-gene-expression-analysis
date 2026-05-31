suppressPackageStartupMessages(library(openxlsx))

runs <- list(
  list(
    run_id       = "1_2_restoration_blockmask_imputed_0",
    cohort       = "1_2",
    mask_type    = "gene_dataset_block",
    imputed_weight = 0.0,
    dir          = "output/phase2b_combat/phase2b_1_2_restoration_blockmask_imputed_0"
  ),
  list(
    run_id       = "2_3_restoration_blockmask_imputed_0",
    cohort       = "2_3",
    mask_type    = "gene_dataset_block",
    imputed_weight = 0.0,
    dir          = "output/phase2b_combat/phase2b_2_3_restoration_blockmask_imputed_0"
  ),
  list(
    run_id       = "1_2_2nd_trim_only",
    cohort       = "1_2",
    mask_type    = "gene_dataset_block",
    imputed_weight = 0.0,
    dir          = "output/phase2b_combat/phase2b_1_2_2nd_trim_only"
  ),
  list(
    run_id       = "2_3_2nd_trim_only",
    cohort       = "2_3",
    mask_type    = "gene_dataset_block",
    imputed_weight = 0.0,
    dir          = "output/phase2b_combat/phase2b_2_3_2nd_trim_only"
  ),
  list(
    run_id       = "1_2_all_datasets",
    cohort       = "1_2",
    mask_type    = "gene_dataset_block",
    imputed_weight = 0.0,
    dir          = "output/phase2b_combat/phase2b_1_2_all_datasets"
  ),
  list(
    run_id       = "2_3_all_datasets",
    cohort       = "2_3",
    mask_type    = "gene_dataset_block",
    imputed_weight = 0.0,
    dir          = "output/phase2b_combat/phase2b_2_3_all_datasets"
  ),
  list(
    run_id       = "1_2_balanced",
    cohort       = "1_2",
    mask_type    = "gene_dataset_block",
    imputed_weight = 0.0,
    dir          = "output/phase2b_combat/phase2b_1_2_balanced"
  ),
  list(
    run_id       = "2_3_balanced",
    cohort       = "2_3",
    mask_type    = "gene_dataset_block",
    imputed_weight = 0.0,
    dir          = "output/phase2b_combat/phase2b_2_3_balanced"
  ),
  list(
    run_id       = "1_2_all_datasets_batch_in_limma",
    cohort       = "1_2",
    mask_type    = "gene_dataset_block",
    imputed_weight = 0.0,
    dir          = "output/phase2b_batch_in_limma/phase2b_1_2_all_datasets_batch_in_limma"
  ),
  list(
    run_id       = "2_3_all_datasets_batch_in_limma",
    cohort       = "2_3",
    mask_type    = "gene_dataset_block",
    imputed_weight = 0.0,
    dir          = "output/phase2b_batch_in_limma/phase2b_2_3_all_datasets_batch_in_limma"
  ),
  list(
    run_id       = "1_2_balanced_ruv",
    cohort       = "1_2",
    mask_type    = "gene_dataset_block",
    imputed_weight = 0.0,
    dir          = "output/phase2b_ruv/phase2b_1_2_balanced_ruv"
  ),
  list(
    run_id       = "2_3_balanced_ruv",
    cohort       = "2_3",
    mask_type    = "gene_dataset_block",
    imputed_weight = 0.0,
    dir          = "output/phase2b_ruv/phase2b_2_3_balanced_ruv"
  ),
  list(
    run_id       = "1_2_all_datasets_ruv",
    cohort       = "1_2",
    mask_type    = "gene_dataset_block",
    imputed_weight = 0.0,
    dir          = "output/phase2b_ruv/phase2b_1_2_all_datasets_ruv"
  ),
  list(
    run_id       = "2_3_all_datasets_ruv",
    cohort       = "2_3",
    mask_type    = "gene_dataset_block",
    imputed_weight = 0.0,
    dir          = "output/phase2b_ruv/phase2b_2_3_all_datasets_ruv"
  )
)

fdr_cutoff   <- 0.05
logfc_cutoff <- 1.0
method_to_file_default <- c(
  none       = "difexp_none_combat.tsv",
  softimpute = "difexp_softimpute_combat.tsv",
  knn        = "difexp_knn_combat.tsv",
  missmda    = "difexp_missmda_combat.tsv",
  sample_knn = "difexp_sample_knn_combat.tsv"
)

find_de_file <- function(run_dir, imp_method) {
  candidates <- list.files(run_dir, pattern = paste0("^difexp_", imp_method, "_.*\\.tsv$"))
  if (length(candidates) == 1) return(candidates[1])
  if (length(candidates) > 1) return(candidates[1])
  if (imp_method %in% names(method_to_file_default)) return(method_to_file_default[[imp_method]])
  NA_character_
}

empty_gene_sets <- function() list(
  fdr_logfc1 = character(0),
  fdr_only   = character(0),
  up_only    = character(0),
  down_only  = character(0)
)

load_de_table <- function(tsv_path) {
  if (!file.exists(tsv_path)) return(NULL)
  read.table(tsv_path, header = TRUE, sep = "\t",
             stringsAsFactors = FALSE, check.names = FALSE,
             quote = "")
}

load_gene_sets <- function(tsv_path) {
  d <- load_de_table(tsv_path)
  if (is.null(d)) return(empty_gene_sets())
  genes <- as.character(d$gene)
  fdr_mask   <- !is.na(d$adj.P.Val) & d$adj.P.Val < fdr_cutoff
  logfc_mask <- !is.na(d$logFC) & abs(d$logFC) > logfc_cutoff
  up_mask    <- !is.na(d$logFC) & d$logFC > 0
  down_mask  <- !is.na(d$logFC) & d$logFC < 0
  list(
    fdr_logfc1 = genes[fdr_mask & logfc_mask],
    fdr_only   = genes[fdr_mask],
    up_only    = genes[fdr_mask & up_mask],
    down_only  = genes[fdr_mask & down_mask]
  )
}

fdr_summary <- function(tsv_path) {
  d <- load_de_table(tsv_path)
  if (is.null(d) || !"adj.P.Val" %in% colnames(d)) {
    return(list(median_fdr = NA_real_, mean_neglog10_fdr = NA_real_,
                pct_fdr_below_05 = NA_real_, pct_fdr_below_01 = NA_real_))
  }
  fdr <- d$adj.P.Val[!is.na(d$adj.P.Val)]
  nlog10 <- -log10(pmax(fdr, .Machine$double.xmin))
  list(
    median_fdr         = median(fdr),
    mean_neglog10_fdr  = mean(nlog10),
    pct_fdr_below_05   = 100 * mean(fdr < 0.05),
    pct_fdr_below_01   = 100 * mean(fdr < 0.01)
  )
}

fmt_mean_sd <- function(x, digits = 3) {
  if (length(x) == 0 || all(is.na(x))) return(NA_character_)
  sprintf(paste0("%.", digits, "f \u00B1 %.", digits, "f"),
          mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}

methods_order <- c("none", "softimpute", "knn", "missmda", "sample_knn")

summary_rows <- list()

for (r in runs) {
  val_csv <- file.path(r$dir, "imputation_validation.csv")
  deg_csv <- file.path(r$dir, "method_comparison.csv")

  val <- tryCatch(read.csv(val_csv, stringsAsFactors = FALSE),
                  error = function(e) NULL)
  deg <- tryCatch(read.csv(deg_csv, stringsAsFactors = FALSE),
                  error = function(e) NULL)

  val_methods <- if (!is.null(val)) unique(val$method) else character()
  deg_methods <- if (!is.null(deg)) unique(deg$imputation) else character()
  all_methods <- unique(c(val_methods, deg_methods))
  if (length(all_methods) == 0) next
  all_methods <- all_methods[order(match(all_methods, methods_order))]

  run_gene_sets <- list()
  for (m in all_methods) {
    de_fname <- find_de_file(r$dir, m)
    run_gene_sets[[m]] <- if (!is.na(de_fname)) {
      load_gene_sets(file.path(r$dir, de_fname))
    } else {
      empty_gene_sets()
    }
  }
  none_sets <- if ("none" %in% names(run_gene_sets)) {
    run_gene_sets[["none"]]
  } else {
    NULL
  }

  for (m in all_methods) {
    vm <- if (!is.null(val)) val[val$method == m, , drop = FALSE] else NULL
    dm <- if (!is.null(deg)) deg[deg$imputation == m, , drop = FALSE] else NULL

    m_sets <- run_gene_sets[[m]]

    if (is.null(none_sets)) {
      added_logfc1   <- NA_integer_
      removed_logfc1 <- NA_integer_
      added_fdr      <- NA_integer_
      removed_fdr    <- NA_integer_
    } else {
      added_logfc1   <- length(setdiff(m_sets$fdr_logfc1, none_sets$fdr_logfc1))
      removed_logfc1 <- length(setdiff(none_sets$fdr_logfc1, m_sets$fdr_logfc1))
      added_fdr      <- length(setdiff(m_sets$fdr_only, none_sets$fdr_only))
      removed_fdr    <- length(setdiff(none_sets$fdr_only, m_sets$fdr_only))
    }

    de_fname <- find_de_file(r$dir, m)
    de_tsv <- if (!is.na(de_fname)) file.path(r$dir, de_fname) else ""
    fdr_stats <- fdr_summary(de_tsv)

    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      run_id           = r$run_id,
      cohort           = r$cohort,
      mask_type        = r$mask_type,
      imputed_weight   = r$imputed_weight,
      method           = m,
      n_repeats        = if (!is.null(vm)) nrow(vm) else NA_integer_,
      n_masked_per_rep = if (!is.null(vm) && nrow(vm) > 0) round(mean(vm$n_masked)) else NA_integer_,
      pearson_mean     = if (!is.null(vm) && nrow(vm) > 0) mean(vm$correlation) else NA_real_,
      pearson_sd       = if (!is.null(vm) && nrow(vm) > 0) sd(vm$correlation)   else NA_real_,
      pearson          = if (!is.null(vm) && nrow(vm) > 0) fmt_mean_sd(vm$correlation) else NA_character_,
      rmse_mean        = if (!is.null(vm) && nrow(vm) > 0) mean(vm$rmse) else NA_real_,
      rmse_sd          = if (!is.null(vm) && nrow(vm) > 0) sd(vm$rmse)   else NA_real_,
      rmse             = if (!is.null(vm) && nrow(vm) > 0) fmt_mean_sd(vm$rmse) else NA_character_,
      mae_mean         = if (!is.null(vm) && nrow(vm) > 0) mean(vm$mae) else NA_real_,
      mae_sd           = if (!is.null(vm) && nrow(vm) > 0) sd(vm$mae)   else NA_real_,
      mae              = if (!is.null(vm) && nrow(vm) > 0) fmt_mean_sd(vm$mae) else NA_character_,
      n_genes                 = if (!is.null(dm) && nrow(dm) > 0) dm$n_genes[1]       else NA_integer_,
      n_sig_fdr_logfc1        = if (!is.null(dm) && nrow(dm) > 0) dm$n_significant[1] else NA_integer_,
      n_up_fdr_logfc1         = if (!is.null(dm) && nrow(dm) > 0) dm$n_up[1]          else NA_integer_,
      n_down_fdr_logfc1       = if (!is.null(dm) && nrow(dm) > 0) dm$n_down[1]        else NA_integer_,
      n_added_fdr_logfc1      = added_logfc1,
      n_removed_fdr_logfc1    = removed_logfc1,
      n_sig_fdr_only          = length(m_sets$fdr_only),
      n_up_fdr_only           = length(m_sets$up_only),
      n_down_fdr_only         = length(m_sets$down_only),
      n_added_fdr_only        = added_fdr,
      n_removed_fdr_only      = removed_fdr,
      median_fdr              = fdr_stats$median_fdr,
      mean_neglog10_fdr       = fdr_stats$mean_neglog10_fdr,
      pct_fdr_below_05        = fdr_stats$pct_fdr_below_05,
      pct_fdr_below_01        = fdr_stats$pct_fdr_below_01,
      stringsAsFactors = FALSE
    )
  }
}

summary_df <- do.call(rbind, summary_rows)
rownames(summary_df) <- NULL

cat("Summary rows assembled:\n")
print(summary_df[, c("run_id", "method", "pearson",
                     "n_sig_fdr_logfc1",
                     "n_added_fdr_logfc1", "n_removed_fdr_logfc1",
                     "n_sig_fdr_only",
                     "n_added_fdr_only", "n_removed_fdr_only")])

legend_df <- data.frame(
  Column = c(
    "run_id", "cohort", "mask_type", "imputed_weight", "method",
    "n_repeats", "n_masked_per_rep",
    "pearson_mean/sd/pearson", "rmse_mean/sd/rmse", "mae_mean/sd/mae",
    "n_genes",
    "n_sig_fdr_logfc1", "n_up_fdr_logfc1", "n_down_fdr_logfc1",
    "n_added_fdr_logfc1", "n_removed_fdr_logfc1",
    "n_sig_fdr_only", "n_up_fdr_only", "n_down_fdr_only",
    "n_added_fdr_only", "n_removed_fdr_only",
    "median_fdr", "mean_neglog10_fdr",
    "pct_fdr_below_05", "pct_fdr_below_01"
  ),
  Meaning = c(
    "Pipeline run identifier (matches output/phase2b_<run_id>/ directory).",
    "Trimester contrast: 1_2 = First vs Second, 2_3 = Second vs Term.",
    "Held-out mask strategy: random_cells draws uniform cells; gene_dataset_block masks every cell for a (gene, dataset) pair.",
    "Weight given to imputed cells in limma fit (1.0 = full, 0.0 = hard-mask from DE).",
    "Imputation method: none / softimpute / knn / missmda / sample_knn (all followed by ComBat).",
    "Number of independent leave-out repeats aggregated.",
    "Mean number of held-out cells per repeat.",
    "Pearson correlation between imputed and ground-truth held-out cells (mean, SD, and formatted 'mean \u00B1 SD').",
    "Root-mean-square error in log2 units (mean, SD, formatted).",
    "Mean absolute error in log2 units (mean, SD, formatted).",
    "Number of protein-coding genes entering the DE model.",
    "Significant DEGs at adj.P.Val < 0.05 AND |logFC| > 1 (article headline cutoff, from method_comparison.csv).",
    "Up-regulated subset of n_sig_fdr_logfc1.",
    "Down-regulated subset of n_sig_fdr_logfc1.",
    "DEGs gained by this method vs 'none' in the same run (FDR+|logFC|>1 cutoff): |method_set \\ none_set|.",
    "DEGs lost by this method vs 'none' in the same run (FDR+|logFC|>1 cutoff): |none_set \\ method_set|.",
    "Significant DEGs at adj.P.Val < 0.05 only (no logFC cutoff), recomputed from difexp_<method>_combat.tsv.",
    "Up-regulated (logFC > 0) subset of n_sig_fdr_only.",
    "Down-regulated (logFC < 0) subset of n_sig_fdr_only.",
    "DEGs gained by this method vs 'none' in the same run (FDR-only cutoff).",
    "DEGs lost by this method vs 'none' in the same run (FDR-only cutoff).",
    "Median adj.P.Val across ALL genes (lower = globally more significant).",
    "Mean -log10(adj.P.Val) across ALL genes (higher = globally more significant).",
    "Percentage of all genes with adj.P.Val < 0.05.",
    "Percentage of all genes with adj.P.Val < 0.01 (stricter cutoff to gauge depth of significance)."
  ),
  stringsAsFactors = FALSE
)

wb <- createWorkbook()
addWorksheet(wb, "runs_summary")
writeData(wb, "runs_summary", summary_df, withFilter = TRUE)
freezePane(wb, "runs_summary", firstActiveRow = 2)
setColWidths(wb, "runs_summary",
             cols = seq_len(ncol(summary_df)),
             widths = "auto")

addWorksheet(wb, "legend")
writeData(wb, "legend", legend_df)
setColWidths(wb, "legend", cols = 1:2, widths = c(26, 110))

out_path <- "articles/imputation_article/phase2b_runs_summary.xlsx"
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
saveWorkbook(wb, out_path, overwrite = TRUE)

cat(sprintf("\nWrote: %s (%d rows)\n", out_path, nrow(summary_df)))
