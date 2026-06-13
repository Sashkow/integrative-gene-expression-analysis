options(bitmapType = "cairo")

script_dir <- "scripts/integrative_analysis/phase2b_direct_merge"
source(file.path(script_dir, "imputation.R"))
library(yaml)

config_path <- file.path(script_dir,
  "config_yehor_sashko",
  "config_phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko.yaml")
config <- read_yaml(config_path)

output_dir <- "output/article_validation/mean_imputation_benchmark"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load expression data (same as run_phase2b.R) ---

cat("Loading datasets...\n")
exprs_list <- list()
for (ds in config$files$datasets) {
  fname <- config$files$file_map[[ds]]
  if (is.null(fname)) fname <- paste0(ds, config$files$suffix, ".tsv")
  fpath <- file.path(config$paths$mapped_data, fname)
  exprs <- read.delim(fpath, row.names = 1, check.names = FALSE)
  exprs_list[[ds]] <- exprs
  cat("  ", ds, ":", nrow(exprs), "genes x", ncol(exprs), "samples\n")
}

# --- Apply sample filter ---

pheno <- read.delim(config$paths$phenodata, stringsAsFactors = FALSE)
col_map <- config$paths$column_map
names(pheno)[names(pheno) == col_map$arraydatafile_exprscolumnnames] <- "sample_id"
names(pheno)[names(pheno) == col_map$secondaryaccession] <- "dataset_id"
names(pheno)[names(pheno) == col_map$`Gestational.Age.Category`] <- "trimester"
names(pheno)[names(pheno) == col_map$Diagnosis] <- "condition"

pheno$condition <- tolower(trimws(pheno$condition))
allowed_conditions <- tolower(config$sample_filter$condition)
allowed_trimesters <- config$sample_filter$trimester
pheno_filt <- pheno[pheno$condition %in% allowed_conditions &
                    pheno$trimester %in% allowed_trimesters, ]
allowed_samples <- pheno_filt$sample_id

for (ds in names(exprs_list)) {
  keep_cols <- intersect(colnames(exprs_list[[ds]]), allowed_samples)
  exprs_list[[ds]] <- exprs_list[[ds]][, keep_cols, drop = FALSE]
}
exprs_list <- exprs_list[sapply(exprs_list, ncol) > 0]
cat("After filter:", length(exprs_list), "datasets,",
    sum(sapply(exprs_list, ncol)), "samples\n\n")

# --- Run CV for each mask type x method ---

methods <- c("softimpute", "gene_mean", "batch_mean")
mask_types <- c("random_cells", "gene_dataset_block", "gene_dataset_block_hard")
n_repeats <- 5

all_results <- list()

for (mt in mask_types) {
  cat("\n", strrep("=", 60), "\n")
  cat("Mask type:", mt, "\n")
  cat(strrep("=", 60), "\n")

  res <- validate_imputation(
    exprs_list,
    min_datasets = 1L,
    leave_out_fraction = 0.1,
    n_repeats = n_repeats,
    methods = methods,
    mask_type = mt,
    rank_max = 30
  )
  all_results[[mt]] <- res
}

results_df <- do.call(rbind, all_results)
rownames(results_df) <- NULL

write.csv(results_df, file.path(output_dir, "cv_all_methods.csv"),
          row.names = FALSE)

# --- Summary table ---

cat("\n\n", strrep("=", 70), "\n")
cat("COMPARISON TABLE\n")
cat(strrep("=", 70), "\n\n")

summary_rows <- list()
for (mt in mask_types) {
  for (m in methods) {
    sub <- results_df[results_df$mask_type == mt & results_df$method == m &
                      results_df$converged, ]
    if (nrow(sub) == 0) {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        mask_type = mt, method = m,
        r_mean = NA, r_sd = NA,
        rmse_mean = NA, rmse_sd = NA,
        mae_mean = NA, mae_sd = NA,
        n_converged = 0L, n_total = nrow(results_df[results_df$mask_type == mt & results_df$method == m, ]),
        stringsAsFactors = FALSE)
      next
    }
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      mask_type = mt, method = m,
      r_mean = mean(sub$correlation), r_sd = sd(sub$correlation),
      rmse_mean = mean(sub$rmse), rmse_sd = sd(sub$rmse),
      mae_mean = mean(sub$mae), mae_sd = sd(sub$mae),
      n_converged = nrow(sub),
      n_total = nrow(results_df[results_df$mask_type == mt & results_df$method == m, ]),
      stringsAsFactors = FALSE)
  }
}
summary_df <- do.call(rbind, summary_rows)

cat(sprintf("%-30s %-12s %12s %12s %12s %5s\n",
            "Mask", "Method", "r", "RMSE", "MAE", "n"))
cat(strrep("-", 85), "\n")
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  cat(sprintf("%-30s %-12s %5.3f±%.3f %5.3f±%.3f %5.3f±%.3f %3d/%d\n",
              r$mask_type, r$method,
              r$r_mean, r$r_sd, r$rmse_mean, r$rmse_sd,
              r$mae_mean, r$mae_sd, r$n_converged, r$n_total))
}

write.csv(summary_df, file.path(output_dir, "cv_summary.csv"),
          row.names = FALSE)
cat("\nAll outputs saved to:", output_dir, "\n")
