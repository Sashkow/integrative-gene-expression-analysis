#!/usr/bin/env Rscript

library(ggplot2)
library(gridExtra)

result_key <- "none_combat_ref"

runs <- list(
  balanced         = "output/dissertation/validation_runs/cv_soncin_mikheev/balanced",
  soncin_1st_only  = "output/dissertation/validation_runs/cv_soncin_mikheev/soncin_1st_only",
  soncin_2nd_only  = "output/dissertation/validation_runs/cv_soncin_mikheev/soncin_2nd_only",
  mikheev_1st_only = "output/dissertation/validation_runs/cv_soncin_mikheev/mikheev_1st_only",
  mikheev_2nd_only = "output/dissertation/validation_runs/cv_soncin_mikheev/mikheev_2nd_only"
)

output_dir <- "output/dissertation/validation_runs/cv_soncin_mikheev/comparison"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------
# Load DE results and expression matrices
# ------------------------------------------------------------------
de_list   <- list()
exprs_list <- list()

for (run_name in names(runs)) {
  de_file <- file.path(runs[[run_name]], paste0("difexp_", result_key, ".tsv"))
  if (!file.exists(de_file)) {
    cat("WARNING: missing", de_file, "- skipping\n")
    next
  }
  de_list[[run_name]] <- read.delim(de_file, stringsAsFactors = FALSE)

  exprs_file <- file.path(runs[[run_name]], paste0("exprs_", result_key, ".tsv"))
  if (file.exists(exprs_file)) {
    exprs_list[[run_name]] <- read.delim(exprs_file, row.names = 1, check.names = FALSE)
  }
}

if (!"balanced" %in% names(de_list)) stop("Balanced reference run not found")

ref_de <- de_list[["balanced"]]
fdr_thresh  <- 0.05
logfc_thresh <- 1.0
ref_sig <- ref_de$gene[ref_de$adj.P.Val < fdr_thresh & abs(ref_de$logFC) > logfc_thresh]

cat("\n============================================================\n")
cat("  Cross-Validation: Soncin x Mikheev Batch-Biology Confounding\n")
cat("============================================================\n\n")
cat("Reference (balanced): ", nrow(ref_de), " genes tested, ",
    length(ref_sig), " DEGs\n\n", sep = "")

# ------------------------------------------------------------------
# 1. DEG overlap analysis
# ------------------------------------------------------------------
cat("=== 1. DEG Overlap with Balanced Reference ===\n\n")

overlap_rows <- list()

for (run_name in setdiff(names(de_list), "balanced")) {
  cv_de  <- de_list[[run_name]]
  cv_sig <- cv_de$gene[cv_de$adj.P.Val < fdr_thresh & abs(cv_de$logFC) > logfc_thresh]

  shared   <- intersect(ref_sig, cv_sig)
  lost     <- setdiff(ref_sig, cv_sig)
  gained   <- setdiff(cv_sig, ref_sig)
  jaccard  <- length(shared) / length(union(ref_sig, cv_sig))

  cat(sprintf("  %-20s: %3d DEGs | shared=%3d  lost=%3d  gained=%3d  Jaccard=%.3f\n",
              run_name, length(cv_sig), length(shared), length(lost), length(gained), jaccard))

  overlap_rows[[run_name]] <- data.frame(
    run       = run_name,
    n_deg     = length(cv_sig),
    shared    = length(shared),
    lost      = length(lost),
    gained    = length(gained),
    jaccard   = jaccard,
    stringsAsFactors = FALSE
  )
}
cat("\n")

overlap_df <- do.call(rbind, overlap_rows)
write.csv(overlap_df, file.path(output_dir, "deg_overlap.csv"), row.names = FALSE)

# ------------------------------------------------------------------
# 2. logFC correlation: how much does the effect size shift?
# ------------------------------------------------------------------
cat("=== 2. logFC Correlation with Balanced Reference ===\n\n")

logfc_rows <- list()
logfc_plots <- list()

for (run_name in setdiff(names(de_list), "balanced")) {
  cv_de <- de_list[[run_name]]
  common_genes <- intersect(ref_de$gene, cv_de$gene)

  ref_lfc <- ref_de$logFC[match(common_genes, ref_de$gene)]
  cv_lfc  <- cv_de$logFC[match(common_genes, cv_de$gene)]

  r <- cor(ref_lfc, cv_lfc, use = "complete.obs")
  rmse <- sqrt(mean((ref_lfc - cv_lfc)^2, na.rm = TRUE))
  mean_shift <- mean(cv_lfc - ref_lfc, na.rm = TRUE)

  cat(sprintf("  %-20s: r=%.4f  RMSE=%.4f  mean_shift=%.4f  (%d genes)\n",
              run_name, r, rmse, mean_shift, length(common_genes)))

  logfc_rows[[run_name]] <- data.frame(
    run = run_name, pearson_r = r, rmse = rmse,
    mean_shift = mean_shift, n_genes = length(common_genes),
    stringsAsFactors = FALSE
  )

  plot_df <- data.frame(balanced = ref_lfc, cv = cv_lfc)
  is_deg <- common_genes %in% ref_sig
  plot_df$deg <- ifelse(is_deg, "DEG in balanced", "not DEG")

  p <- ggplot(plot_df, aes(x = balanced, y = cv, colour = deg)) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = c("DEG in balanced" = "red", "not DEG" = "grey60")) +
    labs(title = run_name, x = "logFC (balanced)", y = paste0("logFC (", run_name, ")")) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")
  logfc_plots[[run_name]] <- p
}
cat("\n")

logfc_df <- do.call(rbind, logfc_rows)
write.csv(logfc_df, file.path(output_dir, "logfc_correlation.csv"), row.names = FALSE)

if (length(logfc_plots) > 0) {
  pdf(file.path(output_dir, "logfc_scatter.pdf"), width = 12, height = 10)
  grid.arrange(grobs = logfc_plots, ncol = 2)
  dev.off()
}

# ------------------------------------------------------------------
# 3. Biology absorbed as batch effect
# ------------------------------------------------------------------
cat("=== 3. Biology Absorbed as Batch Effect ===\n\n")
cat("Genes significant in balanced but NOT in leave-out run suggest\n")
cat("that removing one trimester from a dataset caused ComBat to\n")
cat("absorb real biology as batch effect.\n\n")

absorbed_list <- list()

for (run_name in setdiff(names(de_list), "balanced")) {
  cv_de  <- de_list[[run_name]]
  cv_sig <- cv_de$gene[cv_de$adj.P.Val < fdr_thresh & abs(cv_de$logFC) > logfc_thresh]
  lost   <- setdiff(ref_sig, cv_sig)

  if (length(lost) == 0) {
    cat(sprintf("  %-20s: no lost DEGs\n", run_name))
    next
  }

  lost_de <- cv_de[cv_de$gene %in% lost, ]
  lost_ref <- ref_de[ref_de$gene %in% lost, ]
  merged_lost <- merge(
    lost_ref[, c("gene", "logFC", "adj.P.Val")],
    lost_de[, c("gene", "logFC", "adj.P.Val")],
    by = "gene", suffixes = c(".balanced", paste0(".", run_name))
  )
  merged_lost$logFC_attenuation <- abs(merged_lost$logFC.balanced) - abs(merged_lost[[paste0("logFC.", run_name)]])

  n_attenuated <- sum(merged_lost$logFC_attenuation > 0)
  n_flipped    <- sum(sign(merged_lost$logFC.balanced) != sign(merged_lost[[paste0("logFC.", run_name)]]))
  mean_atten   <- mean(merged_lost$logFC_attenuation)

  cat(sprintf("  %-20s: %3d lost DEGs | %d attenuated, %d sign-flipped, mean attenuation=%.3f\n",
              run_name, nrow(merged_lost), n_attenuated, n_flipped, mean_atten))

  absorbed_list[[run_name]] <- merged_lost
  write.csv(merged_lost,
            file.path(output_dir, paste0("absorbed_as_batch_", run_name, ".csv")),
            row.names = FALSE)
}
cat("\n")

# ------------------------------------------------------------------
# 4. Batch effect leaking as biology
# ------------------------------------------------------------------
cat("=== 4. Batch Effect Leaking as False Biology ===\n\n")
cat("Genes significant in leave-out but NOT in balanced suggest\n")
cat("that confounding introduced false positives.\n\n")

leaked_list <- list()

for (run_name in setdiff(names(de_list), "balanced")) {
  cv_de  <- de_list[[run_name]]
  cv_sig <- cv_de$gene[cv_de$adj.P.Val < fdr_thresh & abs(cv_de$logFC) > logfc_thresh]
  gained <- setdiff(cv_sig, ref_sig)

  if (length(gained) == 0) {
    cat(sprintf("  %-20s: no gained DEGs\n", run_name))
    next
  }

  gained_de  <- cv_de[cv_de$gene %in% gained, ]
  gained_ref <- ref_de[ref_de$gene %in% gained, ]
  merged_gained <- merge(
    gained_ref[, c("gene", "logFC", "adj.P.Val")],
    gained_de[, c("gene", "logFC", "adj.P.Val")],
    by = "gene", suffixes = c(".balanced", paste0(".", run_name))
  )

  n_flipped <- sum(sign(merged_gained$logFC.balanced) != sign(merged_gained[[paste0("logFC.", run_name)]]))

  cat(sprintf("  %-20s: %3d gained DEGs | %d sign-flipped vs balanced\n",
              run_name, nrow(merged_gained), n_flipped))

  leaked_list[[run_name]] <- merged_gained
  write.csv(merged_gained,
            file.path(output_dir, paste0("leaked_as_biology_", run_name, ".csv")),
            row.names = FALSE)
}
cat("\n")

# ------------------------------------------------------------------
# 5. Expression matrix comparison (PCA)
# ------------------------------------------------------------------
if (length(exprs_list) >= 2 && "balanced" %in% names(exprs_list)) {
  cat("=== 5. Expression Matrix PCA Comparison ===\n\n")

  phenodata <- read.csv("data/phenodata/samples.csv", stringsAsFactors = FALSE)

  pdf(file.path(output_dir, "pca_comparison.pdf"), width = 14, height = 10)

  for (run_name in names(exprs_list)) {
    mat <- as.matrix(exprs_list[[run_name]])
    mat <- mat[apply(mat, 1, var, na.rm = TRUE) > 0, , drop = FALSE]
    pca <- prcomp(t(mat), scale. = TRUE)
    pca_df <- data.frame(
      PC1 = pca$x[, 1], PC2 = pca$x[, 2],
      sample = rownames(pca$x)
    )
    pca_df <- merge(pca_df,
                    phenodata[, c("arraydatafile_exprscolumnnames",
                                  "secondaryaccession", "Gestational.Age.Category")],
                    by.x = "sample", by.y = "arraydatafile_exprscolumnnames",
                    all.x = TRUE)
    var_explained <- round(summary(pca)$importance[2, 1:2] * 100, 1)

    p <- ggplot(pca_df, aes(x = PC1, y = PC2,
                             colour = Gestational.Age.Category,
                             shape = secondaryaccession)) +
      geom_point(size = 3) +
      labs(title = paste("PCA:", run_name),
           x = paste0("PC1 (", var_explained[1], "%)"),
           y = paste0("PC2 (", var_explained[2], "%)"),
           colour = "Trimester", shape = "Dataset") +
      theme_minimal(base_size = 12)
    print(p)
  }
  dev.off()
  cat("  Saved pca_comparison.pdf\n\n")
}

# ------------------------------------------------------------------
# 6. Per-gene expression shift between runs
# ------------------------------------------------------------------
if ("balanced" %in% names(exprs_list) && length(exprs_list) > 1) {
  cat("=== 6. Per-Gene Mean Expression Shift ===\n\n")

  ref_means <- rowMeans(exprs_list[["balanced"]], na.rm = TRUE)

  for (run_name in setdiff(names(exprs_list), "balanced")) {
    cv_mat <- exprs_list[[run_name]]
    common <- intersect(names(ref_means), rownames(cv_mat))
    cv_means <- rowMeans(cv_mat[common, , drop = FALSE], na.rm = TRUE)
    r <- cor(ref_means[common], cv_means, use = "complete.obs")
    max_shift <- max(abs(ref_means[common] - cv_means), na.rm = TRUE)
    cat(sprintf("  %-20s: mean expr correlation=%.6f  max shift=%.4f\n",
                run_name, r, max_shift))
  }
  cat("\n")
}

# ------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------
cat("=== Summary ===\n\n")

summary_df <- merge(overlap_df, logfc_df, by = "run")
write.csv(summary_df, file.path(output_dir, "cv_summary.csv"), row.names = FALSE)

cat("Results saved to:", output_dir, "\n")
cat("Key files:\n")
cat("  cv_summary.csv            - overview metrics per run\n")
cat("  deg_overlap.csv           - DEG overlap with balanced reference\n")
cat("  logfc_correlation.csv     - logFC correlation metrics\n")
cat("  logfc_scatter.pdf         - logFC scatter plots\n")
cat("  pca_comparison.pdf        - PCA per run\n")
cat("  absorbed_as_batch_*.csv   - DEGs lost (biology -> batch)\n")
cat("  leaked_as_biology_*.csv   - DEGs gained (batch -> biology)\n")
