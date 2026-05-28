#!/usr/bin/env Rscript
#
# Sex-stratified DE analysis using existing softimpute + combat_ref matrices.
# Produces 4 DE comparisons and PCA plot for article Results paragraph.
#
# Usage:
#   Rscript one_off_scripts/sex_stratified_de_analysis.R 6ds
#   Rscript one_off_scripts/sex_stratified_de_analysis.R 8ds
#   Rscript one_off_scripts/sex_stratified_de_analysis.R 9ds
#   Rscript one_off_scripts/sex_stratified_de_analysis.R 9ds known_only

suppressPackageStartupMessages({
  library(limma)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || !args[1] %in% c("6ds", "8ds", "9ds", "prater")) {
  stop("Usage: Rscript sex_stratified_de_analysis.R [6ds|8ds|9ds|prater] [known_only]")
}
pack <- args[1]
known_only <- length(args) >= 2 && args[2] == "known_only"

if (pack == "6ds") {
  exprs_file <- "output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko/exprs_softimpute_combat_ref.tsv"
  pheno_file <- "data/mapped/yehor/8_yehor_preprocessed_datasets_enriched_sashko/phenodata_placenta_1_2_enriched_sashko.tsv"
  output_dir <- "output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko/sex_stratified"
} else if (pack == "8ds") {
  exprs_file <- "output/yehor_sashko/phase2b_1_2_yehor_8ds_enriched_sashko/exprs_softimpute_combat_ref.tsv"
  pheno_file <- "data/mapped/yehor/8_yehor_preprocessed_datasets_enriched_sashko/phenodata_placenta_1_2_enriched_sashko.tsv"
  output_dir <- "output/yehor_sashko/phase2b_1_2_yehor_8ds_enriched_sashko/sex_stratified"
} else if (pack == "prater") {
  exprs_file <- "output/yehor_sashko/phase2b_1_2_yehor_6ds_trimmed_to_prater/exprs_softimpute_combat_ref.tsv"
  pheno_file <- "data/mapped/yehor/8_yehor_preprocessed_datasets_enriched_sashko/phenodata_placenta_1_2_trimmed_to_prater.tsv"
  output_dir <- "output/yehor_sashko/phase2b_1_2_yehor_6ds_trimmed_to_prater/sex_stratified"
} else {
  exprs_file <- "output/yehor_sashko/phase2b_1_2_yehor_9ds_enriched_sashko/exprs_softimpute_combat_ref.tsv"
  pheno_file <- "data/mapped/yehor/9_yehor_preprocessed_datasets_enriched_sashko/phenodata_placenta_1_2_enriched_sashko.tsv"
  output_dir <- "output/yehor_sashko/phase2b_1_2_yehor_9ds_enriched_sashko/sex_stratified"
}
if (known_only) {
  output_dir <- paste0(output_dir, "_known_only")
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

run_label <- if (known_only) paste0(pack, " known_only") else pack
cat(sprintf("=== Sex-stratified DE analysis (%s) ===\n\n", run_label))

# --- Load data ---
exprs <- read.delim(exprs_file, row.names = 1, check.names = FALSE)
pheno <- read.delim(pheno_file, stringsAsFactors = FALSE)
samples_csv <- read.csv("data/phenodata/samples.csv", stringsAsFactors = FALSE)

# Filter phenodata to samples present in normalized matrix
pheno <- pheno[pheno$sample_id %in% colnames(exprs), ]
pheno <- pheno[pheno$condition == "healthy" &
               pheno$trimester %in% c("First trimester", "Second trimester"), ]

# Determine known vs predicted sex for all samples
samples_csv$gsm <- gsub("\\.(CEL|cel)$", "", samples_csv$arraydatafile_exprscolumnnames)
samples_lookup <- setNames(samples_csv$Fetus.Sex, samples_csv$gsm)
pheno$sex_source <- ifelse(
  pheno$sample_id %in% names(samples_lookup) &
    !is.na(samples_lookup[pheno$sample_id]) &
    nzchar(trimws(samples_lookup[pheno$sample_id])) &
    !tolower(trimws(samples_lookup[pheno$sample_id])) %in% c("na", "n/a", "unknown", "_", "none"),
  "known", "predicted"
)

if (known_only) {
  n_before <- nrow(pheno)
  pheno <- pheno[pheno$sex_source == "known", ]
  cat(sprintf("known_only filter: %d → %d samples (dropped %d with predicted sex)\n\n",
              n_before, nrow(pheno), n_before - nrow(pheno)))
}

exprs <- exprs[, pheno$sample_id]

cat(sprintf("Loaded: %d genes x %d samples\n", nrow(exprs), ncol(exprs)))
cat(sprintf("  First trimester: %d\n", sum(pheno$trimester == "First trimester")))
cat(sprintf("  Second trimester: %d\n\n", sum(pheno$trimester == "Second trimester")))

cat("=== Fetal sex: known vs predicted ===\n\n")

sex_counts <- as.data.frame.matrix(
  table(paste(pheno$trimester, pheno$fetux_sex_estimate, sep = " | "),
        pheno$sex_source)
)
print(sex_counts)
cat("\n")

for (tri in c("First trimester", "Second trimester")) {
  tri_short <- ifelse(tri == "First trimester", "1T", "2T")
  sub <- pheno[pheno$trimester == tri, ]
  for (sx in c("f", "m")) {
    n_known <- sum(sub$fetux_sex_estimate == sx & sub$sex_source == "known")
    n_pred  <- sum(sub$fetux_sex_estimate == sx & sub$sex_source == "predicted")
    sx_label <- ifelse(sx == "f", "female", "male")
    cat(sprintf("  %s %s: %d known + %d predicted = %d total\n",
                tri_short, sx_label, n_known, n_pred, n_known + n_pred))
  }
}
cat("\n")

# --- DE helper ---
fdr_thresh <- 0.05
logfc_thresh <- 1.0

run_de <- function(exprs_mat, pdata, group_col, baseline, contrast, label) {
  keep <- pdata[[group_col]] %in% c(baseline, contrast)
  de_exprs <- exprs_mat[, pdata$sample_id[keep]]
  de_group <- factor(pdata[[group_col]][keep], levels = c(baseline, contrast))

  cat(sprintf("\n--- %s ---\n", label))
  cat(sprintf("  Baseline (%s): %d samples\n", baseline, sum(de_group == baseline)))
  cat(sprintf("  Contrast (%s): %d samples\n", contrast, sum(de_group == contrast)))

  design <- model.matrix(~ de_group)
  fit <- lmFit(de_exprs, design)
  fit <- eBayes(fit)

  results <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  results$gene <- rownames(results)
  results <- results[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

  sig <- results[!is.na(results$adj.P.Val) & !is.na(results$logFC) &
                   results$adj.P.Val < fdr_thresh & abs(results$logFC) > logfc_thresh, ]

  cat(sprintf("  Genes tested: %d\n", nrow(results)))
  cat(sprintf("  Significant (FDR < %.2f, |logFC| > %.1f): %d\n",
              fdr_thresh, logfc_thresh, nrow(sig)))
  cat(sprintf("    Up-regulated: %d\n", sum(sig$logFC > 0)))
  cat(sprintf("    Down-regulated: %d\n", sum(sig$logFC < 0)))

  list(results = results, significant = sig, label = label,
       n_baseline = sum(de_group == baseline),
       n_contrast = sum(de_group == contrast))
}

# --- Run 4 DE comparisons ---
de1 <- run_de(exprs, pheno[pheno$trimester == "First trimester", ],
              "fetux_sex_estimate", "f", "m",
              "1T: male vs female")

de2 <- run_de(exprs, pheno[pheno$trimester == "Second trimester", ],
              "fetux_sex_estimate", "f", "m",
              "2T: male vs female")

de3 <- run_de(exprs, pheno[pheno$fetux_sex_estimate == "m", ],
              "trimester", "First trimester", "Second trimester",
              "Males: 1T vs 2T")

de4 <- run_de(exprs, pheno[pheno$fetux_sex_estimate == "f", ],
              "trimester", "First trimester", "Second trimester",
              "Females: 1T vs 2T")

# --- Save DE tables ---
de_list <- list(
  "1t_male_vs_female" = de1,
  "2t_male_vs_female" = de2,
  "males_1t_vs_2t"    = de3,
  "females_1t_vs_2t"  = de4
)

for (key in names(de_list)) {
  de <- de_list[[key]]
  write.table(de$results,
              file.path(output_dir, paste0("difexp_", key, ".tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(de$significant,
              file.path(output_dir, paste0("difexp_significant_", key, ".tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# --- PCA plot ---
cat("\n\n=== PCA ===\n")
pc <- prcomp(t(exprs), center = TRUE, scale. = FALSE)
var_pct <- 100 * summary(pc)$importance[2, ]

group_labels <- paste0(
  ifelse(pheno$trimester == "First trimester", "1T", "2T"), "-",
  ifelse(pheno$fetux_sex_estimate == "f", "F", "M")
)
group_factor <- factor(group_labels, levels = c("1T-F", "1T-M", "2T-F", "2T-M"))

group_colors <- c("1T-F" = "#E78AC3", "1T-M" = "#66C2A5",
                   "2T-F" = "#FC8D62", "2T-M" = "#8DA0CB")
group_shapes <- c("1T-F" = 16, "1T-M" = 16, "2T-F" = 17, "2T-M" = 17)

png(file.path(output_dir, "pca_trimester_sex.png"),
    width = 1000, height = 800, res = 120)
par(mar = c(5, 5, 3, 2))
plot(pc$x[, 1], pc$x[, 2],
     col = group_colors[group_labels],
     pch = group_shapes[group_labels],
     cex = 1.5,
     xlab = sprintf("PC1 (%.1f%%)", var_pct[1]),
     ylab = sprintf("PC2 (%.1f%%)", var_pct[2]),
     main = sprintf("PCA — %s softimpute+combat_ref (%d samples)", run_label, ncol(exprs)))
legend("topright",
       legend = sprintf("%s (n=%d)", levels(group_factor), table(group_factor)),
       col = group_colors[levels(group_factor)],
       pch = group_shapes[levels(group_factor)],
       pt.cex = 1.5, cex = 0.9)
dev.off()
cat(sprintf("PCA plot saved: %s/pca_trimester_sex.png\n", output_dir))

# --- Summary ---
summary_lines <- c(
  sprintf("Sex-stratified DE analysis — %s pack", run_label),
  sprintf("Date: %s", Sys.Date()),
  sprintf("Total samples: %d", ncol(exprs)),
  sprintf("Genes in normalized matrix: %d", nrow(exprs)),
  "",
  "=== Fetal sex counts ===",
  sprintf("First trimester: %d female (%d known + %d predicted), %d male (%d known + %d predicted)",
          sum(pheno$trimester == "First trimester" & pheno$fetux_sex_estimate == "f"),
          sum(pheno$trimester == "First trimester" & pheno$fetux_sex_estimate == "f" & pheno$sex_source == "known"),
          sum(pheno$trimester == "First trimester" & pheno$fetux_sex_estimate == "f" & pheno$sex_source == "predicted"),
          sum(pheno$trimester == "First trimester" & pheno$fetux_sex_estimate == "m"),
          sum(pheno$trimester == "First trimester" & pheno$fetux_sex_estimate == "m" & pheno$sex_source == "known"),
          sum(pheno$trimester == "First trimester" & pheno$fetux_sex_estimate == "m" & pheno$sex_source == "predicted")),
  sprintf("Second trimester: %d female (%d known + %d predicted), %d male (%d known + %d predicted)",
          sum(pheno$trimester == "Second trimester" & pheno$fetux_sex_estimate == "f"),
          sum(pheno$trimester == "Second trimester" & pheno$fetux_sex_estimate == "f" & pheno$sex_source == "known"),
          sum(pheno$trimester == "Second trimester" & pheno$fetux_sex_estimate == "f" & pheno$sex_source == "predicted"),
          sum(pheno$trimester == "Second trimester" & pheno$fetux_sex_estimate == "m"),
          sum(pheno$trimester == "Second trimester" & pheno$fetux_sex_estimate == "m" & pheno$sex_source == "known"),
          sum(pheno$trimester == "Second trimester" & pheno$fetux_sex_estimate == "m" & pheno$sex_source == "predicted")),
  "",
  sprintf("=== DE results (FDR < %.2f, |logFC| > %.1f) ===", fdr_thresh, logfc_thresh)
)

for (key in names(de_list)) {
  de <- de_list[[key]]
  summary_lines <- c(summary_lines,
    sprintf("%s: %d baseline vs %d contrast → %d significant (%d up, %d down)",
            de$label, de$n_baseline, de$n_contrast,
            nrow(de$significant),
            sum(de$significant$logFC > 0),
            sum(de$significant$logFC < 0)))
}

summary_lines <- c(summary_lines, "",
  "=== Article paragraph (fill in [N]) ===",
  sprintf("Imputation of missing values allowed us to integrate data from all %d available samples.", ncol(exprs)),
  sprintf(paste0(
    "In addition, fetal sex was predicted for samples with missing metadata: ",
    "%d male and %d female fetuses were identified among first-trimester samples, ",
    "in addition to the previously known %d male and %d female samples; ",
    "and %d male and %d female fetuses were identified among second-trimester samples, ",
    "in addition to the previously known %d male and %d female samples."),
    sum(pheno$trimester == "First trimester" & pheno$fetux_sex_estimate == "m" & pheno$sex_source == "predicted"),
    sum(pheno$trimester == "First trimester" & pheno$fetux_sex_estimate == "f" & pheno$sex_source == "predicted"),
    sum(pheno$trimester == "First trimester" & pheno$fetux_sex_estimate == "m" & pheno$sex_source == "known"),
    sum(pheno$trimester == "First trimester" & pheno$fetux_sex_estimate == "f" & pheno$sex_source == "known"),
    sum(pheno$trimester == "Second trimester" & pheno$fetux_sex_estimate == "m" & pheno$sex_source == "predicted"),
    sum(pheno$trimester == "Second trimester" & pheno$fetux_sex_estimate == "f" & pheno$sex_source == "predicted"),
    sum(pheno$trimester == "Second trimester" & pheno$fetux_sex_estimate == "m" & pheno$sex_source == "known"),
    sum(pheno$trimester == "Second trimester" & pheno$fetux_sex_estimate == "f" & pheno$sex_source == "known")),
  sprintf(paste0(
    "Differential expression analysis identified %d genes between first-trimester male ",
    "and first-trimester female samples, %d genes between second-trimester male and ",
    "second-trimester female samples, %d genes between first- and second-trimester male samples, ",
    "and %d genes between first- and second-trimester female samples."),
    nrow(de1$significant), nrow(de2$significant),
    nrow(de3$significant), nrow(de4$significant)),
  "PCA plots showed that, after preprocessing, samples clustered primarily according to gestational age and fetal sex."
)

writeLines(summary_lines, file.path(output_dir, "summary.txt"))
cat(sprintf("\nSummary saved: %s/summary.txt\n", output_dir))

cat("\n=== ARTICLE PARAGRAPH ===\n\n")
for (line in summary_lines[grep("^Imputation|^In addition|^Differential|^PCA", summary_lines)]) {
  cat(line, "\n")
}
cat("\nDone.\n")
