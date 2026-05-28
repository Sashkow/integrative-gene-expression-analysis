#!/usr/bin/env Rscript

library(ggplot2)
library(ggrepel)
library(sva)

input_dir <- "data/mapped/yehor/2026_04_27_9_preprocessed_datasets"
output_dir <- "data/mapped/yehor/8_yehor_preprocessed_datasets_enriched_sashko"
plot_dir <- "output/one_off/pca_3_yehor_gse37653/combat_3batch"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
exprs <- read.table(
  file.path(input_dir, "GSE37653_entrez_protein_coding.tsv"),
  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE
)
cat(sprintf("Expression matrix: %d genes x %d samples\n", nrow(exprs), ncol(exprs)))

pdata <- read.table(
  file.path(output_dir, "phenodata_placenta_1_2_enriched_sashko.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)
pdata <- pdata[pdata$dataset_id == "GSE37653", ]
pdata <- pdata[match(colnames(exprs), pdata$sample_id), ]

pdata$Batch3 <- ifelse(
  grepl("India", pdata$Batch), "India", pdata$Batch
)

cat("\nOriginal batch distribution:\n")
print(table(pdata$Batch))
cat("\nMerged 3-batch distribution:\n")
print(table(pdata$Batch3))
cat("\nSex distribution per merged batch:\n")
print(table(pdata$Batch3, pdata$fetux_sex_estimate))

batch <- as.factor(pdata$Batch3)
sex <- as.factor(pdata$fetux_sex_estimate)
mod <- model.matrix(~ sex)

cat("\nRunning ComBat (no ref.batch, 3 batches, sex covariate)...\n")
exprs_combat <- ComBat(dat = as.matrix(exprs), batch = batch, mod = mod)
cat(sprintf("ComBat output: %d genes x %d samples\n", nrow(exprs_combat), ncol(exprs_combat)))

make_pca_plot <- function(pca, pdata_row, var_exp, color_var, color_label, title_prefix, shape_var = NULL, shape_label = NULL) {
  pc <- data.frame(
    PC1 = pca$x[, 1], PC2 = pca$x[, 2],
    color = pdata_row[[color_var]],
    sample = pdata_row$sample_id
  )
  if (!is.null(shape_var)) pc$shape <- pdata_row[[shape_var]]

  aes_map <- if (!is.null(shape_var)) {
    aes(x = PC1, y = PC2, color = color, shape = shape, label = sample)
  } else {
    aes(x = PC1, y = PC2, color = color, label = sample)
  }

  p <- ggplot(pc, aes_map) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text_repel(size = 2.5, max.overlaps = 20, segment.size = 0.2) +
    labs(
      title = title_prefix,
      x = sprintf("PC1 (%.1f%%)", var_exp[1]),
      y = sprintf("PC2 (%.1f%%)", var_exp[2]),
      color = color_label
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      legend.text = element_text(size = 9)
    )

  if (!is.null(shape_label)) p <- p + labs(shape = shape_label)
  p
}

save_plot <- function(p, name) {
  svg(file.path(plot_dir, paste0(name, ".svg")), width = 10, height = 8)
  print(p)
  dev.off()
  png(file.path(plot_dir, paste0(name, ".png")), width = 10, height = 8, units = "in", res = 300)
  print(p)
  dev.off()
  cat("Saved:", name, "\n")
}

cat("\n=== PCA BEFORE ComBat ===\n")
pca_before <- prcomp(t(exprs), center = TRUE, scale. = FALSE)
var_before <- summary(pca_before)$importance[2, 1:5] * 100
for (i in 1:5) cat(sprintf("  PC%d: %.1f%%\n", i, var_before[i]))

cat("\n=== PCA AFTER ComBat ===\n")
pca_after <- prcomp(t(exprs_combat), center = TRUE, scale. = FALSE)
var_after <- summary(pca_after)$importance[2, 1:5] * 100
for (i in 1:5) cat(sprintf("  PC%d: %.1f%%\n", i, var_after[i]))

pdata$Batch <- as.character(pdata$Batch)
pdata$Batch3 <- as.character(pdata$Batch3)

save_plot(
  make_pca_plot(pca_before, pdata, var_before, "Batch", "Scan Batch", "Before ComBat (by batch)", "fetux_sex_estimate", "Sex"),
  "before_combat_batch"
)
save_plot(
  make_pca_plot(pca_before, pdata, var_before, "fetux_sex_estimate", "Sex", "Before ComBat (by sex)", "Batch", "Batch"),
  "before_combat_sex"
)
save_plot(
  make_pca_plot(pca_after, pdata, var_after, "Batch", "Scan Batch", "After ComBat (by batch)", "fetux_sex_estimate", "Sex"),
  "after_combat_batch"
)
save_plot(
  make_pca_plot(pca_after, pdata, var_after, "fetux_sex_estimate", "Sex", "After ComBat (by sex)", "Batch", "Batch"),
  "after_combat_sex"
)

write.table(exprs_combat, file.path(plot_dir, "GSE37653_combat_3batch.tsv"),
            sep = "\t", quote = FALSE)
cat("\nCorrected matrix saved to:", file.path(plot_dir, "GSE37653_combat_3batch.tsv"), "\n")
cat("Plots saved to:", plot_dir, "\n")
