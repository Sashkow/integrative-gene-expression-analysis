options(bitmapType = "cairo")

base_dir <- "output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko"
output_dir <- "articles/imputation_article/data/logfc_stability"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

si <- read.delim(file.path(base_dir, "difexp_softimpute_combat_ref.tsv"), stringsAsFactors = FALSE)
io <- read.delim(file.path(base_dir, "difexp_none_combat_ref.tsv"), stringsAsFactors = FALSE)

si$gene <- as.character(si$gene)
io$gene <- as.character(io$gene)

shared <- intersect(si$gene, io$gene)
idx_si <- match(shared, si$gene)
idx_io <- match(shared, io$gene)

comp <- data.frame(
  gene = shared,
  logFC_softimpute = si$logFC[idx_si],
  logFC_intersect  = io$logFC[idx_io],
  stringsAsFactors = FALSE
)
comp$diff <- comp$logFC_softimpute - comp$logFC_intersect
comp$abs_diff <- abs(comp$diff)

r_val <- cor(comp$logFC_softimpute, comp$logFC_intersect)

cat(sprintf("Shared genes: %d\n", nrow(comp)))
cat(sprintf("Pearson r: %.4f\n", r_val))
cat(sprintf("Mean diff: %+.4f\n", mean(comp$diff)))
cat(sprintf("Median abs diff: %.4f\n", median(comp$abs_diff)))
cat(sprintf("P95 abs diff: %.4f\n", quantile(comp$abs_diff, 0.95)))
cat(sprintf("P99 abs diff: %.4f\n", quantile(comp$abs_diff, 0.99)))
cat(sprintf("Max abs diff: %.4f\n", max(comp$abs_diff)))

write.csv(comp, file.path(output_dir, "logfc_stability.csv"), row.names = FALSE)

summary_txt <- sprintf(
"logFC stability: softImpute vs intersection-only (shared %d genes)

Pearson r:        %.4f
Mean difference:  %+.4f
Median |diff|:    %.4f
SD |diff|:        %.4f
P75 |diff|:       %.4f
P90 |diff|:       %.4f
P95 |diff|:       %.4f
P99 |diff|:       %.4f
Max |diff|:       %.4f

Threshold counts:
  |diff| > 0.01: %d (%.1f%%)
  |diff| > 0.05: %d (%.1f%%)
  |diff| > 0.10: %d (%.1f%%)
  |diff| > 0.20: %d (%.1f%%)
  |diff| > 0.50: %d (%.1f%%)
",
nrow(comp), r_val,
mean(comp$diff), median(comp$abs_diff), sd(comp$abs_diff),
quantile(comp$abs_diff, 0.75), quantile(comp$abs_diff, 0.90),
quantile(comp$abs_diff, 0.95), quantile(comp$abs_diff, 0.99),
max(comp$abs_diff),
sum(comp$abs_diff > 0.01), 100 * mean(comp$abs_diff > 0.01),
sum(comp$abs_diff > 0.05), 100 * mean(comp$abs_diff > 0.05),
sum(comp$abs_diff > 0.10), 100 * mean(comp$abs_diff > 0.10),
sum(comp$abs_diff > 0.20), 100 * mean(comp$abs_diff > 0.20),
sum(comp$abs_diff > 0.50), 100 * mean(comp$abs_diff > 0.50))

writeLines(summary_txt, file.path(output_dir, "logfc_stability.txt"))
cat("\nSaved to:", output_dir, "\n")
