#!/usr/bin/env Rscript

library(ggplot2)
library(ggrepel)

yehor_dir <- "data/mapped/yehor/2026_04_27_9_preprocessed_datasets"

datasets <- c("GSE9984", "GSE22490", "GSE122214", "GSE37653")
cat("Loading datasets:", paste(datasets, collapse = ", "), "\n")

exprs_list <- list()
for (ds in datasets) {
  f <- file.path(yehor_dir, paste0(ds, "_entrez_protein_coding.tsv"))
  exprs_list[[ds]] <- read.table(f, header = TRUE, row.names = 1,
                                  sep = "\t", check.names = FALSE)
  cat(sprintf("  %s: %d genes x %d samples\n", ds,
              nrow(exprs_list[[ds]]), ncol(exprs_list[[ds]])))
}

pdata <- read.table(file.path(yehor_dir, "phenodata_placenta_1_2.tsv"),
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)

common_genes <- Reduce(intersect, lapply(exprs_list, rownames))
cat(sprintf("\nCommon genes across all datasets: %d\n", length(common_genes)))

merged <- do.call(cbind, unname(lapply(exprs_list, function(e) e[common_genes, ])))
cat(sprintf("Merged matrix: %d genes x %d samples\n", nrow(merged), ncol(merged)))

sample_ids <- colnames(merged)
pdata_matched <- pdata[match(sample_ids, pdata$sample_id), ]

matched <- sum(!is.na(pdata_matched$sample_id))
cat(sprintf("Phenodata matched: %d / %d samples\n", matched, length(sample_ids)))

if (matched < length(sample_ids)) {
  missing <- sample_ids[is.na(pdata_matched$sample_id)]
  cat("Unmatched samples:", paste(missing, collapse = ", "), "\n")
}

pca <- prcomp(t(merged), center = TRUE, scale. = FALSE)
var_exp <- summary(pca)$importance[2, 1:min(10, ncol(pca$x))] * 100
cat("\nVariance explained:\n")
for (i in seq_along(var_exp)) {
  cat(sprintf("  PC%d: %.1f%%\n", i, var_exp[i]))
}

output_dir <- "output/one_off/pca_3_yehor_gse37653"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

pc_data <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  dataset = pdata_matched$dataset_id,
  trimester = pdata_matched$trimester,
  sample = pdata_matched$sample_id
)

p1 <- ggplot(pc_data, aes(x = PC1, y = PC2, color = dataset)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(type = "norm", level = 0.95, linetype = 2) +
  labs(
    title = "PCA: 3 Yehor datasets + GSE37653 (colored by dataset)",
    x = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y = sprintf("PC2 (%.1f%%)", var_exp[2]),
    color = "Dataset"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))

p2 <- ggplot(pc_data, aes(x = PC1, y = PC2, color = trimester)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(type = "norm", level = 0.95, linetype = 2) +
  labs(
    title = "PCA: 3 Yehor datasets + GSE37653 (colored by trimester)",
    x = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y = sprintf("PC2 (%.1f%%)", var_exp[2]),
    color = "Trimester"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))

p3 <- ggplot(pc_data, aes(x = PC1, y = PC2, color = dataset, label = sample)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_text_repel(size = 2, max.overlaps = 20, segment.size = 0.2) +
  labs(
    title = "PCA: 3 Yehor datasets + GSE37653 (with sample labels)",
    x = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y = sprintf("PC2 (%.1f%%)", var_exp[2]),
    color = "Dataset"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))

for (plot_info in list(
  list(p = p1, name = "pca_by_dataset"),
  list(p = p2, name = "pca_by_trimester"),
  list(p = p3, name = "pca_sample_labels")
)) {
  svg_file <- file.path(output_dir, paste0(plot_info$name, ".svg"))
  svg(svg_file, width = 10, height = 8)
  print(plot_info$p)
  dev.off()

  png_file <- file.path(output_dir, paste0(plot_info$name, ".png"))
  png(png_file, width = 10, height = 8, units = "in", res = 300)
  print(plot_info$p)
  dev.off()

  cat("Saved:", svg_file, "\n")
}

cat("\n=== Separate PCA for GSE37653 with batch metadata ===\n")

geo_meta <- data.frame(
  sample_id = paste0("GSM9250", 71:95),
  title = c(
    "Placenta_villi_weeks6_1","Placenta_villi_weeks6_2","Placenta_villi_weeks6_3",
    "Placenta_villi_weeks6_4","Placenta_villi_weeks6_5","Placenta_villi_weeks6_6",
    "Placenta_villi_weeks6_7","Placenta_villi_weeks6_8",
    "Placenta_villi_weeks7_1","Placenta_villi_weeks7_2","Placenta_villi_weeks7_3",
    "Placenta_villi_weeks7_4","Placenta_villi_weeks7_5","Placenta_villi_weeks7_6",
    "Placenta_villi_weeks7_7","Placenta_villi_weeks7_8","Placenta_villi_weeks7_9",
    "Placenta_villi_weeks8_1","Placenta_villi_weeks8_2","Placenta_villi_weeks8_3",
    "Placenta_villi_weeks8_4","Placenta_villi_weeks8_5","Placenta_villi_weeks8_6",
    "Placenta_villi_weeks8_7","Placenta_villi_weeks8_8"
  ),
  gest_week = c(rep(6,8), rep(7,9), rep(8,8)),
  maternal_age = c(42,26,25,27,26,30,21,26, 26,24,26,25,26,25,25,24,22, 30,26,28,23,28,26,20,20),
  slide = c(
    rep("35825101_S03", 6), "123_S01","123_S01",
    rep("35825101_S03", 6), "123_S01","123_S04","123_S04",
    rep("35825101_S03", 6), "123_S01","123_S04"
  ),
  array_design = c(
    rep("070925", 6), "US80303135","US80303135",
    rep("070925", 6), "US80303135","US80303135","US80303135",
    rep("070925", 6), "US80303135","US80303135"
  ),
  stringsAsFactors = FALSE
)
exprs_37653 <- exprs_list[["GSE37653"]]
pdata_37653 <- pdata[pdata$dataset_id == "GSE37653", ]
pdata_37653 <- pdata_37653[match(colnames(exprs_37653), pdata_37653$sample_id), ]

geo_37653 <- geo_meta[match(colnames(exprs_37653), geo_meta$sample_id), ]

pca_37653 <- prcomp(t(exprs_37653), center = TRUE, scale. = FALSE)
var_37653 <- summary(pca_37653)$importance[2, 1:min(10, ncol(pca_37653$x))] * 100
cat("GSE37653 variance explained:\n")
for (i in seq_along(var_37653)) cat(sprintf("  PC%d: %.1f%%\n", i, var_37653[i]))

pc_37653 <- data.frame(
  PC1 = pca_37653$x[, 1],
  PC2 = pca_37653$x[, 2],
  trimester = pdata_37653$trimester,
  sex = pdata_37653$fetux_sex_estimate,
  sample = pdata_37653$sample_id,
  slide = geo_37653$slide,
  array_design = geo_37653$array_design,
  gest_week = as.factor(geo_37653$gest_week),
  maternal_age = geo_37653$maternal_age
)

p4 <- ggplot(pc_37653, aes(x = PC1, y = PC2, color = trimester, label = sample)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.size = 0.2) +
  labs(
    title = "PCA: GSE37653 only (colored by trimester)",
    x = sprintf("PC1 (%.1f%%)", var_37653[1]),
    y = sprintf("PC2 (%.1f%%)", var_37653[2]),
    color = "Trimester"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))

p5 <- ggplot(pc_37653, aes(x = PC1, y = PC2, color = sex, label = sample)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.size = 0.2) +
  labs(
    title = "PCA: GSE37653 only (colored by sex estimate)",
    x = sprintf("PC1 (%.1f%%)", var_37653[1]),
    y = sprintf("PC2 (%.1f%%)", var_37653[2]),
    color = "Sex"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))

p6 <- ggplot(pc_37653, aes(x = PC1, y = PC2, color = slide, label = sample)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.size = 0.2) +
  labs(
    title = "PCA: GSE37653 (colored by slide/subarray)",
    x = sprintf("PC1 (%.1f%%)", var_37653[1]),
    y = sprintf("PC2 (%.1f%%)", var_37653[2]),
    color = "Slide"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))

p7 <- ggplot(pc_37653, aes(x = PC1, y = PC2, color = array_design, label = sample)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.size = 0.2) +
  labs(
    title = "PCA: GSE37653 (colored by array design)",
    x = sprintf("PC1 (%.1f%%)", var_37653[1]),
    y = sprintf("PC2 (%.1f%%)", var_37653[2]),
    color = "Array Design"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))

p8 <- ggplot(pc_37653, aes(x = PC1, y = PC2, color = gest_week, label = sample)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.size = 0.2) +
  labs(
    title = "PCA: GSE37653 (colored by gestational week)",
    x = sprintf("PC1 (%.1f%%)", var_37653[1]),
    y = sprintf("PC2 (%.1f%%)", var_37653[2]),
    color = "Gest. Week"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))

p9 <- ggplot(pc_37653, aes(x = PC1, y = PC2, color = maternal_age, label = sample)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.size = 0.2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    title = "PCA: GSE37653 (colored by maternal age)",
    x = sprintf("PC1 (%.1f%%)", var_37653[1]),
    y = sprintf("PC2 (%.1f%%)", var_37653[2]),
    color = "Maternal Age"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))

p10 <- ggplot(pc_37653, aes(x = PC1, y = PC2, color = sex, shape = slide, label = sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.size = 0.2) +
  scale_shape_manual(values = c("123_S01" = 15, "123_S04" = 17, "35825101_S03" = 16)) +
  labs(
    title = "PCA: GSE37653 (sex + slide batch)",
    x = sprintf("PC1 (%.1f%%)", var_37653[1]),
    y = sprintf("PC2 (%.1f%%)", var_37653[2]),
    color = "Sex",
    shape = "Slide"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))

for (plot_info in list(
  list(p = p4, name = "pca_gse37653_trimester"),
  list(p = p5, name = "pca_gse37653_sex"),
  list(p = p6, name = "pca_gse37653_slide"),
  list(p = p7, name = "pca_gse37653_array_design"),
  list(p = p8, name = "pca_gse37653_gest_week"),
  list(p = p9, name = "pca_gse37653_maternal_age"),
  list(p = p10, name = "pca_gse37653_sex_slide")
)) {
  svg_file <- file.path(output_dir, paste0(plot_info$name, ".svg"))
  svg(svg_file, width = 10, height = 8)
  print(plot_info$p)
  dev.off()

  png_file <- file.path(output_dir, paste0(plot_info$name, ".png"))
  png(png_file, width = 10, height = 8, units = "in", res = 300)
  print(plot_info$p)
  dev.off()

  cat("Saved:", svg_file, "\n")
}

scan_date <- c(
  GSM925071="2009-12-05", GSM925072="2009-12-05", GSM925073="2009-12-05", GSM925074="2009-12-05",
  GSM925075="2010-02-22", GSM925076="2010-02-22",
  GSM925077="2010-07-16", GSM925078="2010-07-16",
  GSM925079="2009-12-05", GSM925080="2009-12-05", GSM925081="2009-12-05", GSM925082="2009-12-05",
  GSM925083="2010-02-22", GSM925084="2010-02-22",
  GSM925085="2010-07-16",
  GSM925086="no date", GSM925087="no date",
  GSM925088="2009-12-05", GSM925089="2009-12-05", GSM925090="2009-12-05", GSM925091="2009-12-05",
  GSM925092="2010-02-22", GSM925093="2010-02-22",
  GSM925094="2010-07-16",
  GSM925095="no date"
)

scan_env <- c(
  GSM925071="Dec 2009 Singapore v2.5", GSM925072="Dec 2009 Singapore v2.5",
  GSM925073="Dec 2009 Singapore v2.5", GSM925074="Dec 2009 Singapore v2.5",
  GSM925075="Feb 2010 India v2.5",     GSM925076="Feb 2010 India v2.5",
  GSM925077="Jul 2010 India v2.6",     GSM925078="Jul 2010 India v2.6",
  GSM925079="Dec 2009 Singapore v2.5", GSM925080="Dec 2009 Singapore v2.5",
  GSM925081="Dec 2009 Singapore v2.5", GSM925082="Dec 2009 Singapore v2.5",
  GSM925083="Feb 2010 India v2.5",     GSM925084="Feb 2010 India v2.5",
  GSM925085="Jul 2010 India v2.6",
  GSM925086="No header",               GSM925087="No header",
  GSM925088="Dec 2009 Singapore v2.5", GSM925089="Dec 2009 Singapore v2.5",
  GSM925090="Dec 2009 Singapore v2.5", GSM925091="Dec 2009 Singapore v2.5",
  GSM925092="Feb 2010 India v2.5",     GSM925093="Feb 2010 India v2.5",
  GSM925094="Jul 2010 India v2.6",
  GSM925095="No header"
)

pc_37653$scan_date <- scan_date[pc_37653$sample]
pc_37653$scan_env <- scan_env[pc_37653$sample]

p_scan <- ggplot(pc_37653, aes(x = PC1, y = PC2, color = scan_env, shape = sex, label = sample)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.size = 0.2) +
  labs(
    title = "PCA: GSE37653 (colored by scan environment)",
    x = sprintf("PC1 (%.1f%%)", var_37653[1]),
    y = sprintf("PC2 (%.1f%%)", var_37653[2]),
    color = "Scan Environment", shape = "Sex"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    legend.text = element_text(size = 9)
  )

for (plot_info in list(list(p = p_scan, name = "pca_gse37653_scan_env"))) {
  svg_file <- file.path(output_dir, paste0(plot_info$name, ".svg"))
  svg(svg_file, width = 10, height = 8)
  print(plot_info$p)
  dev.off()
  png_file <- file.path(output_dir, paste0(plot_info$name, ".png"))
  png(png_file, width = 10, height = 8, units = "in", res = 300)
  print(plot_info$p)
  dev.off()
  cat("Saved:", svg_file, "\n")
}

cat("\n=== PCA after removing sex variation ===\n")
sex_factor <- as.factor(pdata_37653$fetux_sex_estimate)
design_sex <- model.matrix(~ sex_factor)
fit_sex <- limma::lmFit(exprs_37653, design_sex)
exprs_no_sex <- exprs_37653 - fit_sex$coefficients[, 2] %*% t(design_sex[, 2])

pca_nosex <- prcomp(t(exprs_no_sex), center = TRUE, scale. = FALSE)
var_nosex <- summary(pca_nosex)$importance[2, 1:min(10, ncol(pca_nosex$x))] * 100
cat("Variance explained (sex removed):\n")
for (i in seq_along(var_nosex)) cat(sprintf("  PC%d: %.1f%%\n", i, var_nosex[i]))

pc_nosex <- data.frame(
  PC1 = pca_nosex$x[, 1],
  PC2 = pca_nosex$x[, 2],
  sex = pdata_37653$fetux_sex_estimate,
  slide = geo_37653$slide,
  array_design = geo_37653$array_design,
  gest_week = as.factor(geo_37653$gest_week),
  maternal_age = geo_37653$maternal_age,
  sample = pdata_37653$sample_id
)

p_ns1 <- ggplot(pc_nosex, aes(x = PC1, y = PC2, color = slide, shape = sex, label = sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.size = 0.2) +
  labs(
    title = "PCA after removing sex (colored by slide)",
    x = sprintf("PC1 (%.1f%%)", var_nosex[1]),
    y = sprintf("PC2 (%.1f%%)", var_nosex[2]),
    color = "Slide", shape = "Sex"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))

p_ns2 <- ggplot(pc_nosex, aes(x = PC1, y = PC2, color = gest_week, shape = slide, label = sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.size = 0.2) +
  scale_shape_manual(values = c("123_S01" = 15, "123_S04" = 17, "35825101_S03" = 16)) +
  labs(
    title = "PCA after removing sex (colored by gest. week)",
    x = sprintf("PC1 (%.1f%%)", var_nosex[1]),
    y = sprintf("PC2 (%.1f%%)", var_nosex[2]),
    color = "Gest. Week", shape = "Slide"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))

p_ns3 <- ggplot(pc_nosex, aes(x = PC1, y = PC2, color = maternal_age, shape = slide, label = sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.size = 0.2) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_shape_manual(values = c("123_S01" = 15, "123_S04" = 17, "35825101_S03" = 16)) +
  labs(
    title = "PCA after removing sex (colored by maternal age)",
    x = sprintf("PC1 (%.1f%%)", var_nosex[1]),
    y = sprintf("PC2 (%.1f%%)", var_nosex[2]),
    color = "Maternal Age", shape = "Slide"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))

for (plot_info in list(
  list(p = p_ns1, name = "pca_gse37653_nosex_slide"),
  list(p = p_ns2, name = "pca_gse37653_nosex_gestweek"),
  list(p = p_ns3, name = "pca_gse37653_nosex_maternal_age")
)) {
  svg_file <- file.path(output_dir, paste0(plot_info$name, ".svg"))
  svg(svg_file, width = 10, height = 8)
  print(plot_info$p)
  dev.off()

  png_file <- file.path(output_dir, paste0(plot_info$name, ".png"))
  png(png_file, width = 10, height = 8, units = "in", res = 300)
  print(plot_info$p)
  dev.off()

  cat("Saved:", svg_file, "\n")
}

pc_nosex$scan_date <- scan_date[pc_nosex$sample]
pc_nosex$scan_env <- scan_env[pc_nosex$sample]

p_ns_env <- ggplot(pc_nosex, aes(x = PC1, y = PC2, color = scan_env, shape = sex, label = sample)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.size = 0.2) +
  labs(
    title = "PCA after removing sex (colored by scan environment)",
    x = sprintf("PC1 (%.1f%%)", var_nosex[1]),
    y = sprintf("PC2 (%.1f%%)", var_nosex[2]),
    color = "Scan Environment", shape = "Sex"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    legend.text = element_text(size = 9)
  )

for (plot_info in list(list(p = p_ns_env, name = "pca_gse37653_nosex_scan_env"))) {
  svg_file <- file.path(output_dir, paste0(plot_info$name, ".svg"))
  svg(svg_file, width = 10, height = 8)
  print(plot_info$p)
  dev.off()
  png_file <- file.path(output_dir, paste0(plot_info$name, ".png"))
  png(png_file, width = 10, height = 8, units = "in", res = 300)
  print(plot_info$p)
  dev.off()
  cat("Saved:", svg_file, "\n")
}

cat("\nAll plots saved to:", output_dir, "\n")
