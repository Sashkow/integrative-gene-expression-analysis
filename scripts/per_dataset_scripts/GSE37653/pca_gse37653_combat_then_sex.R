#!/usr/bin/env Rscript

library(ggplot2)
library(ggrepel)
library(sva)

yehor_dir <- "data/mapped/yehor/2026_04_27_9_preprocessed_datasets"
output_dir <- "output/one_off/pca_3_yehor_gse37653/combat_then_sex"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

exprs <- read.table(
  file.path(yehor_dir, "GSE37653_entrez_protein_coding.tsv"),
  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE
)

pdata <- read.table(
  file.path(yehor_dir, "phenodata_placenta_1_2.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)
pdata <- pdata[pdata$dataset_id == "GSE37653", ]
pdata <- pdata[match(colnames(exprs), pdata$sample_id), ]

batch <- as.factor(pdata$scan_batch)

cat("=== Step 1: ComBat WITHOUT sex covariate ===\n")
exprs_combat <- ComBat(
  dat = as.matrix(exprs), batch = batch, mod = NULL,
  ref.batch = "Dec 2009 Singapore v2.5"
)

cat("\n=== Step 2: Check sex-linked genes after batch correction ===\n")

# XIST (Entrez 7503) - female marker, silences one X chromosome
# RPS4Y1 (Entrez 6192) - Y-chromosome ribosomal protein
# DDX3Y (Entrez 8653) - Y-chromosome DEAD-box helicase
# EIF1AY (Entrez 9086) - Y-chromosome translation initiation factor
# KDM5D (Entrez 8284) - Y-chromosome lysine demethylase
sex_genes <- c(XIST = "7503", RPS4Y1 = "6192", DDX3Y = "8653", EIF1AY = "9086", KDM5D = "8284")

cat("\nSex-linked gene expression (after batch correction):\n")
for (gene_name in names(sex_genes)) {
  eid <- sex_genes[gene_name]
  if (eid %in% rownames(exprs_combat)) {
    vals <- exprs_combat[eid, ]
    cat(sprintf("  %s (Entrez %s): mean=%.2f, range=[%.2f, %.2f]\n",
                gene_name, eid, mean(vals), min(vals), max(vals)))
  } else {
    cat(sprintf("  %s (Entrez %s): not found\n", gene_name, eid))
  }
}

# Re-estimate sex from Y-linked genes
y_genes <- sex_genes[c("RPS4Y1", "DDX3Y", "EIF1AY", "KDM5D")]
y_found <- y_genes[y_genes %in% rownames(exprs_combat)]
cat(sprintf("\nUsing %d Y-linked genes for sex estimation: %s\n",
            length(y_found), paste(names(y_found), collapse = ", ")))

if (length(y_found) > 0) {
  y_score_raw <- colMeans(as.matrix(exprs[y_found, ]))
  y_score_combat <- colMeans(exprs_combat[y_found, ])

  # Approach 2: global threshold on post-ComBat data
  sex_postcombat <- ifelse(y_score_combat > median(y_score_combat), "m", "f")

  # Approach 3: per-batch threshold on raw data, then ComBat with sex covariate
  sex_perbatch <- character(length(y_score_raw))
  names(sex_perbatch) <- names(y_score_raw)
  for (b in levels(batch)) {
    idx <- which(pdata$scan_batch == b)
    scores_b <- y_score_raw[idx]
    sex_perbatch[idx] <- ifelse(scores_b > median(scores_b), "m", "f")
  }

  pdata$sex_postcombat <- sex_postcombat
  pdata$sex_perbatch <- sex_perbatch

  cat("\n=== Step 3: ComBat with per-batch sex estimate as covariate ===\n")
  sex_pb_factor <- as.factor(sex_perbatch)
  mod_perbatch <- model.matrix(~ sex_pb_factor)
  exprs_combat_pb <- ComBat(
    dat = as.matrix(exprs), batch = batch, mod = mod_perbatch,
    ref.batch = "Dec 2009 Singapore v2.5"
  )

  cat("\nY-gene score per sample — 3 approaches:\n")
  cat(sprintf("%-12s  %6s  %6s  %8s  %10s  %10s  %-25s\n",
              "Sample", "Raw", "ComBat", "Original", "PostComBat", "PerBatch", "Batch"))
  for (i in seq_along(y_score_raw)) {
    s <- names(y_score_raw)[i]
    cat(sprintf("%-12s  %6.2f  %6.2f  %8s  %10s  %10s  %-25s\n",
                s, y_score_raw[i], y_score_combat[i],
                pdata$fetux_sex_estimate[i], sex_postcombat[i], sex_perbatch[i],
                pdata$scan_batch[i]))
  }

  cat("\n--- Agreement between approaches ---\n")
  cat(sprintf("Original vs PostComBat:  %d / %d differ\n",
              sum(pdata$fetux_sex_estimate != sex_postcombat), nrow(pdata)))
  cat(sprintf("Original vs PerBatch:    %d / %d differ\n",
              sum(pdata$fetux_sex_estimate != sex_perbatch), nrow(pdata)))
  cat(sprintf("PostComBat vs PerBatch:  %d / %d differ\n",
              sum(sex_postcombat != sex_perbatch), nrow(pdata)))

  cat("\nSex distribution per batch — Original:\n")
  print(table(pdata$scan_batch, pdata$fetux_sex_estimate))
  cat("\nSex distribution per batch — PostComBat:\n")
  print(table(pdata$scan_batch, pdata$sex_postcombat))
  cat("\nSex distribution per batch — PerBatch:\n")
  print(table(pdata$scan_batch, pdata$sex_perbatch))
}

save_plot <- function(p, name) {
  svg(file.path(output_dir, paste0(name, ".svg")), width = 10, height = 8)
  print(p)
  dev.off()
  png(file.path(output_dir, paste0(name, ".png")), width = 10, height = 8, units = "in", res = 300)
  print(p)
  dev.off()
  cat("Saved:", name, "\n")
}

make_sex_pca <- function(data, sex_col, sex_label, title_text, var_exp) {
  ggplot(data, aes(x = PC1, y = PC2, color = .data[[sex_col]],
                   shape = batch, label = sample)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text_repel(size = 2.5, max.overlaps = 20, segment.size = 0.2) +
    labs(
      title = title_text,
      x = sprintf("PC1 (%.1f%%)", var_exp[1]),
      y = sprintf("PC2 (%.1f%%)", var_exp[2]),
      color = sex_label, shape = "Batch"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          legend.text = element_text(size = 9))
}

# --- Approach 1: ComBat with original sex (from earlier script) ---
cat("\n=== PCA — Approach 1: ComBat + original sex covariate ===\n")
mod_orig <- model.matrix(~ as.factor(pdata$fetux_sex_estimate))
exprs_combat_orig <- ComBat(
  dat = as.matrix(exprs), batch = batch, mod = mod_orig,
  ref.batch = "Dec 2009 Singapore v2.5"
)
pca1 <- prcomp(t(exprs_combat_orig), center = TRUE, scale. = FALSE)
var1 <- summary(pca1)$importance[2, 1:5] * 100
cat("Variance explained:\n")
for (i in 1:5) cat(sprintf("  PC%d: %.1f%%\n", i, var1[i]))

pc1 <- data.frame(PC1 = pca1$x[, 1], PC2 = pca1$x[, 2],
                  sex = pdata$fetux_sex_estimate,
                  batch = pdata$scan_batch, sample = pdata$sample_id)
save_plot(make_sex_pca(pc1, "sex", "Sex",
                       "1. ComBat + original sex covariate", var1),
          "pca_1_combat_original_sex")

# --- Approach 2: ComBat without sex, re-estimate globally ---
cat("\n=== PCA — Approach 2: ComBat no sex, global re-estimate ===\n")
pca2 <- prcomp(t(exprs_combat), center = TRUE, scale. = FALSE)
var2 <- summary(pca2)$importance[2, 1:5] * 100
cat("Variance explained:\n")
for (i in 1:5) cat(sprintf("  PC%d: %.1f%%\n", i, var2[i]))

pc2 <- data.frame(PC1 = pca2$x[, 1], PC2 = pca2$x[, 2],
                  sex = pdata$sex_postcombat,
                  batch = pdata$scan_batch, sample = pdata$sample_id)
save_plot(make_sex_pca(pc2, "sex", "Sex",
                       "2. ComBat no sex prot., global re-estimate", var2),
          "pca_2_combat_nosex_global")

# --- Approach 3: per-batch sex estimate, then ComBat with sex covariate ---
cat("\n=== PCA — Approach 3: per-batch sex estimate, then ComBat + sex ===\n")
pca3 <- prcomp(t(exprs_combat_pb), center = TRUE, scale. = FALSE)
var3 <- summary(pca3)$importance[2, 1:5] * 100
cat("Variance explained:\n")
for (i in 1:5) cat(sprintf("  PC%d: %.1f%%\n", i, var3[i]))

pc3 <- data.frame(PC1 = pca3$x[, 1], PC2 = pca3$x[, 2],
                  sex = pdata$sex_perbatch,
                  batch = pdata$scan_batch, sample = pdata$sample_id)
save_plot(make_sex_pca(pc3, "sex", "Sex",
                       "3. Per-batch sex est. → ComBat + sex covariate", var3),
          "pca_3_perbatch_then_combat")

# --- Summary comparison ---
cat("\n=== Variance explained comparison ===\n")
cat(sprintf("  Approach 1 (original sex):     PC1=%.1f%%  PC2=%.1f%%\n", var1[1], var1[2]))
cat(sprintf("  Approach 2 (no sex, global):   PC1=%.1f%%  PC2=%.1f%%\n", var2[1], var2[2]))
cat(sprintf("  Approach 3 (per-batch→ComBat): PC1=%.1f%%  PC2=%.1f%%\n", var3[1], var3[2]))

cat("\n=== Step 4: Compare Y-gene scores & 3 approaches ===\n")
diag <- data.frame(
  sample = names(y_score_raw),
  y_raw = as.numeric(y_score_raw),
  y_combat = as.numeric(y_score_combat),
  batch = pdata$scan_batch,
  sex_original = pdata$fetux_sex_estimate,
  sex_postcombat = pdata$sex_postcombat,
  sex_perbatch = pdata$sex_perbatch,
  stringsAsFactors = FALSE
)
diag <- diag[order(diag$y_raw), ]

cat("\nY-gene scores sorted by raw score (3 sex estimates):\n")
cat(sprintf("%-12s  %6s  %6s  %5s %5s %5s  %-25s\n",
            "Sample", "Raw", "ComBat", "Orig", "Post", "PerB", "Batch"))
for (i in seq_len(nrow(diag))) {
  agree <- diag$sex_original[i] == diag$sex_postcombat[i] &&
           diag$sex_postcombat[i] == diag$sex_perbatch[i]
  flag <- ifelse(agree, "", " ***")
  cat(sprintf("%-12s  %6.2f  %6.2f  %5s %5s %5s  %-25s%s\n",
              diag$sample[i], diag$y_raw[i], diag$y_combat[i],
              diag$sex_original[i], diag$sex_postcombat[i], diag$sex_perbatch[i],
              diag$batch[i], flag))
}
cat("(*** = approaches disagree)\n")

cat("\nPer-batch median Y-gene thresholds (raw):\n")
for (b in levels(batch)) {
  idx <- which(pdata$scan_batch == b)
  cat(sprintf("  %-30s  median=%.2f  n=%d\n", b, median(y_score_raw[idx]), length(idx)))
}

p3 <- ggplot(diag, aes(x = y_raw, y = y_combat, color = sex_original,
                        shape = batch, label = sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 2.5, max.overlaps = 25, segment.size = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    title = "Y-gene mean score: Raw vs ComBat (no sex protection)",
    x = "Y-gene score (raw)", y = "Y-gene score (after ComBat)",
    color = "Sex (original)", shape = "Batch"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
save_plot(p3, "y_score_raw_vs_combat")

diag$batch_short <- gsub(" v2\\.[56]", "", diag$batch)
p_strip <- ggplot(diag, aes(x = batch_short, y = y_raw, color = sex_perbatch)) +
  geom_jitter(width = 0.15, size = 3, alpha = 0.8) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5,
               color = "black", linewidth = 0.4) +
  labs(title = "Y-gene scores by batch (raw) — per-batch sex estimation",
       x = "Batch", y = "Y-gene mean score", color = "Sex (per-batch)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        axis.text.x = element_text(angle = 15, hjust = 1))
save_plot(p_strip, "y_score_strip_perbatch")

cat("\nAll plots saved to:", output_dir, "\n")

# ==========================================================================
# Step 5: Save approach 3 matrix and update phenodata
# ==========================================================================

cat("\n=== Step 5: Save corrected matrix and update phenodata ===\n")

# 5a. Rename old matrix, save approach 3 corrected matrix
old_path <- file.path(yehor_dir, "GSE37653_entrez_protein_coding.tsv")
backup_path <- file.path(yehor_dir, "GSE37653_entrez_protein_coding_old.tsv")
if (!file.exists(backup_path)) {
  file.rename(old_path, backup_path)
  cat("Renamed original matrix to _old.tsv\n")
} else {
  cat("Backup _old.tsv already exists, skipping rename\n")
}
write.table(exprs_combat_pb, old_path, sep = "\t", quote = FALSE)
cat(sprintf("Saved approach 3 matrix: %d genes x %d samples\n",
            nrow(exprs_combat_pb), ncol(exprs_combat_pb)))

# 5b. Build per-batch sex lookup (m/f)
sex_lookup <- setNames(sex_perbatch, names(sex_perbatch))
sex_lookup_MF <- ifelse(sex_lookup == "m", "Male", "Female")

# 5c. Update phenodata_placenta_1_2.tsv and _fix_sex.tsv
for (pheno_file in c("phenodata_placenta_1_2.tsv", "phenodata_placenta_1_2_fix_sex.tsv")) {
  pheno_path <- file.path(yehor_dir, pheno_file)
  ph <- read.table(pheno_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  changed <- 0
  for (i in seq_len(nrow(ph))) {
    sid <- ph$sample_id[i]
    if (sid %in% names(sex_lookup)) {
      old_val <- ph$fetux_sex_estimate[i]
      new_val <- sex_lookup[sid]
      if (old_val != new_val) {
        ph$fetux_sex_estimate[i] <- new_val
        changed <- changed + 1
      }
    }
  }
  write.table(ph, pheno_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Updated %s: %d sex values changed\n", pheno_file, changed))
}

# 5d. Update samples.csv
samples_path <- "data/phenodata/samples.csv"
samples <- read.csv(samples_path, stringsAsFactors = FALSE, check.names = FALSE)
changed_csv <- 0
for (i in seq_len(nrow(samples))) {
  sid <- samples[i, 1]
  if (sid %in% names(sex_lookup_MF)) {
    new_sex <- sex_lookup_MF[sid]
    old_est <- samples$Estimated.Fetus.Sex[i]
    if (!is.na(old_est) && old_est != "_" && old_est != new_sex) {
      samples$Estimated.Fetus.Sex[i] <- new_sex
      samples$estimated_sex[i] <- new_sex
      samples$Combined.Fetus.Sex[i] <- new_sex
      changed_csv <- changed_csv + 1
    }
  }
}
write.csv(samples, samples_path, row.names = FALSE)
cat(sprintf("Updated samples.csv: %d sex values changed\n", changed_csv))
