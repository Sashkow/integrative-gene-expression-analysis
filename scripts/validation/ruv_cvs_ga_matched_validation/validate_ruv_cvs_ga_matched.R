#!/usr/bin/env Rscript
#' Validate RUV k=2 correction for GA-matched CVS vs abortion (3 datasets).
#'
#' Three validation tests:
#'   1. Permutation test — shuffle group labels, re-run RUV+limma, count DEGs
#'   2. P-value histogram — check for uniform + spike (good) vs skewed (bad)
#'   3. Negative control gene check — control genes should not be significant
#'
#' Usage: Rscript validate_ruv_cvs_ga_matched.R [n_permutations]
#'   Default: 200 permutations

suppressPackageStartupMessages({
  library(limma)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

args <- commandArgs(trailingOnly = TRUE)
n_perm <- if (length(args) >= 1) as.integer(args[1]) else 200L
ruv_k <- 2L

out_dir <- "output/validation/ruv_cvs_ga_matched_validation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fdr_cutoff <- 0.05
logfc_cutoff <- 1.0

# ── Load data ────────────────────────────────────────────────────────────────

cvs_expr <- read.delim("data/mapped/cvs/GSE12767_entrez_protein_coding.tsv",
                        row.names = 1, check.names = FALSE)
abort1_expr <- read.delim("data/mapped/GSE93520.tsv",
                           row.names = 1, check.names = FALSE)
abort2_expr <- read.delim("data/mapped/GSE100051.tsv",
                           row.names = 1, check.names = FALSE)

pheno <- read.csv("data/phenodata/samples_cvs_ga_matched.csv", stringsAsFactors = FALSE)

cvs_pheno <- pheno[pheno$secondaryaccession == "GSE12767" &
                    pheno$Diagnosis == "Healthy", ]
cvs_expr <- cvs_expr[, colnames(cvs_expr) %in% cvs_pheno$arraydatafile_exprscolumnnames]

abort1_pheno <- pheno[pheno$secondaryaccession == "GSE93520" &
                       pheno$Diagnosis == "Healthy" &
                       pheno$Gestational.Age %in% c("10"), ]
abort1_expr <- abort1_expr[, colnames(abort1_expr) %in% abort1_pheno$arraydatafile_exprscolumnnames]

abort2_pheno <- pheno[pheno$secondaryaccession == "GSE100051" &
                       pheno$Diagnosis == "Healthy" &
                       pheno$Gestational.Age %in% c("10", "11", "12"), ]
abort2_expr <- abort2_expr[, colnames(abort2_expr) %in% abort2_pheno$arraydatafile_exprscolumnnames]

shared <- Reduce(intersect, list(
  rownames(cvs_expr), rownames(abort1_expr), rownames(abort2_expr)
))
merged <- cbind(cvs_expr[shared, ], abort1_expr[shared, ], abort2_expr[shared, ])

gene_types <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db, keys = rownames(merged),
                        columns = "GENETYPE", keytype = "ENTREZID")
)
pc <- gene_types$ENTREZID[gene_types$GENETYPE == "protein-coding" & !is.na(gene_types$GENETYPE)]
merged <- merged[rownames(merged) %in% pc, ]
vars <- apply(merged, 1, var)
merged <- merged[vars > 0.01, ]

n_cvs <- ncol(cvs_expr)
n_abort1 <- ncol(abort1_expr)
n_abort2 <- ncol(abort2_expr)
n_abort <- n_abort1 + n_abort2
n_total <- n_cvs + n_abort

control_genes <- readLines(
  "scripts/integrative_analysis/phase2b_direct_merge/config_ruv_cvs_ga_matched/ruv_control_genes_ga_matched.txt")
control_idx <- which(rownames(merged) %in% control_genes)

cat(sprintf("Merged matrix: %d genes x %d samples (%d CVS, %d abortion [%d GSE93520 + %d GSE100051])\n",
            nrow(merged), ncol(merged), n_cvs, n_abort, n_abort1, n_abort2))
cat(sprintf("Control genes in matrix: %d\n", length(control_idx)))

# ── Helper: run RUV k + limma for a given group vector ───────────────────────

run_ruv_limma <- function(exprs, group_vec, ctl_idx, k = 2) {
  Y_c <- t(exprs[ctl_idx, , drop = FALSE])
  Y_c <- scale(Y_c, center = TRUE, scale = FALSE)
  svd_c <- svd(Y_c, nu = k, nv = 0)
  W <- svd_c$u[, seq_len(k), drop = FALSE]

  design <- model.matrix(~ group_vec + W)
  fit <- lmFit(exprs, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  tt
}

# ── Real result ──────────────────────────────────────────────────────────────

real_group <- factor(c(rep("CVS", n_cvs), rep("Abortion", n_abort)),
                     levels = c("Abortion", "CVS"))

cat(sprintf("\nRunning real analysis (RUV k=%d)...\n", ruv_k))
real_de <- run_ruv_limma(as.matrix(merged), real_group, control_idx, k = ruv_k)
real_sig <- sum(real_de$adj.P.Val < fdr_cutoff & abs(real_de$logFC) > logfc_cutoff)
real_up <- sum(real_de$adj.P.Val < fdr_cutoff & real_de$logFC > logfc_cutoff)
real_down <- sum(real_de$adj.P.Val < fdr_cutoff & real_de$logFC < -logfc_cutoff)
cat(sprintf("Real DEGs: %d (up=%d, down=%d)\n", real_sig, real_up, real_down))

# ══════════════════════════════════════════════════════════════════════════════
# TEST 1: Permutation test
# ══════════════════════════════════════════════════════════════════════════════

cat(sprintf("\n=== Test 1: Permutation test (%d permutations) ===\n", n_perm))

set.seed(42)
perm_degs <- integer(n_perm)
perm_up <- integer(n_perm)
perm_down <- integer(n_perm)

merged_mat <- as.matrix(merged)

for (i in seq_len(n_perm)) {
  if (i %% 20 == 0) cat(sprintf("  Permutation %d / %d\n", i, n_perm))
  perm_labels <- sample(real_group)
  perm_de <- run_ruv_limma(merged_mat, perm_labels, control_idx, k = ruv_k)
  sig <- perm_de$adj.P.Val < fdr_cutoff & abs(perm_de$logFC) > logfc_cutoff
  perm_degs[i] <- sum(sig)
  perm_up[i] <- sum(sig & perm_de$logFC > 0)
  perm_down[i] <- sum(sig & perm_de$logFC < 0)
}

perm_p <- (sum(perm_degs >= real_sig) + 1) / (n_perm + 1)

cat(sprintf("\nReal DEGs: %d\n", real_sig))
cat(sprintf("Permutation DEGs: mean=%.1f, median=%d, max=%d, sd=%.1f\n",
            mean(perm_degs), median(perm_degs), max(perm_degs), sd(perm_degs)))
cat(sprintf("Empirical p-value: %.4f\n", perm_p))
cat(sprintf("Fold enrichment: %.1fx over permutation mean\n", real_sig / max(mean(perm_degs), 1)))

png(file.path(out_dir, "permutation_test.png"),
    width = 10, height = 6, units = "in", res = 200)
par(mar = c(5, 5, 4, 2))
hist(perm_degs, breaks = 30, col = "grey80", border = "grey50",
     main = sprintf("Permutation test: RUV k=%d CVS vs Abortion, GA-matched (%d perms)", ruv_k, n_perm),
     xlab = "Number of DEGs (FDR<0.05, |logFC|>1)",
     ylab = "Frequency",
     xlim = c(0, max(real_sig * 1.1, max(perm_degs) * 1.1)))
abline(v = real_sig, col = "red", lwd = 2.5, lty = 1)
text(real_sig, par("usr")[4] * 0.9,
     sprintf("Real: %d DEGs\np = %.4f", real_sig, perm_p),
     col = "red", pos = 4, cex = 0.9, font = 2)
abline(v = mean(perm_degs), col = "blue", lwd = 1.5, lty = 2)
text(mean(perm_degs), par("usr")[4] * 0.7,
     sprintf("Permutation mean: %.0f", mean(perm_degs)),
     col = "blue", pos = 4, cex = 0.8)
dev.off()
cat("Saved: permutation_test.png\n")

png(file.path(out_dir, "permutation_test_direction.png"),
    width = 12, height = 5, units = "in", res = 200)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

hist(perm_up, breaks = 25, col = "#FFCDD2", border = "#E57373",
     main = "Permutation: UP-regulated DEGs",
     xlab = "Number of UP DEGs",
     xlim = c(0, max(real_up * 1.2, max(perm_up) * 1.1)))
abline(v = real_up, col = "red", lwd = 2.5)
text(real_up, par("usr")[4] * 0.85,
     sprintf("Real: %d", real_up), col = "red", pos = 4, font = 2)

hist(perm_down, breaks = 25, col = "#C8E6C9", border = "#81C784",
     main = "Permutation: DOWN-regulated DEGs",
     xlab = "Number of DOWN DEGs",
     xlim = c(0, max(real_down * 1.2, max(perm_down) * 1.1)))
abline(v = real_down, col = "red", lwd = 2.5)
text(real_down, par("usr")[4] * 0.85,
     sprintf("Real: %d", real_down), col = "red", pos = 4, font = 2)

dev.off()
cat("Saved: permutation_test_direction.png\n")

# ══════════════════════════════════════════════════════════════════════════════
# TEST 2: P-value histogram
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== Test 2: P-value histogram ===\n")

png(file.path(out_dir, "pvalue_histogram.png"),
    width = 10, height = 6, units = "in", res = 200)
par(mar = c(5, 5, 4, 2))
hist(real_de$P.Value, breaks = 50, col = "steelblue", border = "white",
     main = sprintf("P-value distribution (RUV k=%d, CVS vs Abortion, GA-matched)", ruv_k),
     xlab = "Raw p-value", ylab = "Frequency",
     sub = sprintf("%d genes tested", nrow(real_de)))
abline(h = nrow(real_de) / 50, col = "red", lwd = 1.5, lty = 2)
text(0.5, nrow(real_de) / 50,
     "Expected under null (uniform)", col = "red", pos = 3, cex = 0.8)
dev.off()
cat("Saved: pvalue_histogram.png\n")

n_small <- sum(real_de$P.Value < 0.05)
cat(sprintf("P < 0.05: %d / %d (%.1f%%, expected: 5%%)\n",
            n_small, nrow(real_de), 100 * n_small / nrow(real_de)))

pi0 <- min(1, 2 * mean(real_de$P.Value > 0.5))
cat(sprintf("Estimated pi0 (proportion true nulls): %.3f\n", pi0))
cat(sprintf("Estimated true positives: ~%d\n", round(nrow(real_de) * (1 - pi0))))

# ══════════════════════════════════════════════════════════════════════════════
# TEST 3: Negative control gene check
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== Test 3: Negative control gene check ===\n")

control_de <- real_de[rownames(real_de) %in% control_genes, ]
noncontrol_de <- real_de[!rownames(real_de) %in% control_genes, ]

ctl_sig <- sum(control_de$adj.P.Val < fdr_cutoff & abs(control_de$logFC) > logfc_cutoff)
nonctl_sig <- sum(noncontrol_de$adj.P.Val < fdr_cutoff & abs(noncontrol_de$logFC) > logfc_cutoff)

cat(sprintf("Control genes tested: %d\n", nrow(control_de)))
cat(sprintf("Control genes significant: %d (%.1f%%)\n",
            ctl_sig, 100 * ctl_sig / nrow(control_de)))
cat(sprintf("Non-control genes tested: %d\n", nrow(noncontrol_de)))
cat(sprintf("Non-control genes significant: %d (%.1f%%)\n",
            nonctl_sig, 100 * nonctl_sig / nrow(noncontrol_de)))

png(file.path(out_dir, "control_gene_pvalues.png"),
    width = 10, height = 6, units = "in", res = 200)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

hist(control_de$P.Value, breaks = 30, col = "#B3E5FC", border = "#4FC3F7",
     main = sprintf("P-values: control genes (n=%d)", nrow(control_de)),
     xlab = "Raw p-value", ylab = "Frequency")
abline(h = nrow(control_de) / 30, col = "red", lty = 2)
mtext(sprintf("Significant: %d (%.1f%%)", ctl_sig, 100 * ctl_sig / nrow(control_de)),
      side = 3, line = 0, cex = 0.8, col = "red")

hist(noncontrol_de$P.Value, breaks = 30, col = "#FFCCBC", border = "#FF8A65",
     main = sprintf("P-values: non-control genes (n=%d)", nrow(noncontrol_de)),
     xlab = "Raw p-value", ylab = "Frequency")
abline(h = nrow(noncontrol_de) / 30, col = "red", lty = 2)
mtext(sprintf("Significant: %d (%.1f%%)", nonctl_sig, 100 * nonctl_sig / nrow(noncontrol_de)),
      side = 3, line = 0, cex = 0.8, col = "red")

dev.off()
cat("Saved: control_gene_pvalues.png\n")

png(file.path(out_dir, "control_gene_logfc.png"),
    width = 8, height = 6, units = "in", res = 200)
par(mar = c(5, 5, 4, 2))
boxplot(list("Control genes" = control_de$logFC,
             "Non-control genes" = noncontrol_de$logFC),
        col = c("#B3E5FC", "#FFCCBC"),
        main = "logFC distribution: control vs non-control genes",
        ylab = "logFC (CVS vs Abortion)",
        outline = FALSE)
abline(h = 0, lty = 2, col = "grey50")
abline(h = c(-1, 1), lty = 3, col = "red")
dev.off()
cat("Saved: control_gene_logfc.png\n")

ks <- ks.test(control_de$P.Value, "punif")
cat(sprintf("KS test (control p-values vs uniform): D=%.3f, p=%.4f\n",
            ks$statistic, ks$p.value))
if (ks$p.value > 0.05) {
  cat("  -> Control gene p-values are consistent with uniform (GOOD)\n")
} else {
  cat("  -> Control gene p-values deviate from uniform (CONCERNING)\n")
}

# ── Summary ──────────────────────────────────────────────────────────────────

cat("\n")
cat("════════════════════════════════════════════════════════════════\n")
cat(sprintf("  RUV k=%d GA-matched Validation Summary\n", ruv_k))
cat("════════════════════════════════════════════════════════════════\n")
cat(sprintf("  Datasets: GSE12767 (CVS, %d) + GSE93520 (%d) + GSE100051 (%d)\n",
            n_cvs, n_abort1, n_abort2))
cat(sprintf("  1. Permutation test:   real=%d DEGs, perm mean=%.0f, p=%.4f\n",
            real_sig, mean(perm_degs), perm_p))
cat(sprintf("  2. Pi0 estimate:       %.3f (%.0f estimated true positives)\n",
            pi0, nrow(real_de) * (1 - pi0)))
cat(sprintf("  3. Control gene leak:  %d / %d (%.1f%%) significant\n",
            ctl_sig, nrow(control_de), 100 * ctl_sig / nrow(control_de)))
cat(sprintf("     KS test p-value:    %.4f\n", ks$p.value))
cat("════════════════════════════════════════════════════════════════\n")

sink(file.path(out_dir, "validation_summary.txt"))
cat(sprintf("RUV k=%d GA-matched CVS vs Abortion Validation\n", ruv_k))
cat(sprintf("Date: %s\n\n", Sys.time()))
cat(sprintf("Samples: %d CVS (GSE12767) vs %d abortion (%d GSE93520 + %d GSE100051)\n",
            n_cvs, n_abort, n_abort1, n_abort2))
cat(sprintf("Genes tested: %d\n", nrow(real_de)))
cat(sprintf("Control genes: %d\n\n", length(control_idx)))
cat(sprintf("Real DEGs: %d (up=%d, down=%d)\n\n", real_sig, real_up, real_down))
cat(sprintf("1. Permutation test (%d permutations)\n", n_perm))
cat(sprintf("   Real DEGs: %d\n", real_sig))
cat(sprintf("   Permutation mean: %.1f (sd=%.1f, max=%d)\n", mean(perm_degs), sd(perm_degs), max(perm_degs)))
cat(sprintf("   Empirical p-value: %.4f\n", perm_p))
cat(sprintf("   Fold enrichment: %.1fx\n\n", real_sig / max(mean(perm_degs), 1)))
cat(sprintf("2. P-value distribution\n"))
cat(sprintf("   P < 0.05: %d / %d (%.1f%%)\n", n_small, nrow(real_de), 100 * n_small / nrow(real_de)))
cat(sprintf("   Pi0 (true null proportion): %.3f\n", pi0))
cat(sprintf("   Estimated true positives: %d\n\n", round(nrow(real_de) * (1 - pi0))))
cat(sprintf("3. Negative control gene check\n"))
cat(sprintf("   Control significant: %d / %d (%.1f%%)\n", ctl_sig, nrow(control_de), 100 * ctl_sig / nrow(control_de)))
cat(sprintf("   Non-control significant: %d / %d (%.1f%%)\n", nonctl_sig, nrow(noncontrol_de), 100 * nonctl_sig / nrow(noncontrol_de)))
cat(sprintf("   KS test (control vs uniform): D=%.3f, p=%.4f\n", ks$statistic, ks$p.value))
sink()
cat("Saved: validation_summary.txt\n")
