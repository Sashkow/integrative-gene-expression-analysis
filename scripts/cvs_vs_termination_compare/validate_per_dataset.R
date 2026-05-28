#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(limma)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript validate_per_dataset.R <dataset> [n_permutations] [logfc_cutoff]")
dataset      <- args[1]
n_perm       <- if (length(args) >= 2) as.integer(args[2]) else 200L
logfc_cutoff <- if (length(args) >= 3) as.numeric(args[3]) else 1.0

fdr_cutoff <- 0.05
ruv_k      <- 1L

config_dir <- "scripts/integrative_analysis/phase2b_direct_merge/config_cvs_vs_termination_compare"
suffix     <- if (logfc_cutoff == 0) "_fdr_only" else ""
out_dir    <- sprintf("output/validation/cvs_vs_termination_compare/cvs_vs_%s%s", dataset, suffix)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

thresh_label <- if (logfc_cutoff > 0) {
  sprintf("FDR<%.2f, |logFC|>%.1f", fdr_cutoff, logfc_cutoff)
} else {
  sprintf("FDR<%.2f only", fdr_cutoff)
}

# ── Load CVS data ────────────────────────────────────────────────────────────
cvs_expr <- read.delim("data/mapped/cvs/GSE12767_entrez_protein_coding.tsv",
                        row.names = 1, check.names = FALSE)
pheno <- read.csv("data/phenodata/samples_cvs_ga_matched.csv",
                   stringsAsFactors = FALSE, quote = "\"")

cvs_pheno <- pheno[pheno$secondaryaccession == "GSE12767" &
                    pheno$Diagnosis == "Healthy", ]
cvs_expr <- cvs_expr[, colnames(cvs_expr) %in% cvs_pheno$arraydatafile_exprscolumnnames]

# ── Load termination dataset ─────────────────────────────────────────────────
file_map <- list(
  GSE100051 = "data/mapped/GSE100051.tsv",
  GSE93520  = "data/mapped/GSE93520.tsv",
  GSE28551  = "data/mapped/GSE28551.tsv"
)
ga_filter <- list(
  GSE100051 = c("10", "11", "12"),
  GSE93520  = c("10"),
  GSE28551  = NULL
)

abort_expr <- read.delim(file_map[[dataset]], row.names = 1, check.names = FALSE)

if (!is.null(ga_filter[[dataset]])) {
  abort_pheno <- pheno[pheno$secondaryaccession == dataset &
                        pheno$Diagnosis == "Healthy" &
                        pheno$Gestational.Age %in% ga_filter[[dataset]], ]
} else {
  abort_pheno <- pheno[pheno$secondaryaccession == dataset &
                        pheno$Diagnosis == "Healthy" &
                        pheno$Gestational.Age.Category == "First Trimester", ]
}
abort_expr <- abort_expr[, colnames(abort_expr) %in% abort_pheno$arraydatafile_exprscolumnnames]

shared <- intersect(rownames(cvs_expr), rownames(abort_expr))
merged <- cbind(cvs_expr[shared, ], abort_expr[shared, ])

gene_types <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db, keys = rownames(merged),
                        columns = "GENETYPE", keytype = "ENTREZID")
)
pc <- gene_types$ENTREZID[gene_types$GENETYPE == "protein-coding" & !is.na(gene_types$GENETYPE)]
merged <- merged[rownames(merged) %in% pc, ]
vars <- apply(merged, 1, var)
merged <- merged[vars > 0.01, ]

n_cvs   <- ncol(cvs_expr)
n_abort <- ncol(abort_expr)

control_genes <- readLines(file.path(config_dir,
  sprintf("control_genes_cvs_%s.txt", dataset)))
control_idx <- which(rownames(merged) %in% control_genes)

cat(sprintf("=== CVS vs %s validation ===\n", dataset))
cat(sprintf("Matrix: %d genes x %d samples (%d CVS, %d termination)\n",
            nrow(merged), ncol(merged), n_cvs, n_abort))
cat(sprintf("Control genes: %d\n", length(control_idx)))

# ── RUV + limma helper ───────────────────────────────────────────────────────
run_ruv_limma <- function(exprs, group_vec, ctl_idx, k = 1) {
  Y_c <- t(exprs[ctl_idx, , drop = FALSE])
  Y_c <- scale(Y_c, center = TRUE, scale = FALSE)
  svd_c <- svd(Y_c, nu = k, nv = 0)
  W <- svd_c$u[, seq_len(k), drop = FALSE]
  design <- model.matrix(~ group_vec + W)
  fit <- lmFit(exprs, design)
  fit <- eBayes(fit)
  topTable(fit, coef = 2, number = Inf, sort.by = "none")
}

# ── Real result ──────────────────────────────────────────────────────────────
real_group <- factor(c(rep("CVS", n_cvs), rep("Abortion", n_abort)),
                     levels = c("Abortion", "CVS"))

real_de <- run_ruv_limma(as.matrix(merged), real_group, control_idx, k = ruv_k)
real_sig  <- sum(real_de$adj.P.Val < fdr_cutoff & abs(real_de$logFC) > logfc_cutoff)
real_up   <- sum(real_de$adj.P.Val < fdr_cutoff & real_de$logFC > logfc_cutoff)
real_down <- sum(real_de$adj.P.Val < fdr_cutoff & real_de$logFC < -logfc_cutoff)
cat(sprintf("Real DEGs: %d (up=%d, down=%d)\n", real_sig, real_up, real_down))

# ── Permutation test ─────────────────────────────────────────────────────────
cat(sprintf("\nPermutation test (%d permutations)...\n", n_perm))
set.seed(42)
perm_degs <- integer(n_perm)
merged_mat <- as.matrix(merged)

for (i in seq_len(n_perm)) {
  if (i %% 50 == 0) cat(sprintf("  %d / %d\n", i, n_perm))
  perm_labels <- sample(real_group)
  perm_de <- run_ruv_limma(merged_mat, perm_labels, control_idx, k = ruv_k)
  perm_degs[i] <- sum(perm_de$adj.P.Val < fdr_cutoff & abs(perm_de$logFC) > logfc_cutoff)
}

perm_p <- (sum(perm_degs >= real_sig) + 1) / (n_perm + 1)
cat(sprintf("Permutation mean: %.1f, max: %d\n", mean(perm_degs), max(perm_degs)))
cat(sprintf("Empirical p-value: %.4f\n", perm_p))
cat(sprintf("Fold enrichment: %.1fx\n", real_sig / max(mean(perm_degs), 1)))

# ── Permutation plot ─────────────────────────────────────────────────────────
png(file.path(out_dir, "permutation_test.png"),
    width = 10, height = 6, units = "in", res = 200)
par(mar = c(5, 5, 4, 2))
hist(perm_degs, breaks = 30, col = "grey80", border = "grey50",
     main = sprintf("Permutation: CVS vs %s (RUV k=%d, %d perms)",
                    dataset, ruv_k, n_perm),
     xlab = sprintf("Number of DEGs (%s)", thresh_label),
     xlim = c(0, max(real_sig * 1.1, max(perm_degs) * 1.1, 10)))
abline(v = real_sig, col = "red", lwd = 2.5)
text(real_sig, par("usr")[4] * 0.9,
     sprintf("Real: %d DEGs\np = %.4f", real_sig, perm_p),
     col = "red", pos = 4, cex = 0.9, font = 2)
dev.off()

# ── P-value histogram ────────────────────────────────────────────────────────
png(file.path(out_dir, "pvalue_histogram.png"),
    width = 10, height = 6, units = "in", res = 200)
par(mar = c(5, 5, 4, 2))
hist(real_de$P.Value, breaks = 50, col = "steelblue", border = "white",
     main = sprintf("P-value distribution: CVS vs %s (RUV k=%d)", dataset, ruv_k),
     xlab = "Raw p-value")
abline(h = nrow(real_de) / 50, col = "red", lwd = 1.5, lty = 2)
dev.off()

n_small <- sum(real_de$P.Value < 0.05)
pi0 <- min(1, 2 * mean(real_de$P.Value > 0.5))
cat(sprintf("P < 0.05: %d / %d (%.1f%%)\n", n_small, nrow(real_de),
            100 * n_small / nrow(real_de)))
cat(sprintf("Pi0: %.3f, estimated true positives: ~%d\n",
            pi0, round(nrow(real_de) * (1 - pi0))))

# ── Control gene check ───────────────────────────────────────────────────────
control_de    <- real_de[rownames(real_de) %in% control_genes, ]
noncontrol_de <- real_de[!rownames(real_de) %in% control_genes, ]
ctl_sig    <- sum(control_de$adj.P.Val < fdr_cutoff & abs(control_de$logFC) > logfc_cutoff)
nonctl_sig <- sum(noncontrol_de$adj.P.Val < fdr_cutoff & abs(noncontrol_de$logFC) > logfc_cutoff)

cat(sprintf("Control significant: %d / %d (%.1f%%)\n",
            ctl_sig, nrow(control_de), 100 * ctl_sig / nrow(control_de)))
cat(sprintf("Non-control significant: %d / %d (%.1f%%)\n",
            nonctl_sig, nrow(noncontrol_de), 100 * nonctl_sig / nrow(noncontrol_de)))

png(file.path(out_dir, "control_gene_pvalues.png"),
    width = 10, height = 6, units = "in", res = 200)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
hist(control_de$P.Value, breaks = 30, col = "#B3E5FC", border = "#4FC3F7",
     main = sprintf("Control genes (n=%d)", nrow(control_de)),
     xlab = "Raw p-value")
abline(h = nrow(control_de) / 30, col = "red", lty = 2)
hist(noncontrol_de$P.Value, breaks = 30, col = "#FFCCBC", border = "#FF8A65",
     main = sprintf("Non-control genes (n=%d)", nrow(noncontrol_de)),
     xlab = "Raw p-value")
abline(h = nrow(noncontrol_de) / 30, col = "red", lty = 2)
dev.off()

ks <- ks.test(control_de$P.Value, "punif")
cat(sprintf("KS test (control vs uniform): D=%.3f, p=%.4f\n", ks$statistic, ks$p.value))

# ── Write summary ────────────────────────────────────────────────────────────
sink(file.path(out_dir, "validation_summary.txt"))
cat(sprintf("CVS vs %s — RUV k=%d Validation (%s)\n", dataset, ruv_k, thresh_label))
cat(sprintf("Date: %s\n\n", Sys.time()))
cat(sprintf("Threshold: %s\n", thresh_label))
cat(sprintf("Samples: %d CVS vs %d termination (%s)\n", n_cvs, n_abort, dataset))
cat(sprintf("Genes tested: %d\n", nrow(real_de)))
cat(sprintf("Control genes: %d\n\n", length(control_idx)))
cat(sprintf("Real DEGs: %d (up=%d, down=%d)\n\n", real_sig, real_up, real_down))
cat(sprintf("1. Permutation test (%d permutations)\n", n_perm))
cat(sprintf("   Permutation mean: %.1f (sd=%.1f, max=%d)\n",
            mean(perm_degs), sd(perm_degs), max(perm_degs)))
cat(sprintf("   Empirical p-value: %.4f\n", perm_p))
cat(sprintf("   Fold enrichment: %.1fx\n\n", real_sig / max(mean(perm_degs), 1)))
cat(sprintf("2. P-value distribution\n"))
cat(sprintf("   P < 0.05: %d / %d (%.1f%%)\n", n_small, nrow(real_de),
            100 * n_small / nrow(real_de)))
cat(sprintf("   Pi0: %.3f\n", pi0))
cat(sprintf("   Estimated true positives: %d\n\n", round(nrow(real_de) * (1 - pi0))))
cat(sprintf("3. Negative control gene check\n"))
cat(sprintf("   Control significant: %d / %d (%.1f%%)\n",
            ctl_sig, nrow(control_de), 100 * ctl_sig / nrow(control_de)))
cat(sprintf("   Non-control significant: %d / %d (%.1f%%)\n",
            nonctl_sig, nrow(noncontrol_de), 100 * nonctl_sig / nrow(noncontrol_de)))
cat(sprintf("   KS test: D=%.3f, p=%.4f\n", ks$statistic, ks$p.value))
sink()
cat(sprintf("\nSaved to: %s\n", out_dir))
