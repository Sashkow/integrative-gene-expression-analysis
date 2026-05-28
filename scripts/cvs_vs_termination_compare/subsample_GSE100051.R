#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(limma)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

out_dir <- "output/cvs_vs_termination_compare/cvs_vs_GSE100051/subsampling"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fdr_cutoff   <- 0.05
logfc_cutoff <- 1.0
ruv_k        <- 1L
n_iter       <- 30L

# ── Load data ────────────────────────────────────────────────────────────────

cvs_expr <- read.delim("data/mapped/cvs/GSE12767_entrez_protein_coding.tsv",
                        row.names = 1, check.names = FALSE)
abort_expr <- read.delim("data/mapped/GSE100051.tsv",
                          row.names = 1, check.names = FALSE)

pheno <- read.csv("data/phenodata/samples_cvs_ga_matched.csv",
                   stringsAsFactors = FALSE, quote = "\"")

cvs_pheno <- pheno[pheno$secondaryaccession == "GSE12767" &
                    pheno$Diagnosis == "Healthy", ]
cvs_samples <- cvs_pheno$arraydatafile_exprscolumnnames
cvs_expr <- cvs_expr[, colnames(cvs_expr) %in% cvs_samples]

abort_pheno <- pheno[pheno$secondaryaccession == "GSE100051" &
                      pheno$Diagnosis == "Healthy" &
                      pheno$Gestational.Age %in% c("10", "11", "12"), ]

week_samples <- list(
  "10" = abort_pheno$arraydatafile_exprscolumnnames[abort_pheno$Gestational.Age == "10"],
  "11" = abort_pheno$arraydatafile_exprscolumnnames[abort_pheno$Gestational.Age == "11"],
  "12" = abort_pheno$arraydatafile_exprscolumnnames[abort_pheno$Gestational.Age == "12"]
)
all_abort_samples <- abort_pheno$arraydatafile_exprscolumnnames

cat(sprintf("CVS: %d samples\n", length(cvs_samples)))
for (w in names(week_samples))
  cat(sprintf("GA %s: %d samples (%s)\n", w, length(week_samples[[w]]),
              paste(week_samples[[w]], collapse = ", ")))

control_genes <- readLines(
  "scripts/integrative_analysis/phase2b_direct_merge/config_cvs_vs_termination_compare/control_genes_cvs_GSE100051.txt")

# ── Helper: run RUV k=1 + limma on given samples ────────────────────────────

run_ruv_de <- function(cvs_mat, abort_mat, ctl_genes) {
  shared <- intersect(rownames(cvs_mat), rownames(abort_mat))
  merged <- cbind(cvs_mat[shared, , drop = FALSE],
                  abort_mat[shared, , drop = FALSE])

  gene_types <- suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db, keys = rownames(merged),
                          columns = "GENETYPE", keytype = "ENTREZID")
  )
  pc <- gene_types$ENTREZID[gene_types$GENETYPE == "protein-coding" & !is.na(gene_types$GENETYPE)]
  merged <- merged[rownames(merged) %in% pc, ]
  vars <- apply(merged, 1, var)
  merged <- merged[vars > 0.01, ]

  ctl_idx <- which(rownames(merged) %in% ctl_genes)
  if (length(ctl_idx) < 50) return(NULL)

  n_cvs <- ncol(cvs_mat)
  n_ab  <- ncol(abort_mat)
  group <- factor(c(rep("CVS", n_cvs), rep("Abort", n_ab)),
                  levels = c("Abort", "CVS"))

  Y_c <- t(as.matrix(merged[ctl_idx, , drop = FALSE]))
  Y_c <- scale(Y_c, center = TRUE, scale = FALSE)
  svd_c <- svd(Y_c, nu = ruv_k, nv = 0)
  W <- svd_c$u[, seq_len(ruv_k), drop = FALSE]

  design <- model.matrix(~ group + W)
  fit <- eBayes(lmFit(as.matrix(merged), design))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")

  sig <- tt$adj.P.Val < fdr_cutoff & abs(tt$logFC) > logfc_cutoff
  list(
    de       = tt,
    sig_genes = rownames(tt)[sig],
    n_sig    = sum(sig),
    n_up     = sum(sig & tt$logFC > 0),
    n_down   = sum(sig & tt$logFC < 0),
    n_tested = nrow(tt)
  )
}

jaccard <- function(a, b) {
  inter <- length(intersect(a, b))
  union_n <- length(union(a, b))
  if (union_n == 0) return(0)
  inter / union_n
}

# ── Reference: full dataset ──────────────────────────────────────────────────

cat("\n=== Reference: full GSE100051 (15 samples) ===\n")
ref_result <- run_ruv_de(
  cvs_expr,
  abort_expr[, colnames(abort_expr) %in% all_abort_samples],
  control_genes
)
cat(sprintf("Reference DEGs: %d (up=%d, down=%d) from %d genes\n",
            ref_result$n_sig, ref_result$n_up, ref_result$n_down, ref_result$n_tested))

# ══════════════════════════════════════════════════════════════════════════════
# TEST A: Per-week comparisons (deterministic)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== Test A: Per-week comparisons ===\n")

week_results <- list()
for (w in c("10", "11", "12")) {
  samples <- week_samples[[w]]
  abort_sub <- abort_expr[, colnames(abort_expr) %in% samples]
  res <- run_ruv_de(cvs_expr, abort_sub, control_genes)

  jac <- jaccard(res$sig_genes, ref_result$sig_genes)

  shared_genes <- intersect(rownames(res$de), rownames(ref_result$de))
  lfc_cor <- cor(res$de[shared_genes, "logFC"],
                 ref_result$de[shared_genes, "logFC"],
                 method = "spearman")

  shared_sig <- intersect(res$sig_genes, ref_result$sig_genes)
  if (length(shared_sig) > 1) {
    same_dir <- sum(sign(res$de[shared_sig, "logFC"]) ==
                    sign(ref_result$de[shared_sig, "logFC"]))
    dir_pct <- 100 * same_dir / length(shared_sig)
  } else {
    dir_pct <- NA
  }

  week_results[[w]] <- res
  cat(sprintf("  GA %s: %d DEGs (up=%d, down=%d), Jaccard=%.3f, logFC rho=%.3f, dir=%.1f%%\n",
              w, res$n_sig, res$n_up, res$n_down, jac, lfc_cor,
              ifelse(is.na(dir_pct), 0, dir_pct)))
}

# Cross-week agreement
cat("\n  Cross-week agreement:\n")
weeks <- c("10", "11", "12")
for (i in 1:(length(weeks) - 1)) {
  for (j in (i + 1):length(weeks)) {
    a <- week_results[[weeks[i]]]
    b <- week_results[[weeks[j]]]
    jac_ab <- jaccard(a$sig_genes, b$sig_genes)
    shared_sig <- intersect(a$sig_genes, b$sig_genes)
    shared_all <- intersect(rownames(a$de), rownames(b$de))
    lfc_cor <- cor(a$de[shared_all, "logFC"],
                   b$de[shared_all, "logFC"],
                   method = "spearman")
    if (length(shared_sig) > 1) {
      same_dir <- sum(sign(a$de[shared_sig, "logFC"]) ==
                      sign(b$de[shared_sig, "logFC"]))
      dir_pct <- 100 * same_dir / length(shared_sig)
    } else {
      dir_pct <- NA
    }
    cat(sprintf("    GA%s vs GA%s: Jaccard=%.3f, overlap=%d, logFC rho=%.3f, dir=%.1f%%\n",
                weeks[i], weeks[j], jac_ab, length(shared_sig), lfc_cor,
                ifelse(is.na(dir_pct), 0, dir_pct)))
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# TEST B: Week-balanced subsets (N per week = 1, 2, 3, 4, 5)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== Test B: Week-balanced subsets ===\n")

balanced_results <- data.frame()

for (n_per_week in 1:5) {
  if (n_per_week == 5) {
    n_iter_local <- 1L
  } else {
    n_iter_local <- n_iter
  }

  for (iter in seq_len(n_iter_local)) {
    set.seed(iter * 1000 + n_per_week)

    selected <- character(0)
    for (w in c("10", "11", "12")) {
      pool <- week_samples[[w]]
      selected <- c(selected, sample(pool, n_per_week))
    }

    abort_sub <- abort_expr[, colnames(abort_expr) %in% selected]
    res <- run_ruv_de(cvs_expr, abort_sub, control_genes)

    if (is.null(res)) {
      balanced_results <- rbind(balanced_results, data.frame(
        n_per_week = n_per_week, n_total = 3 * n_per_week,
        iter = iter, n_deg = NA, n_up = NA, n_down = NA,
        jaccard = NA, logfc_rho = NA, dir_pct = NA, status = "failed"
      ))
      next
    }

    jac <- jaccard(res$sig_genes, ref_result$sig_genes)
    shared_all <- intersect(rownames(res$de), rownames(ref_result$de))
    lfc_cor <- cor(res$de[shared_all, "logFC"],
                   ref_result$de[shared_all, "logFC"],
                   method = "spearman")

    shared_sig <- intersect(res$sig_genes, ref_result$sig_genes)
    if (length(shared_sig) > 1) {
      same_dir <- sum(sign(res$de[shared_sig, "logFC"]) ==
                      sign(ref_result$de[shared_sig, "logFC"]))
      dir_pct <- 100 * same_dir / length(shared_sig)
    } else {
      dir_pct <- NA
    }

    balanced_results <- rbind(balanced_results, data.frame(
      n_per_week = n_per_week, n_total = 3 * n_per_week,
      iter = iter, n_deg = res$n_sig, n_up = res$n_up, n_down = res$n_down,
      jaccard = jac, logfc_rho = lfc_cor, dir_pct = dir_pct, status = "ok"
    ))
  }

  ok <- balanced_results[balanced_results$n_per_week == n_per_week &
                          balanced_results$status == "ok", ]
  cat(sprintf("  N=%d/week (total=%d): median DEGs=%.0f, median Jaccard=%.3f, median rho=%.3f\n",
              n_per_week, 3 * n_per_week,
              median(ok$n_deg, na.rm = TRUE),
              median(ok$jaccard, na.rm = TRUE),
              median(ok$logfc_rho, na.rm = TRUE)))
}

write.csv(balanced_results, file.path(out_dir, "test_balanced_subsets.csv"),
          row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# TEST C: Random size sweep (any week mix)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== Test C: Random size sweep ===\n")

size_results <- data.frame()

for (n_abort in c(3, 5, 8, 10, 12, 15)) {
  if (n_abort == 15) {
    n_iter_local <- 1L
  } else {
    n_iter_local <- n_iter
  }

  for (iter in seq_len(n_iter_local)) {
    set.seed(iter * 3000 + n_abort)

    selected <- sample(all_abort_samples, n_abort)
    abort_sub <- abort_expr[, colnames(abort_expr) %in% selected]
    res <- run_ruv_de(cvs_expr, abort_sub, control_genes)

    if (is.null(res)) {
      size_results <- rbind(size_results, data.frame(
        n_abort = n_abort, iter = iter,
        n_deg = NA, n_up = NA, n_down = NA,
        jaccard = NA, logfc_rho = NA, dir_pct = NA,
        n_week10 = NA, n_week11 = NA, n_week12 = NA,
        status = "failed"
      ))
      next
    }

    jac <- jaccard(res$sig_genes, ref_result$sig_genes)
    shared_all <- intersect(rownames(res$de), rownames(ref_result$de))
    lfc_cor <- cor(res$de[shared_all, "logFC"],
                   ref_result$de[shared_all, "logFC"],
                   method = "spearman")

    shared_sig <- intersect(res$sig_genes, ref_result$sig_genes)
    if (length(shared_sig) > 1) {
      same_dir <- sum(sign(res$de[shared_sig, "logFC"]) ==
                      sign(ref_result$de[shared_sig, "logFC"]))
      dir_pct <- 100 * same_dir / length(shared_sig)
    } else {
      dir_pct <- NA
    }

    nw10 <- sum(selected %in% week_samples[["10"]])
    nw11 <- sum(selected %in% week_samples[["11"]])
    nw12 <- sum(selected %in% week_samples[["12"]])

    size_results <- rbind(size_results, data.frame(
      n_abort = n_abort, iter = iter,
      n_deg = res$n_sig, n_up = res$n_up, n_down = res$n_down,
      jaccard = jac, logfc_rho = lfc_cor, dir_pct = dir_pct,
      n_week10 = nw10, n_week11 = nw11, n_week12 = nw12,
      status = "ok"
    ))
  }

  ok <- size_results[size_results$n_abort == n_abort &
                      size_results$status == "ok", ]
  cat(sprintf("  N=%d: median DEGs=%.0f, median Jaccard=%.3f, median rho=%.3f\n",
              n_abort, median(ok$n_deg, na.rm = TRUE),
              median(ok$jaccard, na.rm = TRUE),
              median(ok$logfc_rho, na.rm = TRUE)))
}

write.csv(size_results, file.path(out_dir, "test_random_size_sweep.csv"),
          row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# TEST D: Single-week subsets vs each other
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== Test D: Week-specific results ===\n")

week_de_file <- file.path(out_dir, "per_week_degs.csv")
week_rows <- list()

for (w in c("10", "11", "12")) {
  res <- week_results[[w]]
  sig_genes <- res$sig_genes
  for (g in sig_genes) {
    week_rows[[length(week_rows) + 1]] <- data.frame(
      gene = g, week = w,
      logFC = res$de[g, "logFC"],
      adj.P.Val = res$de[g, "adj.P.Val"],
      stringsAsFactors = FALSE
    )
  }
}
if (length(week_rows) > 0) {
  week_deg_df <- do.call(rbind, week_rows)
  write.csv(week_deg_df, week_de_file, row.names = FALSE)
  cat(sprintf("  Per-week DEGs: GA10=%d, GA11=%d, GA12=%d\n",
              week_results[["10"]]$n_sig,
              week_results[["11"]]$n_sig,
              week_results[["12"]]$n_sig))

  all_week_genes <- unique(week_deg_df$gene)
  in_all_3 <- Reduce(intersect, lapply(week_results, function(x) x$sig_genes))
  in_any_2 <- character(0)
  for (i in 1:2) {
    for (j in (i + 1):3) {
      in_any_2 <- union(in_any_2, intersect(
        week_results[[weeks[i]]]$sig_genes,
        week_results[[weeks[j]]]$sig_genes))
    }
  }
  cat(sprintf("  In all 3 weeks: %d, in >=2 weeks: %d, any week: %d\n",
              length(in_all_3), length(in_any_2), length(all_week_genes)))
}

# ══════════════════════════════════════════════════════════════════════════════
# Summary plots
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== Generating plots ===\n")

# Plot 1: Balanced subset DEG count and Jaccard
png(file.path(out_dir, "balanced_subset_metrics.png"),
    width = 14, height = 5, units = "in", res = 200)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))

ok <- balanced_results[balanced_results$status == "ok", ]
boxplot(n_deg ~ n_per_week, data = ok,
        col = "#B3E5FC", border = "#1565C0",
        xlab = "Samples per week", ylab = "DEG count",
        main = "DEGs by week-balanced subset size")
abline(h = ref_result$n_sig, col = "red", lty = 2, lwd = 1.5)
text(0.6, ref_result$n_sig, sprintf("Full: %d", ref_result$n_sig),
     col = "red", pos = 3, cex = 0.8)

boxplot(jaccard ~ n_per_week, data = ok,
        col = "#C8E6C9", border = "#388E3C",
        xlab = "Samples per week", ylab = "Jaccard vs full",
        main = "Jaccard similarity to full dataset")

boxplot(logfc_rho ~ n_per_week, data = ok,
        col = "#FFE0B2", border = "#F57C00",
        xlab = "Samples per week", ylab = "Spearman rho",
        main = "logFC correlation with full dataset")

dev.off()
cat("  Saved: balanced_subset_metrics.png\n")

# Plot 2: Random size sweep
png(file.path(out_dir, "random_size_sweep.png"),
    width = 14, height = 5, units = "in", res = 200)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))

ok2 <- size_results[size_results$status == "ok", ]
boxplot(n_deg ~ n_abort, data = ok2,
        col = "#B3E5FC", border = "#1565C0",
        xlab = "N termination samples", ylab = "DEG count",
        main = "DEGs by random subset size")
abline(h = ref_result$n_sig, col = "red", lty = 2, lwd = 1.5)

boxplot(jaccard ~ n_abort, data = ok2,
        col = "#C8E6C9", border = "#388E3C",
        xlab = "N termination samples", ylab = "Jaccard vs full",
        main = "Jaccard similarity to full dataset")

boxplot(logfc_rho ~ n_abort, data = ok2,
        col = "#FFE0B2", border = "#F57C00",
        xlab = "N termination samples", ylab = "Spearman rho",
        main = "logFC correlation with full dataset")

dev.off()
cat("  Saved: random_size_sweep.png\n")

# Plot 3: Per-week DEG comparison
png(file.path(out_dir, "per_week_comparison.png"),
    width = 10, height = 6, units = "in", res = 200)
par(mar = c(5, 5, 4, 2))

week_counts <- sapply(week_results, function(x) c(x$n_up, x$n_down))
barplot(week_counts, beside = FALSE,
        col = c("#EF5350", "#42A5F5"),
        names.arg = paste("GA", c(10, 11, 12)),
        ylab = "DEG count", main = "DEGs per gestational week vs CVS (RUV k=1)",
        border = NA)
legend("topright", legend = c("Up in CVS", "Down in CVS"),
       fill = c("#EF5350", "#42A5F5"), border = NA)
abline(h = ref_result$n_sig, col = "grey30", lty = 2)
text(4, ref_result$n_sig, sprintf("Full (15 samples): %d", ref_result$n_sig),
     pos = 3, cex = 0.8)

dev.off()
cat("  Saved: per_week_comparison.png\n")

# ── Text summary ─────────────────────────────────────────────────────────────

sink(file.path(out_dir, "subsampling_summary.txt"))
cat("CVS vs GSE100051 Subsampling Analysis\n")
cat(sprintf("Date: %s\n\n", Sys.time()))

cat("Reference: full GSE100051 (15 samples, GA 10-12)\n")
cat(sprintf("  DEGs: %d (up=%d, down=%d), genes tested: %d\n\n",
            ref_result$n_sig, ref_result$n_up, ref_result$n_down, ref_result$n_tested))

cat("Test A: Per-week comparisons (each week vs all 8 CVS)\n")
for (w in c("10", "11", "12")) {
  res <- week_results[[w]]
  jac <- jaccard(res$sig_genes, ref_result$sig_genes)
  cat(sprintf("  GA %s: %d DEGs (up=%d, down=%d), Jaccard=%.3f\n",
              w, res$n_sig, res$n_up, res$n_down, jac))
}

cat("\nTest B: Week-balanced subsets (median values)\n")
for (n in 1:5) {
  ok <- balanced_results[balanced_results$n_per_week == n &
                          balanced_results$status == "ok", ]
  cat(sprintf("  %d/week (total=%d): DEGs=%.0f [%.0f-%.0f], Jaccard=%.3f, rho=%.3f\n",
              n, 3 * n,
              median(ok$n_deg, na.rm = TRUE),
              quantile(ok$n_deg, 0.25, na.rm = TRUE),
              quantile(ok$n_deg, 0.75, na.rm = TRUE),
              median(ok$jaccard, na.rm = TRUE),
              median(ok$logfc_rho, na.rm = TRUE)))
}

cat("\nTest C: Random size sweep (median values)\n")
for (n in c(3, 5, 8, 10, 12, 15)) {
  ok <- size_results[size_results$n_abort == n &
                      size_results$status == "ok", ]
  cat(sprintf("  N=%d: DEGs=%.0f [%.0f-%.0f], Jaccard=%.3f, rho=%.3f\n",
              n,
              median(ok$n_deg, na.rm = TRUE),
              quantile(ok$n_deg, 0.25, na.rm = TRUE),
              quantile(ok$n_deg, 0.75, na.rm = TRUE),
              median(ok$jaccard, na.rm = TRUE),
              median(ok$logfc_rho, na.rm = TRUE)))
}
sink()
cat("\nSaved: subsampling_summary.txt\n")
cat("Done.\n")
