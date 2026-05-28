#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

out_dir <- "output/cvs_vs_termination_compare"

# ── Load DEG lists ───────────────────────────────────────────────────────────

datasets <- c("GSE100051", "GSE93520", "GSE28551")
deg_lists <- list()
full_de   <- list()

for (ds in datasets) {
  f <- file.path(out_dir, sprintf("cvs_vs_%s", ds), "difexp_significant_none_ruv.tsv")
  de <- read.delim(f, stringsAsFactors = FALSE)
  deg_lists[[ds]] <- rownames(de)
  full_de[[ds]]   <- de
  cat(sprintf("%s: %d DEGs (up=%d, down=%d)\n",
              ds, nrow(de), sum(de$logFC > 0), sum(de$logFC < 0)))
}

# Load 3-dataset reference
ref_file <- "output/phase2b_ruv_cvs_ga_matched/cvs_vs_abortion_ruv_k2_ga_matched/difexp_significant_none_ruv.tsv"
if (file.exists(ref_file)) {
  ref_de <- read.delim(ref_file, stringsAsFactors = FALSE)
  deg_lists[["3ds_k2"]] <- rownames(ref_de)
  full_de[["3ds_k2"]]   <- ref_de
  cat(sprintf("3ds_k2: %d DEGs (up=%d, down=%d)\n",
              nrow(ref_de), sum(ref_de$logFC > 0), sum(ref_de$logFC < 0)))
}

all_names <- names(deg_lists)

# ── Jaccard similarity ───────────────────────────────────────────────────────
jaccard <- function(a, b) {
  inter <- length(intersect(a, b))
  union <- length(union(a, b))
  if (union == 0) return(0)
  inter / union
}

n <- length(all_names)
jac_mat <- matrix(0, n, n, dimnames = list(all_names, all_names))
overlap_mat <- matrix(0, n, n, dimnames = list(all_names, all_names))

for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    jac_mat[i, j] <- jaccard(deg_lists[[i]], deg_lists[[j]])
    overlap_mat[i, j] <- length(intersect(deg_lists[[i]], deg_lists[[j]]))
  }
}

cat("\n=== Jaccard similarity ===\n")
print(round(jac_mat, 3))

cat("\n=== Overlap counts ===\n")
print(overlap_mat)

# ── Direction concordance ────────────────────────────────────────────────────
cat("\n=== Direction concordance (shared DEGs) ===\n")
for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    shared <- intersect(deg_lists[[i]], deg_lists[[j]])
    if (length(shared) < 2) {
      cat(sprintf("%s vs %s: %d shared (too few)\n",
                  all_names[i], all_names[j], length(shared)))
      next
    }
    lfc_i <- full_de[[i]][shared, "logFC"]
    lfc_j <- full_de[[j]][shared, "logFC"]
    same_dir <- sum(sign(lfc_i) == sign(lfc_j))
    cat(sprintf("%s vs %s: %d shared, %d same direction (%.1f%%), cor=%.3f\n",
                all_names[i], all_names[j], length(shared), same_dir,
                100 * same_dir / length(shared),
                cor(lfc_i, lfc_j, method = "spearman")))
  }
}

# ── Venn-style breakdown (3 independent comparisons) ─────────────────────────
a <- deg_lists[["GSE100051"]]
b <- deg_lists[["GSE93520"]]
c <- deg_lists[["GSE28551"]]

only_a <- setdiff(setdiff(a, b), c)
only_b <- setdiff(setdiff(b, a), c)
only_c <- setdiff(setdiff(c, a), b)
ab_only <- setdiff(intersect(a, b), c)
ac_only <- setdiff(intersect(a, c), b)
bc_only <- setdiff(intersect(b, c), a)
abc <- Reduce(intersect, list(a, b, c))

cat("\n=== Venn breakdown (3 independent comparisons) ===\n")
cat(sprintf("GSE100051 only: %d\n", length(only_a)))
cat(sprintf("GSE93520 only:  %d\n", length(only_b)))
cat(sprintf("GSE28551 only:  %d\n", length(only_c)))
cat(sprintf("GSE100051 & GSE93520:  %d\n", length(ab_only)))
cat(sprintf("GSE100051 & GSE28551:  %d\n", length(ac_only)))
cat(sprintf("GSE93520 & GSE28551:   %d\n", length(bc_only)))
cat(sprintf("All three:             %d\n", length(abc)))

# ── Summary table ────────────────────────────────────────────────────────────
summary_df <- data.frame(
  Comparison = c("CVS vs GSE100051", "CVS vs GSE93520", "CVS vs GSE28551", "3-dataset k=2"),
  Platform = c("Illumina HT-12", "Agilent 4x44K", "ABI HG v2", "Mixed (3)"),
  n_CVS = c(8, 8, 8, 8),
  n_Termination = c(15, 5, 16, 20),
  n_genes_tested = c(16434, 16118, 9629, 15712),
  n_DEGs = sapply(all_names, function(x) length(deg_lists[[x]])),
  n_up = sapply(all_names, function(x) sum(full_de[[x]]$logFC > 0)),
  n_down = sapply(all_names, function(x) sum(full_de[[x]]$logFC < 0)),
  stringsAsFactors = FALSE
)
rownames(summary_df) <- NULL

write.csv(summary_df, file.path(out_dir, "summary_table.csv"), row.names = FALSE)
cat("\n=== Summary table ===\n")
print(summary_df)

# ── Top GO overlap ───────────────────────────────────────────────────────────
cat("\n=== Shared GO BP terms (down in CVS, top 10 each) ===\n")
go_terms <- list()
for (ds in datasets) {
  gf <- file.path(out_dir, sprintf("cvs_vs_%s", ds), "go_bp_down_in_cvs.csv")
  if (file.exists(gf)) {
    go <- read.csv(gf, stringsAsFactors = FALSE)
    if ("Description" %in% colnames(go))
      go_terms[[ds]] <- go$Description
  }
}
if (length(go_terms) >= 2) {
  shared_go <- Reduce(intersect, go_terms)
  cat(sprintf("Shared across all datasets with GO: %d terms\n", length(shared_go)))
  if (length(shared_go) > 0)
    cat(paste("  -", head(shared_go, 15)), sep = "\n")
}

# ── FDR-only sensitivity analysis ────────────────────────────────────────────
cat("\n=== FDR-only sensitivity analysis ===\n")
fdr_only_info <- list()

for (ds in datasets) {
  f <- file.path(out_dir, sprintf("cvs_vs_%s", ds), "difexp_none_ruv.tsv")
  de_all <- read.delim(f, stringsAsFactors = FALSE)
  de_fdr <- de_all[!is.na(de_all$adj.P.Val) & de_all$adj.P.Val < 0.05, ]
  min_lfc <- min(abs(de_fdr$logFC))
  n_fdr <- nrow(de_fdr)
  n_logfc <- length(deg_lists[[ds]])
  fdr_only_info[[ds]] <- list(
    n_fdr_only = n_fdr, n_fdr_logfc = n_logfc, min_logfc = min_lfc,
    identical = (n_fdr == n_logfc)
  )
  cat(sprintf("%s: FDR-only=%d, FDR+logFC=%d, identical=%s, min|logFC|=%.2f\n",
              ds, n_fdr, n_logfc, n_fdr == n_logfc, min_lfc))
}

ref_all_file <- "output/phase2b_ruv_cvs_ga_matched/cvs_vs_abortion_ruv_k2_ga_matched/difexp_none_ruv.tsv"
if (file.exists(ref_all_file)) {
  ref_all <- read.delim(ref_all_file, stringsAsFactors = FALSE)
  ref_fdr <- ref_all[!is.na(ref_all$adj.P.Val) & ref_all$adj.P.Val < 0.05, ]
  fdr_only_info[["3ds_k2"]] <- list(
    n_fdr_only = nrow(ref_fdr), n_fdr_logfc = length(deg_lists[["3ds_k2"]]),
    min_logfc = min(abs(ref_fdr$logFC)),
    identical = (nrow(ref_fdr) == length(deg_lists[["3ds_k2"]]))
  )
  cat(sprintf("3ds_k2: FDR-only=%d, FDR+logFC=%d, identical=%s, min|logFC|=%.2f\n",
              nrow(ref_fdr), length(deg_lists[["3ds_k2"]]),
              nrow(ref_fdr) == length(deg_lists[["3ds_k2"]]),
              min(abs(ref_fdr$logFC))))
}

fdr_summary <- data.frame(
  Comparison = summary_df$Comparison,
  n_DEGs_fdr_logfc = summary_df$n_DEGs,
  n_DEGs_fdr_only = sapply(all_names, function(x) fdr_only_info[[x]]$n_fdr_only),
  identical = sapply(all_names, function(x) fdr_only_info[[x]]$identical),
  min_abs_logFC = sapply(all_names, function(x) round(fdr_only_info[[x]]$min_logfc, 2)),
  stringsAsFactors = FALSE
)
rownames(fdr_summary) <- NULL
write.csv(fdr_summary, file.path(out_dir, "fdr_only_comparison.csv"), row.names = FALSE)
cat("\n=== FDR-only comparison ===\n")
print(fdr_summary)

# ── Load FDR-only validation results ────────────────────────────────────────
cat("\n=== FDR-only permutation validation ===\n")
for (ds in datasets) {
  val_file <- sprintf("output/validation/cvs_vs_termination_compare/cvs_vs_%s_fdr_only/validation_summary.txt",
                      ds)
  if (file.exists(val_file)) {
    lines <- readLines(val_file)
    perm_line <- grep("Empirical p-value", lines, value = TRUE)
    fold_line <- grep("Fold enrichment", lines, value = TRUE)
    cat(sprintf("%s: %s, %s\n", ds,
                trimws(perm_line), trimws(fold_line)))
  }
}

# ── Write text summary ───────────────────────────────────────────────────────
sink(file.path(out_dir, "comparison_summary.txt"))
cat("CVS vs Termination — Independent Dataset Comparisons\n")
cat(sprintf("Date: %s\n\n", Sys.time()))

cat("1. DEG counts per comparison\n")
for (nm in all_names) {
  de <- full_de[[nm]]
  cat(sprintf("   %s: %d DEGs (up=%d, down=%d)\n",
              nm, nrow(de), sum(de$logFC > 0), sum(de$logFC < 0)))
}

cat("\n2. Jaccard similarity\n")
for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    cat(sprintf("   %s vs %s: %.3f (%d shared)\n",
                all_names[i], all_names[j],
                jac_mat[i, j], overlap_mat[i, j]))
  }
}

cat("\n3. Venn breakdown (3 independent)\n")
cat(sprintf("   GSE100051 only: %d\n", length(only_a)))
cat(sprintf("   GSE93520 only:  %d\n", length(only_b)))
cat(sprintf("   GSE28551 only:  %d\n", length(only_c)))
cat(sprintf("   GSE100051 & GSE93520:  %d\n", length(ab_only)))
cat(sprintf("   GSE100051 & GSE28551:  %d\n", length(ac_only)))
cat(sprintf("   GSE93520 & GSE28551:   %d\n", length(bc_only)))
cat(sprintf("   All three:             %d\n", length(abc)))

cat("\n4. Direction concordance\n")
for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    shared <- intersect(deg_lists[[i]], deg_lists[[j]])
    if (length(shared) < 2) next
    lfc_i <- full_de[[i]][shared, "logFC"]
    lfc_j <- full_de[[j]][shared, "logFC"]
    same_dir <- sum(sign(lfc_i) == sign(lfc_j))
    cat(sprintf("   %s vs %s: %d/%d same direction (%.1f%%)\n",
                all_names[i], all_names[j], same_dir, length(shared),
                100 * same_dir / length(shared)))
  }
}

cat("\n5. FDR-only sensitivity analysis\n")
cat("   Removing the |logFC| > 1 threshold does not change DEG counts.\n")
cat("   All FDR-significant genes have extreme fold changes:\n")
for (nm in names(fdr_only_info)) {
  info <- fdr_only_info[[nm]]
  cat(sprintf("   %s: FDR-only=%d, FDR+logFC=%d, min|logFC|=%.2f, identical=%s\n",
              nm, info$n_fdr_only, info$n_fdr_logfc, info$min_logfc, info$identical))
}

cat("\n6. FDR-only permutation validation\n")
for (ds in datasets) {
  val_file <- sprintf("output/validation/cvs_vs_termination_compare/cvs_vs_%s_fdr_only/validation_summary.txt",
                      ds)
  if (file.exists(val_file)) {
    lines <- readLines(val_file)
    perm_line <- grep("Empirical p-value", lines, value = TRUE)
    perm_mean_line <- grep("Permutation mean", lines, value = TRUE)
    cat(sprintf("   %s: %s, %s\n", ds, trimws(perm_line), trimws(perm_mean_line)))
  }
}
sink()
cat(sprintf("\nSaved: %s/comparison_summary.txt\n", out_dir))
