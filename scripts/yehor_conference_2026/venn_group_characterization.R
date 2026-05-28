#!/usr/bin/env Rscript
#
# Characterization of Venn diagram groups from sex-stratified DEG analysis.
# Annotates genes, runs GO/Hallmark enrichment, produces interpretation.
#
# Usage:
#   Rscript scripts/yehor_conference_2026/venn_group_characterization.R

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(DOSE)
  library(msigdbr)
  library(ggplot2)
  library(dplyr)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

base_dir <- paste0(
  "articles/yehor_conference_2026/data/",
  "phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko"
)
sex_dir <- file.path(base_dir, "sex_stratified")
output_dir <- file.path(sex_dir, "venn_analysis")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

fdr_thresh <- 0.05
logfc_thresh <- 1.0

# ── 1. Load DEG lists and compute Venn groups ──

cat("=== Loading DE tables ===\n\n")

sig_all <- read.delim(
  file.path(base_dir, "difexp_significant_softimpute_combat_ref.tsv"),
  stringsAsFactors = FALSE
)
sig_m <- read.delim(
  file.path(sex_dir, "difexp_significant_males_1t_vs_2t.tsv"),
  stringsAsFactors = FALSE
)
sig_f <- read.delim(
  file.path(sex_dir, "difexp_significant_females_1t_vs_2t.tsv"),
  stringsAsFactors = FALSE
)

de_all <- read.delim(
  file.path(base_dir, "difexp_softimpute_combat_ref.tsv"),
  stringsAsFactors = FALSE
)
de_m <- read.delim(
  file.path(sex_dir, "difexp_males_1t_vs_2t.tsv"),
  stringsAsFactors = FALSE
)
de_f <- read.delim(
  file.path(sex_dir, "difexp_females_1t_vs_2t.tsv"),
  stringsAsFactors = FALSE
)

genes_all <- as.character(sig_all$gene)
genes_m <- as.character(sig_m$gene)
genes_f <- as.character(sig_f$gene)

cat(sprintf("  1_2 (all):     %d sig DEGs\n", length(genes_all)))
cat(sprintf("  1_2_m (males): %d sig DEGs\n", length(genes_m)))
cat(sprintf("  1_2_f (females): %d sig DEGs\n\n", length(genes_f)))

shared_all_3 <- intersect(intersect(genes_all, genes_m), genes_f)
male_unique <- setdiff(genes_m, union(genes_all, genes_f))
female_unique <- setdiff(genes_f, union(genes_all, genes_m))
combined_unique <- setdiff(genes_all, union(genes_m, genes_f))
male_combined <- setdiff(intersect(genes_m, genes_all), genes_f)
female_combined <- setdiff(intersect(genes_f, genes_all), genes_m)
male_female_only <- setdiff(intersect(genes_m, genes_f), genes_all)

venn_groups <- list(
  shared_all_3 = shared_all_3,
  male_unique = male_unique,
  female_unique = female_unique,
  combined_unique = combined_unique,
  male_combined = male_combined,
  female_combined = female_combined,
  male_female_only = male_female_only
)

cat("=== Venn group sizes ===\n")
for (nm in names(venn_groups)) {
  cat(sprintf("  %-20s %d\n", nm, length(venn_groups[[nm]])))
}
total <- sum(sapply(venn_groups, length))
cat(sprintf("  %-20s %d\n\n", "TOTAL", total))

# ── 2. Annotate genes ──

cat("=== Annotating genes ===\n")

all_genes <- unique(unlist(venn_groups))

ann <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = all_genes,
  columns = c("SYMBOL", "GENENAME", "CHR"),
  keytype = "ENTREZID"
)
ann <- ann[!duplicated(ann$ENTREZID), ]

chr_info <- tryCatch({
  cat("  Querying biomaRt for chromosome bands...\n")
  mart <- biomaRt::useMart("ensembl",
                            dataset = "hsapiens_gene_ensembl")
  bm <- biomaRt::getBM(
    filters = "entrezgene_id",
    values = all_genes,
    attributes = c("entrezgene_id", "chromosome_name", "band"),
    mart = mart
  )
  bm <- bm[!duplicated(bm$entrezgene_id), ]
  bm$entrezgene_id <- as.character(bm$entrezgene_id)
  cat(sprintf("  biomaRt returned %d annotations\n", nrow(bm)))
  bm
}, error = function(e) {
  cat("  biomaRt failed, using org.Hs.eg.db CHR fallback\n")
  NULL
})

get_chr <- function(entrezid) {
  if (!is.null(chr_info)) {
    idx <- match(entrezid, chr_info$entrezgene_id)
    if (!is.na(idx)) return(chr_info$chromosome_name[idx])
  }
  idx <- match(entrezid, ann$ENTREZID)
  if (!is.na(idx)) return(ann$CHR[idx])
  NA_character_
}

get_band <- function(entrezid) {
  if (!is.null(chr_info)) {
    idx <- match(entrezid, chr_info$entrezgene_id)
    if (!is.na(idx)) return(chr_info$band[idx])
  }
  NA_character_
}

build_master <- function(gene_ids, group_name) {
  df <- data.frame(ENTREZID = gene_ids, stringsAsFactors = FALSE)
  df$venn_group <- group_name

  idx_ann <- match(df$ENTREZID, ann$ENTREZID)
  df$SYMBOL <- ann$SYMBOL[idx_ann]
  df$GENENAME <- ann$GENENAME[idx_ann]

  df$chromosome <- sapply(df$ENTREZID, get_chr)
  df$band <- sapply(df$ENTREZID, get_band)
  df$is_sex_chr <- df$chromosome %in% c("X", "Y")

  idx_all <- match(df$ENTREZID, as.character(de_all$gene))
  df$logFC_all <- de_all$logFC[idx_all]
  df$adjP_all <- de_all$adj.P.Val[idx_all]
  df$sig_all <- df$ENTREZID %in% genes_all

  idx_m <- match(df$ENTREZID, as.character(de_m$gene))
  df$logFC_m <- de_m$logFC[idx_m]
  df$adjP_m <- de_m$adj.P.Val[idx_m]
  df$sig_m <- df$ENTREZID %in% genes_m

  idx_f <- match(df$ENTREZID, as.character(de_f$gene))
  df$logFC_f <- de_f$logFC[idx_f]
  df$adjP_f <- de_f$adj.P.Val[idx_f]
  df$sig_f <- df$ENTREZID %in% genes_f

  df$direction_concordant <- sign(df$logFC_m) == sign(df$logFC_f)

  df$near_miss_m <- !is.na(df$adjP_m) &
    df$adjP_m < 0.10 & abs(df$logFC_m) > 0.8 & !df$sig_m
  df$near_miss_f <- !is.na(df$adjP_f) &
    df$adjP_f < 0.10 & abs(df$logFC_f) > 0.8 & !df$sig_f
  df$near_miss_all <- !is.na(df$adjP_all) &
    df$adjP_all < 0.10 & abs(df$logFC_all) > 0.8 & !df$sig_all

  df
}

master <- do.call(rbind, lapply(names(venn_groups), function(nm) {
  if (length(venn_groups[[nm]]) == 0) return(NULL)
  build_master(venn_groups[[nm]], nm)
}))
rownames(master) <- NULL

cat(sprintf("  Master table: %d genes\n\n", nrow(master)))

# ── 3. Save annotation tables ──

cat("=== Saving annotation tables ===\n")

write.table(master,
            file.path(output_dir, "all_groups_annotated.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

for (nm in names(venn_groups)) {
  sub <- master[master$venn_group == nm, ]
  if (nrow(sub) == 0) next
  write.table(
    sub,
    file.path(output_dir, paste0("group_", nm, "_annotated.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  cat(sprintf("  group_%s_annotated.tsv (%d genes)\n", nm, nrow(sub)))
}
cat("\n")

# ── 4. GO enrichment ──

cat("=== GO enrichment (BP) ===\n\n")

universe <- as.character(de_all$gene)

run_go <- function(gene_ids, label) {
  if (length(gene_ids) < 5) {
    cat(sprintf("  %s: %d genes, skipping\n", label, length(gene_ids)))
    return(NULL)
  }
  cat(sprintf("  %s: %d genes ... ", label, length(gene_ids)))
  ego <- enrichGO(
    gene = gene_ids,
    universe = universe,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    keyType = "ENTREZID",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
  )
  n_sig <- sum(ego@result$p.adjust < 0.05)
  cat(sprintf("%d significant terms\n", n_sig))

  if (n_sig > 0) {
    write.csv(
      ego@result[ego@result$p.adjust < 0.05, ],
      file.path(output_dir, paste0("go_bp_", label, ".csv")),
      row.names = FALSE
    )
    n_show <- min(20, n_sig)
    png(file.path(output_dir, paste0("go_bp_", label, "_dotplot.png")),
        width = 10, height = max(4, n_show * 0.4),
        units = "in", res = 200, type = "cairo")
    print(dotplot(ego, showCategory = n_show,
                  title = paste("GO BP:", label)))
    dev.off()
  }
  ego
}

ego_shared <- run_go(shared_all_3, "shared_all_3")
ego_female_unique <- run_go(female_unique, "female_unique")
ego_female_combined <- run_go(female_combined, "female_combined")
ego_male_unique <- run_go(male_unique, "male_unique")
ego_combined_unique <- run_go(combined_unique, "combined_unique")
ego_male_combined <- run_go(male_combined, "male_combined")
ego_male_female <- run_go(male_female_only, "male_female_only")

fu_sub <- master[master$venn_group == "female_unique", ]
fu_up <- fu_sub$ENTREZID[fu_sub$logFC_f > 0]
fu_down <- fu_sub$ENTREZID[fu_sub$logFC_f < 0]
ego_fu_up <- run_go(fu_up, "female_unique_up")
ego_fu_down <- run_go(fu_down, "female_unique_down")
cat("\n")

# compareCluster
cat("  compareCluster across main groups ... ")
cc_list <- list()
if (length(shared_all_3) >= 5) cc_list[["Shared all 3"]] <- shared_all_3
if (length(female_unique) >= 5) cc_list[["Female unique"]] <- female_unique
if (length(female_combined) >= 5) cc_list[["Female+Combined"]] <- female_combined
if (length(male_unique) >= 5) cc_list[["Male unique"]] <- male_unique

if (length(cc_list) >= 2) {
  cc <- compareCluster(
    geneClusters = cc_list,
    fun = "enrichGO",
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    keyType = "ENTREZID",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    readable = TRUE
  )
  n_cc <- nrow(cc@compareClusterResult[
    cc@compareClusterResult$p.adjust < 0.05, ])
  cat(sprintf("%d significant term-group pairs\n", n_cc))
  if (n_cc > 0) {
    png(file.path(output_dir, "go_bp_compare_groups_dotplot.png"),
        width = 14, height = 10, units = "in",
        res = 200, type = "cairo")
    print(dotplot(cc, showCategory = 8,
                  title = "GO BP: Venn group comparison"))
    dev.off()
  }
} else {
  cat("not enough groups with >= 5 genes\n")
}
cat("\n")

# ── 5. msigdbr Hallmark enrichment ──

cat("=== Hallmark enrichment ===\n\n")

h_sets <- msigdbr(species = "Homo sapiens", category = "H")
h_t2g <- h_sets %>%
  dplyr::distinct(gs_name, entrez_gene) %>%
  dplyr::mutate(entrez_gene = as.character(entrez_gene)) %>%
  as.data.frame()

run_hallmark <- function(gene_ids, label) {
  if (length(gene_ids) < 5) {
    cat(sprintf("  %s: %d genes, skipping\n", label, length(gene_ids)))
    return(NULL)
  }
  cat(sprintf("  %s: %d genes ... ", label, length(gene_ids)))
  enr <- enricher(
    gene = gene_ids,
    universe = universe,
    TERM2GENE = h_t2g,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  n_sig <- sum(enr@result$p.adjust < 0.05)
  cat(sprintf("%d significant\n", n_sig))
  if (n_sig > 0) {
    write.csv(
      enr@result[enr@result$p.adjust < 0.05, ],
      file.path(output_dir, paste0("hallmark_", label, ".csv")),
      row.names = FALSE
    )
  }
  enr
}

for (nm in names(venn_groups)) {
  run_hallmark(venn_groups[[nm]], nm)
}
cat("\n")

# ── 6. Plots ──

cat("=== Generating plots ===\n\n")

# logFC scatter: female_unique genes
fu <- master[master$venn_group == "female_unique", ]
if (nrow(fu) > 0) {
  p_scatter <- ggplot(fu, aes(x = logFC_m, y = logFC_f)) +
    geom_point(aes(color = near_miss_m), size = 2.5, alpha = 0.8) +
    scale_color_manual(
      values = c("FALSE" = "grey40", "TRUE" = "orange"),
      labels = c("Not near-miss", "Near-miss in males"),
      name = ""
    ) +
    geom_vline(xintercept = c(-1, 1),
               linetype = "dashed", color = "red", alpha = 0.5) +
    geom_hline(yintercept = c(-1, 1),
               linetype = "dashed", color = "red", alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dotted", color = "blue") +
    geom_text(
      data = fu[!is.na(fu$SYMBOL), ],
      aes(label = SYMBOL), size = 2.3,
      nudge_y = 0.15, check_overlap = TRUE
    ) +
    labs(
      x = "logFC (Males 1T vs 2T)",
      y = "logFC (Females 1T vs 2T)",
      title = "Female-unique DEGs: effect size comparison",
      subtitle = sprintf(
        "%d genes | %d concordant direction | %d near-miss in males",
        nrow(fu),
        sum(fu$direction_concordant, na.rm = TRUE),
        sum(fu$near_miss_m, na.rm = TRUE)
      )
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  png(file.path(output_dir, "female_unique_logfc_comparison.png"),
      width = 9, height = 8, units = "in", res = 150, type = "cairo")
  print(p_scatter)
  dev.off()
  cat("  female_unique_logfc_comparison.png\n")
}

# Chromosome distribution
chr_order <- c(as.character(1:22), "X", "Y")
master_chr <- master[master$chromosome %in% chr_order, ]
master_chr$chromosome <- factor(master_chr$chromosome, levels = chr_order)
master_chr$chr_type <- ifelse(
  master_chr$chromosome %in% c("X", "Y"), "Sex", "Autosome"
)

main_groups <- c("shared_all_3", "female_unique",
                 "female_combined", "male_unique")
master_chr_main <- master_chr[master_chr$venn_group %in% main_groups, ]
master_chr_main$venn_group <- factor(
  master_chr_main$venn_group,
  levels = main_groups,
  labels = c("Shared all 3", "Female unique",
             "Female+Combined", "Male unique")
)

if (nrow(master_chr_main) > 0) {
  p_chr <- ggplot(master_chr_main,
                  aes(x = chromosome, fill = chr_type)) +
    geom_bar() +
    facet_wrap(~ venn_group, scales = "free_y", ncol = 2) +
    scale_fill_manual(
      values = c("Autosome" = "grey60", "Sex" = "firebrick"),
      name = ""
    ) +
    labs(x = "Chromosome", y = "Count",
         title = "Chromosome distribution by Venn group") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      legend.position = "bottom"
    )

  png(file.path(output_dir, "chromosome_distribution_by_group.png"),
      width = 12, height = 7, units = "in", res = 150, type = "cairo")
  print(p_chr)
  dev.off()
  cat("  chromosome_distribution_by_group.png\n")
}
cat("\n")

# ── 7. Biological interpretation ──

cat("=== BIOLOGICAL INTERPRETATION ===\n\n")

interp <- character()
add <- function(...) {
  line <- sprintf(...)
  cat(line, "\n")
  interp <<- c(interp, line)
}

add("--- Male-unique DEGs (%d genes) ---", length(male_unique))
mu <- master[master$venn_group == "male_unique", ]
n_y <- sum(mu$chromosome == "Y", na.rm = TRUE)
n_x <- sum(mu$chromosome == "X", na.rm = TRUE)
add("  Y-linked: %d, X-linked: %d, autosomal: %d",
    n_y, n_x, nrow(mu) - n_y - n_x)
if (n_y > 0) {
  y_genes <- mu$SYMBOL[mu$chromosome == "Y" & !is.na(mu$chromosome)]
  add("  Y-linked genes: %s", paste(y_genes, collapse = ", "))
}
n_conc <- sum(mu$direction_concordant, na.rm = TRUE)
add("  Direction concordant with females: %d/%d", n_conc, nrow(mu))
n_nm_f <- sum(mu$near_miss_f, na.rm = TRUE)
add("  Near-miss in females: %d (power artifact candidates)", n_nm_f)
add("")

add("--- Female-unique DEGs (%d genes) ---", length(female_unique))
n_x_f <- sum(fu$chromosome == "X", na.rm = TRUE)
n_y_f <- sum(fu$chromosome == "Y", na.rm = TRUE)
add("  X-linked: %d, Y-linked: %d, autosomal: %d",
    n_x_f, n_y_f, nrow(fu) - n_x_f - n_y_f)
n_conc_f <- sum(fu$direction_concordant, na.rm = TRUE)
add("  Direction concordant with males: %d/%d (%.0f%%)",
    n_conc_f, nrow(fu),
    100 * n_conc_f / max(1, nrow(fu)))
n_nm_m <- sum(fu$near_miss_m, na.rm = TRUE)
add("  Near-miss in males: %d/%d (%.0f%%) -- likely power artifacts",
    n_nm_m, nrow(fu),
    100 * n_nm_m / max(1, nrow(fu)))
n_true_unique <- nrow(fu) - n_nm_m
add("  Likely truly sex-specific: %d genes", n_true_unique)
n_up_f <- sum(fu$logFC_f > 0)
n_down_f <- sum(fu$logFC_f < 0)
add("  Direction: %d up in 2T, %d down in 2T", n_up_f, n_down_f)
add("")

add("--- Combined-unique DEGs (%d genes) ---", length(combined_unique))
cu <- master[master$venn_group == "combined_unique", ]
n_conc_cu <- sum(cu$direction_concordant, na.rm = TRUE)
add("  Direction concordant M/F: %d/%d", n_conc_cu, nrow(cu))
add("  These genes reached significance only when both sexes")
add("  were pooled (increased power) with sex as covariate.")
if (nrow(cu) > 0 && any(!is.na(cu$SYMBOL))) {
  add("  Genes: %s", paste(na.omit(cu$SYMBOL), collapse = ", "))
}
add("")

add("--- Shared all 3 (%d genes) ---", length(shared_all_3))
add("  Core trimester-responsive genes, robust regardless of sex.")
sh <- master[master$venn_group == "shared_all_3", ]
n_x_sh <- sum(sh$chromosome == "X", na.rm = TRUE)
add("  X-linked among shared: %d", n_x_sh)
add("")

add("--- Female+Combined (%d genes) ---", length(female_combined))
fc <- master[master$venn_group == "female_combined", ]
n_nm_fc <- sum(fc$near_miss_m, na.rm = TRUE)
add("  Near-miss in males: %d/%d (%.0f%%)",
    n_nm_fc, nrow(fc),
    100 * n_nm_fc / max(1, nrow(fc)))
add("  These genes show trimester effect in the combined analysis")
add("  and in females, but not males -- likely a mix of power")
add("  effects (male 2T n=6) and genuinely stronger female response.")
add("")

add("--- Power context ---")
add("  Male samples:   46 (1T) + 6 (2T) = 52 total")
add("  Female samples: 56 (1T) + 9 (2T) = 65 total")
add("  The smaller male 2T group (n=6) reduces statistical power,")
add("  explaining why many genes are female-unique or female+combined")
add("  rather than genuinely sex-specific.")

writeLines(interp, file.path(output_dir, "interpretation.txt"))
cat(sprintf("\nInterpretation saved: %s/interpretation.txt\n", output_dir))

cat("\n=== Done ===\n")
