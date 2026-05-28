#!/usr/bin/env Rscript
#
# 2-way Venn (males vs females 1T→2T) with GO/Hallmark enrichment.
#
# Usage:
#   Rscript scripts/yehor_conference_2026/venn_m_vs_f_characterization.R

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(DOSE)
  library(msigdbr)
  library(VennDiagram)
  library(ggplot2)
  library(dplyr)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

base_dir <- paste0(
  "articles/yehor_conference_2026/data/",
  "phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko"
)
sex_dir <- file.path(base_dir, "sex_stratified")
output_dir <- file.path(sex_dir, "venn_m_vs_f")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

fdr_thresh <- 0.05
logfc_thresh <- 1.0

# ── 1. Load DEG lists ──

cat("=== Loading DE tables ===\n\n")

sig_m <- read.delim(
  file.path(sex_dir, "difexp_significant_males_1t_vs_2t.tsv"),
  stringsAsFactors = FALSE
)
sig_f <- read.delim(
  file.path(sex_dir, "difexp_significant_females_1t_vs_2t.tsv"),
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

genes_m <- as.character(sig_m$gene)
genes_f <- as.character(sig_f$gene)

cat(sprintf("  Males 1T vs 2T:   %d sig DEGs\n", length(genes_m)))
cat(sprintf("  Females 1T vs 2T: %d sig DEGs\n\n", length(genes_f)))

# ── 2. Compute Venn groups ──

shared <- intersect(genes_m, genes_f)
male_only <- setdiff(genes_m, genes_f)
female_only <- setdiff(genes_f, genes_m)

venn_groups <- list(
  shared = shared,
  male_only = male_only,
  female_only = female_only
)

cat("=== Venn group sizes ===\n")
for (nm in names(venn_groups)) {
  cat(sprintf("  %-15s %d\n", nm, length(venn_groups[[nm]])))
}
cat("\n")

# ── 3. Annotate genes ──

cat("=== Annotating genes ===\n")

all_genes <- unique(unlist(venn_groups))

ann <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = all_genes,
  columns = c("SYMBOL", "GENENAME", "CHR"),
  keytype = "ENTREZID"
)
ann <- ann[!duplicated(ann$ENTREZID), ]

build_master <- function(gene_ids, group_name) {
  df <- data.frame(ENTREZID = gene_ids, stringsAsFactors = FALSE)
  df$venn_group <- group_name

  idx_ann <- match(df$ENTREZID, ann$ENTREZID)
  df$SYMBOL <- ann$SYMBOL[idx_ann]
  df$GENENAME <- ann$GENENAME[idx_ann]
  df$chromosome <- ann$CHR[idx_ann]
  df$is_sex_chr <- df$chromosome %in% c("X", "Y")

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

  df
}

master <- do.call(rbind, lapply(names(venn_groups), function(nm) {
  if (length(venn_groups[[nm]]) == 0) return(NULL)
  build_master(venn_groups[[nm]], nm)
}))
rownames(master) <- NULL

cat(sprintf("  Master table: %d genes\n\n", nrow(master)))

# ── 4. Save annotation tables ──

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

# ── 5. Venn diagram with callouts ──

cat("=== Venn diagram ===\n")

mo <- master[master$venn_group == "male_only", ]
fo <- master[master$venn_group == "female_only", ]

mo_n_y <- sum(mo$chromosome == "Y", na.rm = TRUE)
mo_n_x <- sum(mo$chromosome == "X", na.rm = TRUE)
mo_conc <- sum(mo$direction_concordant, na.rm = TRUE)
mo_nm_f <- sum(mo$near_miss_f, na.rm = TRUE)

fo_n_x <- sum(fo$chromosome == "X", na.rm = TRUE)
fo_n_y <- sum(fo$chromosome == "Y", na.rm = TRUE)
fo_conc <- sum(fo$direction_concordant, na.rm = TRUE)
fo_nm_m <- sum(fo$near_miss_m, na.rm = TRUE)
fo_truly <- nrow(fo) - fo_nm_m

venn.plot <- venn.diagram(
  x = list(
    "Males 1T vs 2T" = genes_m,
    "Females 1T vs 2T" = genes_f
  ),
  filename = NULL,
  fill = c("#66C2A5", "#FC8D62"),
  alpha = 0.5,
  cex = 1.8,
  cat.cex = 1.3,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.06, 0.06),
  margin = 0.15,
  main = "DEG overlap: Males vs Females (1T vs 2T)",
  main.cex = 1.4
)
venn.plot <- venn.plot[!sapply(venn.plot, function(x) inherits(x, "rect"))]

venn_path <- file.path(output_dir, "venn_m_vs_f.png")
png(venn_path, width = 14, height = 9, units = "in",
    res = 150, type = "cairo")
grid::grid.newpage()

grid::pushViewport(grid::viewport(
  x = 0.5, y = 0.5, width = 0.50, height = 0.90
))
grid::grid.draw(venn.plot)
grid::popViewport()

dev.off()
cat(sprintf("  Saved: %s\n\n", venn_path))

# ── 6. GO enrichment ──

cat("=== GO enrichment (BP) ===\n\n")

universe <- as.character(de_m$gene)

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

ego_shared <- run_go(shared, "shared")
ego_male_only <- run_go(male_only, "male_only")
ego_female_only <- run_go(female_only, "female_only")

# Split by direction
mo_up <- mo$ENTREZID[mo$logFC_m > 0]
mo_down <- mo$ENTREZID[mo$logFC_m < 0]
fo_up <- fo$ENTREZID[fo$logFC_f > 0]
fo_down <- fo$ENTREZID[fo$logFC_f < 0]

ego_mo_up <- run_go(mo_up, "male_only_up")
ego_mo_down <- run_go(mo_down, "male_only_down")
ego_fo_up <- run_go(fo_up, "female_only_up")
ego_fo_down <- run_go(fo_down, "female_only_down")
cat("\n")

# compareCluster
cat("  compareCluster ... ")
cc_list <- list()
if (length(shared) >= 5) cc_list[["Shared"]] <- shared
if (length(male_only) >= 5) cc_list[["Male-only"]] <- male_only
if (length(female_only) >= 5) cc_list[["Female-only"]] <- female_only

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
    print(dotplot(cc, showCategory = 10,
                  title = "GO BP: Male-only vs Female-only vs Shared"))
    dev.off()
  }
} else {
  cat("not enough groups\n")
}
cat("\n")

# ── 7. Hallmark enrichment ──

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

# ── 8. Plots ──

cat("=== Generating plots ===\n\n")

# logFC scatter: female-only genes
if (nrow(fo) > 0) {
  p_scatter_f <- ggplot(fo, aes(x = logFC_m, y = logFC_f)) +
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
      data = fo[!is.na(fo$SYMBOL), ],
      aes(label = SYMBOL), size = 2.3,
      nudge_y = 0.15, check_overlap = TRUE
    ) +
    labs(
      x = "logFC (Males 1T vs 2T)",
      y = "logFC (Females 1T vs 2T)",
      title = "Female-only DEGs: effect size comparison",
      subtitle = sprintf(
        "%d genes | %d concordant | %d near-miss in males",
        nrow(fo), fo_conc, fo_nm_m
      )
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  png(file.path(output_dir, "female_only_logfc_comparison.png"),
      width = 9, height = 8, units = "in", res = 150, type = "cairo")
  print(p_scatter_f)
  dev.off()
  cat("  female_only_logfc_comparison.png\n")
}

# logFC scatter: male-only genes
if (nrow(mo) > 0) {
  p_scatter_m <- ggplot(mo, aes(x = logFC_m, y = logFC_f)) +
    geom_point(aes(color = near_miss_f), size = 2.5, alpha = 0.8) +
    scale_color_manual(
      values = c("FALSE" = "grey40", "TRUE" = "orange"),
      labels = c("Not near-miss", "Near-miss in females"),
      name = ""
    ) +
    geom_vline(xintercept = c(-1, 1),
               linetype = "dashed", color = "red", alpha = 0.5) +
    geom_hline(yintercept = c(-1, 1),
               linetype = "dashed", color = "red", alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dotted", color = "blue") +
    geom_text(
      data = mo[!is.na(mo$SYMBOL), ],
      aes(label = SYMBOL), size = 2.8,
      nudge_y = 0.1, check_overlap = TRUE
    ) +
    labs(
      x = "logFC (Males 1T vs 2T)",
      y = "logFC (Females 1T vs 2T)",
      title = "Male-only DEGs: effect size comparison",
      subtitle = sprintf(
        "%d genes | %d concordant | %d near-miss in females",
        nrow(mo), mo_conc, mo_nm_f
      )
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  png(file.path(output_dir, "male_only_logfc_comparison.png"),
      width = 9, height = 8, units = "in", res = 150, type = "cairo")
  print(p_scatter_m)
  dev.off()
  cat("  male_only_logfc_comparison.png\n")
}

# Chromosome distribution
chr_order <- c(as.character(1:22), "X", "Y")
master_chr <- master[master$chromosome %in% chr_order, ]
master_chr$chromosome <- factor(master_chr$chromosome, levels = chr_order)
master_chr$chr_type <- ifelse(
  master_chr$chromosome %in% c("X", "Y"), "Sex", "Autosome"
)
master_chr$venn_group <- factor(
  master_chr$venn_group,
  levels = c("shared", "male_only", "female_only"),
  labels = c("Shared", "Male-only", "Female-only")
)

if (nrow(master_chr) > 0) {
  p_chr <- ggplot(master_chr,
                  aes(x = chromosome, fill = chr_type)) +
    geom_bar() +
    facet_wrap(~ venn_group, scales = "free_y", ncol = 3) +
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

  png(file.path(output_dir, "chromosome_distribution.png"),
      width = 14, height = 5, units = "in", res = 150, type = "cairo")
  print(p_chr)
  dev.off()
  cat("  chromosome_distribution.png\n")
}
cat("\n")

# ── 9. Biological interpretation ──

cat("=== BIOLOGICAL INTERPRETATION ===\n\n")

interp <- character()
add <- function(...) {
  line <- sprintf(...)
  cat(line, "\n")
  interp <<- c(interp, line)
}

add("--- Male-only DEGs (%d genes) ---", length(male_only))
add("  Y-linked: %d, X-linked: %d, autosomal: %d",
    mo_n_y, mo_n_x, nrow(mo) - mo_n_y - mo_n_x)
add("  Direction concordant with females: %d/%d (%.0f%%)",
    mo_conc, nrow(mo), 100 * mo_conc / max(1, nrow(mo)))
add("  Near-miss in females: %d/%d (%.0f%%)",
    mo_nm_f, nrow(mo), 100 * mo_nm_f / max(1, nrow(mo)))
mo_truly <- nrow(mo) - mo_nm_f
add("  Likely truly sex-specific: ~%d genes", mo_truly)
add("")

add("--- Female-only DEGs (%d genes) ---", length(female_only))
add("  X-linked: %d, Y-linked: %d, autosomal: %d",
    fo_n_x, fo_n_y, nrow(fo) - fo_n_x - fo_n_y)
add("  Direction concordant with males: %d/%d (%.0f%%)",
    fo_conc, nrow(fo), 100 * fo_conc / max(1, nrow(fo)))
add("  Near-miss in males: %d/%d (%.0f%%)",
    fo_nm_m, nrow(fo), 100 * fo_nm_m / max(1, nrow(fo)))
add("  Likely truly sex-specific: ~%d genes", fo_truly)
fo_up_n <- sum(fo$logFC_f > 0)
fo_down_n <- sum(fo$logFC_f < 0)
add("  Direction: %d up in 2T, %d down in 2T", fo_up_n, fo_down_n)
add("")

add("--- Shared DEGs (%d genes) ---", length(shared))
sh <- master[master$venn_group == "shared", ]
sh_x <- sum(sh$chromosome == "X", na.rm = TRUE)
sh_y <- sum(sh$chromosome == "Y", na.rm = TRUE)
add("  X-linked: %d, Y-linked: %d", sh_x, sh_y)
sh_conc <- sum(sh$direction_concordant, na.rm = TRUE)
add("  Direction concordant: %d/%d (%.0f%%)",
    sh_conc, nrow(sh), 100 * sh_conc / max(1, nrow(sh)))
add("  Core trimester-responsive genes, robust in both sexes.")
add("")

add("--- Power context ---")
add("  Male samples:   46 (1T) + 6 (2T) = 52 total")
add("  Female samples: 56 (1T) + 9 (2T) = 65 total")
add("  Unbalanced 2T groups (6 vs 9) reduce male power,")
add("  contributing to the asymmetry (more female-only than male-only).")

writeLines(interp, file.path(output_dir, "interpretation.txt"))
cat(sprintf("\nInterpretation saved: %s/interpretation.txt\n", output_dir))

cat("\n=== Done ===\n")
