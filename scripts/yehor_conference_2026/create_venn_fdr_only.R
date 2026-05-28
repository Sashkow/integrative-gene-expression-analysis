#!/usr/bin/env Rscript
library(ggplot2)

BASE <- "articles/yehor_conference_2026/data/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko"
SEX_DIR <- file.path(BASE, "sex_stratified")
OUTPUT <- "articles/yehor_conference_2026/data/venn_m_vs_f_fdr_only.png"

de_m <- read.delim(file.path(SEX_DIR, "difexp_males_1t_vs_2t.tsv"))
de_f <- read.delim(file.path(SEX_DIR, "difexp_females_1t_vs_2t.tsv"))

sig_m <- de_m$gene[de_m$adj.P.Val < 0.05]
sig_f <- de_f$gene[de_f$adj.P.Val < 0.05]

both <- length(intersect(sig_m, sig_f))
m_only <- length(setdiff(sig_m, sig_f))
f_only <- length(setdiff(sig_f, sig_m))

cat(sprintf("Males FDR<0.05: %d\n", length(sig_m)))
cat(sprintf("Females FDR<0.05: %d\n", length(sig_f)))
cat(sprintf("Shared: %d\n", both))
cat(sprintf("Males only: %d\n", m_only))
cat(sprintf("Females only: %d\n", f_only))

if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram", repos = "https://cloud.r-project.org")
}
library(VennDiagram)

png(OUTPUT, width = 8, height = 7, units = "in",
    res = 150, type = "cairo")
grid.newpage()
venn <- draw.pairwise.venn(
  area1 = length(sig_m),
  area2 = length(sig_f),
  cross.area = both,
  category = c(
    sprintf("Males 1T→2T\n(FDR<0.05, n=%d)", length(sig_m)),
    sprintf("Females 1T→2T\n(FDR<0.05, n=%d)", length(sig_f))
  ),
  fill = c("#1976D2", "#F57C00"),
  alpha = 0.4,
  cat.cex = 1.2,
  cex = 1.5,
  fontfamily = "sans",
  cat.fontfamily = "sans",
  cat.pos = c(-20, 20),
  cat.dist = 0.10,
  margin = 0.05
)
grid.text("DEGs: Males vs Females 1T→2T (FDR<0.05, no logFC filter)",
          y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()

cat("Saved:", OUTPUT, "\n")
