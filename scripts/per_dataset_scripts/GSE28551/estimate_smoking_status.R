suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(cluster)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

# Smoking-responsive genes grouped by functional category:
#
# --- AHR pathway / xenobiotic metabolism ---
# CYP1B1, AHR, ARNT, TIPARP, EPHX1, NQO1, AKR1B10
# Cigarette smoke PAHs activate AHR, which dimerizes with ARNT and induces
# CYP1B1/TIPARP transcription. EPHX1 is phase I epoxide hydrolase; NQO1 is
# phase II quinone reductase (also Nrf2-regulated); AKR1B10 is an aldo-keto
# reductase induced via Nrf2 oxidative stress response.
#   - Huuskonen et al. 2008 (PMID 17928820): AhR-CYP1A1 is the dominant
#     induced pathway in placentas of smoking mothers.
#   - Penning 2017 (doi:10.1021/acs.chemrestox.6b00319): NQO1, AKR1B10,
#     CYP1B1 regulation by Nrf2/AhR in stress and xenobiotic response.
#
# --- Whole-blood smoking signature ---
# LRRN3, SASH1
# Reproducibly the top DEGs in blood of smokers vs non-smokers across studies.
#   - Huan et al. 2016 (Hum Mol Genet 25:4611): LRRN3 is the #1 smoking gene
#     in a meta-analysis of >10,000 blood transcriptomes.
#   - Beineke et al. 2012 (BMC Med Genomics 5:58): LRRN3 and SASH1 form the
#     core of a whole-blood smoking status signature.
#   - Poussin et al. 2017 (PMID 28085253): LRRN3, SASH1 confirmed across
#     independent teams in the sbv IMPROVER challenge.
#
# --- Inflammatory mediators / oxidative stress ---
# S100A8, S100A9 (calprotectin subunits), CXCL1, CXCL2, IL1B, TGFB1
# S100A8/A9 are DAMPs released by neutrophils under oxidative stress;
# CXCL1/CXCL2 are neutrophil-recruiting chemokines; IL1B is a central
# pro-inflammatory cytokine; TGFB1 mediates tissue remodeling and
# immunosuppression. All elevated in smoking-related inflammation.
#   - Dehghani et al. 2025 (Front Immunol 16:1590290): smoke exposure in
#     pregnancy exacerbates CXCL1/CXCL2-driven neutrophil inflammation
#     and alters placental immune landscape.
#
# --- Placental markers ---
# LEP, CGA, HSD11B2
# LEP (leptin): produced by syncytiotrophoblast, sensitive to cadmium from
# cigarette smoke.
#   - Stasenko et al. 2010 (PMID 19847775): cadmium reduces placental LEP
#     mRNA dose-dependently in trophoblast culture.
# CGA: alpha subunit of hCG, trophoblast-specific hormone marker.
# HSD11B2: 11beta-hydroxysteroid dehydrogenase type 2, converts cortisol to
# cortisone in placenta, protecting the fetus. Smoking impairs this barrier.
#   - Appleton et al. 2013 (PMID 24040322): placental HSD11B2 methylation
#     is susceptible to maternal environmental stressors.
smoking_symbols <- c(
  "CYP1B1", "NQO1", "AKR1B10", "TIPARP", "LRRN3", "SASH1",
  "ARNT", "AHR", "EPHX1", "S100A8", "S100A9",
  "CXCL1", "CXCL2", "IL1B", "TGFB1", "LEP", "CGA", "HSD11B2"
)

mapped <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = smoking_symbols,
  columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL"
))
mapped <- mapped[!is.na(mapped$ENTREZID) & !duplicated(mapped$SYMBOL), ]
sym_to_eid <- setNames(mapped$ENTREZID, mapped$SYMBOL)

d <- read.delim("data/mapped/GSE28551.tsv", row.names = 1)

pdata <- read.csv("data/phenodata/samples.csv", stringsAsFactors = FALSE)
pdata <- pdata[pdata$secondaryaccession == "GSE28551", ]
ft_samples <- pdata$arraydatafile_exprscolumnnames[
  pdata$Gestational.Age.Category == "First Trimester"
]
smoking_eids <- unname(sym_to_eid[sym_to_eid %in% rownames(d)])
ft <- d[smoking_eids, ft_samples]
ft_z <- t(scale(t(ft)))

cat("Samples:", ncol(ft), "\n")
cat("Smoking genes found:", nrow(ft), "of", length(smoking_symbols), "\n\n")

# Clustering: Ward.D2 and k-means both give 7 vs 9
hc <- hclust(dist(t(ft_z)), method = "ward.D2")
groups_ward <- cutree(hc, k = 2)

set.seed(42)
km <- kmeans(t(ft_z), centers = 2, nstart = 25)

cat("Ward.D2 k=2:", table(groups_ward), "\n")
cat("K-means k=2:", table(km$cluster), "\n")
cat("Agreement:", max(sum(groups_ward == km$cluster),
                      sum(groups_ward != km$cluster)), "/ 16\n\n")

# Identify which cluster is smokers by CYP1B1/NQO1/AKR1B10 expression
# (higher in smokers due to AHR pathway activation)
marker_eids <- sym_to_eid[c("CYP1B1", "NQO1", "AKR1B10")]
m1 <- mean(as.numeric(unlist(ft[marker_eids, groups_ward == 1])))
m2 <- mean(as.numeric(unlist(ft[marker_eids, groups_ward == 2])))

if (m1 > m2) {
  smoker_cluster <- 1
} else {
  smoker_cluster <- 2
}

smoker_ids <- names(groups_ward[groups_ward == smoker_cluster])
nonsmoker_ids <- names(groups_ward[groups_ward != smoker_cluster])

cat("Smoker cluster (higher CYP1B1/NQO1/AKR1B10):", smoker_cluster, "\n")
cat("Smokers (", length(smoker_ids), "):", paste(smoker_ids, collapse = ", "), "\n")
cat("Non-smokers (", length(nonsmoker_ids), "):", paste(nonsmoker_ids, collapse = ", "), "\n\n")

# Sitras et al. 2012 (PMID 22442682) reports 9/16 smokers in first trimester
# The 9/7 split matches exactly.
cat("Article reports: 9 smokers, 7 non-smokers in first trimester group\n")
cat("Our clustering:  ", length(smoker_ids), "vs", length(nonsmoker_ids), "\n\n")

# Silhouette width
dm <- dist(t(ft_z))
pm <- pam(dm, k = 2)
cat("PAM silhouette (k=2):", round(pm$silinfo$avg.width, 3), "\n\n")

# Per-gene means
cat("Mean expression per cluster:\n")
cat(sprintf("%-10s  Non-smoker  Smoker   Diff\n", "Gene"))
for (sym in smoking_symbols) {
  eid <- sym_to_eid[sym]
  if (!eid %in% rownames(ft)) next
  m_ns <- mean(as.numeric(ft[eid, nonsmoker_ids]))
  m_sm <- mean(as.numeric(ft[eid, smoker_ids]))
  cat(sprintf("%-10s  %9.2f  %7.2f  %+.2f\n", sym, m_ns, m_sm, m_sm - m_ns))
}

# Validation: re-cluster just the 7 non-smokers — should not split meaningfully
cat("\n=== Validation: clustering 7 non-smokers (should be homogeneous) ===\n")
ft_ns <- ft[, nonsmoker_ids]
ft_ns_z <- t(scale(t(ft_ns)))
hc_ns <- hclust(dist(t(ft_ns_z)), method = "ward.D2")
g_ns <- cutree(hc_ns, k = 2)
cat("Ward k=2 on non-smokers:", table(g_ns), "(expect 1 outlier, not bimodal)\n")
pm_ns <- pam(dist(t(ft_ns_z)), k = 2)
cat("PAM silhouette:", round(pm_ns$silinfo$avg.width, 3),
    "(expect much lower than", round(pm$silinfo$avg.width, 3), ")\n")

# Plot
png("output/GSE28551_smoking_estimation.png", width = 1400, height = 600, res = 120)
par(mfrow = c(1, 3))

# Dendrogram
plot(hc, main = "GSE28551 1st trimester (n=16)\n18 smoking-responsive genes",
     xlab = "", sub = "Ward.D2", cex = 0.7)
rect.hclust(hc, k = 2, border = c("red", "blue"))

# PCA
pc <- prcomp(t(ft_z))
cols <- ifelse(groups_ward == smoker_cluster, "blue", "red")
plot(pc$x[, 1], pc$x[, 2], col = cols, pch = 19, cex = 1.5,
     xlab = paste0("PC1 (", round(100 * summary(pc)$importance[2, 1], 1), "%)"),
     ylab = paste0("PC2 (", round(100 * summary(pc)$importance[2, 2], 1), "%)"),
     main = "PCA — smoking gene signature")
text(pc$x[, 1], pc$x[, 2], labels = gsub("GSM7070", "", rownames(pc$x)),
     pos = 3, cex = 0.7)
legend("topright",
       legend = c(paste0("Smoker (n=", length(smoker_ids), ")"),
                  paste0("Non-smoker (n=", length(nonsmoker_ids), ")")),
       col = c("blue", "red"), pch = 19)

# Heatmap of key markers
heatmap(ft_z, Colv = as.dendrogram(hc), scale = "none",
        col = colorRampPalette(c("navy", "white", "firebrick"))(50),
        labRow = smoking_symbols[smoking_symbols %in% mapped$SYMBOL[match(rownames(ft_z), mapped$ENTREZID)]],
        main = "Smoking genes heatmap", margins = c(8, 8))

dev.off()
cat("\nSaved: output/GSE28551_smoking_estimation.png\n")
