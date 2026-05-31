options(bitmapType = "cairo")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)

output_dir <- "output/article_validation/deg_enrichment"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

base_dir <- "output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko"

# --- Load DEG lists ---

sig_softimpute <- read.delim(file.path(base_dir, "difexp_significant_softimpute_combat_ref.tsv"),
                             stringsAsFactors = FALSE)
sig_intersection <- read.delim(file.path(base_dir, "difexp_significant_none_combat_ref.tsv"),
                               stringsAsFactors = FALSE)

full_genes <- as.character(sig_softimpute$gene)           # 380
intersection_genes <- as.character(sig_intersection$gene)  # 182
gained_genes <- setdiff(full_genes, intersection_genes)    # ~201

cat("Full DEGs:", length(full_genes), "\n")
cat("Intersection DEGs:", length(intersection_genes), "\n")
cat("Gained DEGs:", length(gained_genes), "\n\n")

# --- Background: all tested genes from the softImpute run ---

full_limma <- read.delim(file.path(base_dir, "difexp_softimpute_combat_ref.tsv"),
                         stringsAsFactors = FALSE)
universe <- as.character(full_limma$gene)
cat("Background universe:", length(universe), "genes\n\n")

# --- GO enrichment (Biological Process) ---

run_enrichment <- function(gene_list, name, universe) {
  cat("Running GO BP for:", name, "(", length(gene_list), "genes)\n")
  ego <- enrichGO(gene          = gene_list,
                  universe      = universe,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  readable      = TRUE)
  cat("  Significant GO BP terms:", nrow(as.data.frame(ego)), "\n")

  cat("Running KEGG for:", name, "\n")
  ekegg <- enrichKEGG(gene         = gene_list,
                      universe     = universe,
                      organism     = "hsa",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)
  cat("  Significant KEGG pathways:", nrow(as.data.frame(ekegg)), "\n\n")

  list(go = ego, kegg = ekegg)
}

enrich_full <- run_enrichment(full_genes, "full_380", universe)
enrich_intersection <- run_enrichment(intersection_genes, "intersection_182", universe)
enrich_gained <- run_enrichment(gained_genes, "gained_201", universe)

# --- Save enrichment tables ---

save_enrichment <- function(result, prefix) {
  go_df <- as.data.frame(result$go)
  kegg_df <- as.data.frame(result$kegg)
  if (nrow(go_df) > 0)
    write.csv(go_df, file.path(output_dir, paste0(prefix, "_GO_BP.csv")), row.names = FALSE)
  if (nrow(kegg_df) > 0)
    write.csv(kegg_df, file.path(output_dir, paste0(prefix, "_KEGG.csv")), row.names = FALSE)
}

save_enrichment(enrich_full, "enrichment_full_380")
save_enrichment(enrich_intersection, "enrichment_intersection_182")
save_enrichment(enrich_gained, "enrichment_gained_201")

# --- compareCluster: side-by-side enrichment ---

gene_clusters <- list(
  "Full (380)" = full_genes,
  "Intersection (182)" = intersection_genes,
  "Gained (~201)" = gained_genes
)

cat("Running compareCluster GO BP...\n")
cc_go <- compareCluster(geneCluster = gene_clusters,
                        fun = "enrichGO",
                        OrgDb = org.Hs.eg.db,
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        universe = universe,
                        readable = TRUE)

cat("Running compareCluster KEGG...\n")
cc_kegg <- compareCluster(geneCluster = gene_clusters,
                          fun = "enrichKEGG",
                          organism = "hsa",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          universe = universe)

# --- Dot plots ---

if (nrow(as.data.frame(cc_go)) > 0) {
  p_go <- dotplot(cc_go, showCategory = 15) +
    ggtitle("GO Biological Process: DEG set comparison") +
    theme(axis.text.y = element_text(size = 8))

  ggsave(file.path(output_dir, "fig_compareCluster_GO_BP.pdf"), p_go,
         width = 12, height = 10)
  png(file.path(output_dir, "fig_compareCluster_GO_BP.png"),
      width = 3600, height = 3000, res = 300, type = "cairo")
  print(p_go)
  invisible(dev.off())
  cat("Saved GO BP compareCluster dot plot\n")
}

if (nrow(as.data.frame(cc_kegg)) > 0) {
  p_kegg <- dotplot(cc_kegg, showCategory = 15) +
    ggtitle("KEGG pathways: DEG set comparison") +
    theme(axis.text.y = element_text(size = 8))

  ggsave(file.path(output_dir, "fig_compareCluster_KEGG.pdf"), p_kegg,
         width = 12, height = 8)
  png(file.path(output_dir, "fig_compareCluster_KEGG.png"),
      width = 3600, height = 2400, res = 300, type = "cairo")
  print(p_kegg)
  invisible(dev.off())
  cat("Saved KEGG compareCluster dot plot\n")
}

# --- Individual dot plots for each set ---

plot_individual <- function(result, prefix, title_prefix) {
  go_df <- as.data.frame(result$go)
  if (nrow(go_df) > 0) {
    p <- dotplot(result$go, showCategory = 20) +
      ggtitle(paste0(title_prefix, ": GO Biological Process"))
    ggsave(file.path(output_dir, paste0(prefix, "_GO_BP_dotplot.pdf")), p,
           width = 10, height = 8)
  }
  kegg_df <- as.data.frame(result$kegg)
  if (nrow(kegg_df) > 0) {
    p <- dotplot(result$kegg, showCategory = 20) +
      ggtitle(paste0(title_prefix, ": KEGG pathways"))
    ggsave(file.path(output_dir, paste0(prefix, "_KEGG_dotplot.pdf")), p,
           width = 10, height = 8)
  }
}

plot_individual(enrich_full, "full_380", "Full DEGs (380)")
plot_individual(enrich_intersection, "intersection_182", "Intersection DEGs (182)")
plot_individual(enrich_gained, "gained_201", "Gained DEGs (~201)")

# --- Summary ---

cat("\n=== Enrichment summary ===\n")
cat(sprintf("Full 380: %d GO BP terms, %d KEGG pathways\n",
            nrow(as.data.frame(enrich_full$go)), nrow(as.data.frame(enrich_full$kegg))))
cat(sprintf("Intersection 182: %d GO BP terms, %d KEGG pathways\n",
            nrow(as.data.frame(enrich_intersection$go)), nrow(as.data.frame(enrich_intersection$kegg))))
cat(sprintf("Gained ~201: %d GO BP terms, %d KEGG pathways\n",
            nrow(as.data.frame(enrich_gained$go)), nrow(as.data.frame(enrich_gained$kegg))))

cat("\nAll outputs saved to:", output_dir, "\n")
