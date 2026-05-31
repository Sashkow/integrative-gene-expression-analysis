#!/usr/bin/env Rscript
#
# Apply the restoration STRING + PPI + clustering pipeline to phase2b output
# using the newest STRINGdb version (v12.0) instead of the v11.0 the
# restoration uses. Then compare to restoration's difexp_final.csv.
#
# Mirrors these steps from
# scripts/article_4_trims_compare_restoration/run.R:pipeline_restored():
#   - filter DE by |logFC| >= 1 & adj.P.Val < 0.05
#   - STRINGdb$map on background (expression matrix) to set background
#   - STRINGdb$map on DE genes (takeFirst = FALSE), dedupe STRING_ids
#   - get_subnetwork -> simplify -> remove isolated nodes
#   - fastgreedy.community clustering
#
# Usage:
#   Rscript scripts/integrative_analysis/phase2b_direct_merge/validate_vs_restoration.R \
#     [method]     # default: softimpute_combat

suppressPackageStartupMessages({
  library(STRINGdb)
  library(igraph)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
method <- if (length(args) >= 1) args[1] else "softimpute_combat"

phase2b_dir <- "output/phase2b_combat/phase2b_1_2_restoration"
restoration_dir <- "scripts/article_4_trims_compare_restoration"
out_dir <- file.path(phase2b_dir, paste0("validation_", method))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n=== Validating", method, "with STRINGdb v12.0 ===\n\n")

# ---- 1. Load phase2b outputs --------------------------------------------

de_all <- read.delim(file.path(phase2b_dir, paste0("difexp_", method, ".tsv")),
                     stringsAsFactors = FALSE)
exprs <- read.delim(file.path(phase2b_dir, paste0("exprs_", method, ".tsv")),
                    row.names = 1, check.names = FALSE)
de_all$gene <- as.character(de_all$gene)
cat("phase2b DE (all):", nrow(de_all), "genes\n")
cat("phase2b background:", nrow(exprs), "genes x", ncol(exprs), "samples\n")

# ---- 2. ENTREZID -> SYMBOL mapping --------------------------------------

all_entrez <- unique(c(de_all$gene, rownames(exprs)))
sym_map <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db, keys = all_entrez,
                         columns = c("ENTREZID", "SYMBOL", "GENENAME"),
                         keytype = "ENTREZID")
)
sym_map <- sym_map[!is.na(sym_map$SYMBOL), ]
sym_map <- sym_map[!duplicated(sym_map$ENTREZID), ]
cat("ENTREZID -> SYMBOL:", nrow(sym_map), "/", length(all_entrez), "mapped\n")

de_all <- merge(de_all, sym_map, by.x = "gene", by.y = "ENTREZID")
colnames(de_all)[colnames(de_all) == "gene"] <- "ENTREZID"
de_all <- de_all[!is.na(de_all$SYMBOL) & de_all$SYMBOL != "", ]
de_all <- de_all[!duplicated(de_all$SYMBOL), ]

exprs$ENTREZID <- rownames(exprs)
exprs <- merge(exprs, sym_map[, c("ENTREZID", "SYMBOL")], by = "ENTREZID")
exprs <- exprs[!duplicated(exprs$SYMBOL), ]
rownames(exprs) <- exprs$SYMBOL

# ---- 3. Filter DE by |logFC| >= 1 & FDR < 0.05 --------------------------

de_filt <- de_all[abs(de_all$logFC) >= 1 & de_all$adj.P.Val < 0.05, ]
cat("DE with |logFC|>=1 & FDR<0.05:", nrow(de_filt), "\n\n")

# ---- 4. STRINGdb v12 setup ----------------------------------------------

cache_dir <- file.path(phase2b_dir, "stringdb_cache_v12")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
string_db <- STRINGdb$new(
  version = "12.0", species = 9606,
  score_threshold = 100, input_directory = cache_dir
)
cat("STRINGdb version: 12.0\n")

# ---- 5. Map background (set STRING background) -------------------------

bg_mapped <- string_db$map(exprs, "SYMBOL", removeUnmappedRows = TRUE)
backgroundV <- bg_mapped$STRING_id
string_db$set_background(backgroundV)
cat("Background mapped to STRING:", nrow(bg_mapped), "/", nrow(exprs), "\n")

# ---- 6. Map DE genes ----------------------------------------------------

de_mapped <- string_db$map(
  de_filt, "SYMBOL", removeUnmappedRows = TRUE, takeFirst = FALSE
)
de_mapped <- de_mapped[!duplicated(de_mapped$STRING_id), ]
cat("DE mapped to STRING:", nrow(de_mapped), "/", nrow(de_filt), "\n")

# ---- 7. PPI subnetwork, remove isolated nodes ---------------------------

G <- igraph::simplify(string_db$get_subnetwork(de_mapped$STRING_id))
iso <- which(igraph::degree(G) == 0)
G <- igraph::delete_vertices(G, iso)
cat("PPI subnetwork nodes:", length(igraph::V(G)),
    " edges:", length(igraph::E(G)), "\n")

de_mapped$on_graph <- ifelse(
  de_mapped$STRING_id %in% names(igraph::V(G)), 1, 0
)
de_mapped <- de_mapped[de_mapped$on_graph == 1, ]

# ---- 8. Fastgreedy clustering -------------------------------------------

fg <- igraph::cluster_fast_greedy(G, merges = TRUE, modularity = TRUE)
cat("Modularity:", igraph::modularity(fg), "\n")
de_mapped$cluster <- igraph::membership(fg)[de_mapped$STRING_id]
de_mapped$updown <- ifelse(de_mapped$logFC > 0, "green", "red")
cat("Final genes post-STRING pipeline:", nrow(de_mapped), "\n\n")

# ---- 9. Save --------------------------------------------------------

final_file <- file.path(out_dir, "difexp_final.csv")
write.table(de_mapped, final_file, row.names = FALSE, sep = ",")
cat("Saved:", final_file, "\n\n")

# ---- 10. Compare to restoration reference ------------------------------

ref <- read.csv(file.path(restoration_dir, "original_results", "1_2",
                          "difexp_final.csv"), stringsAsFactors = FALSE)
cat("=== Comparison to restoration reference ===\n")
cat("Reference:", nrow(ref), "rows,",
    length(unique(ref$SYMBOL)), "unique SYMBOLs\n")
cat("Phase2b (post-STRING v12):", nrow(de_mapped), "rows,",
    length(unique(de_mapped$SYMBOL)), "unique SYMBOLs\n\n")

ref_sym <- unique(ref$SYMBOL)
new_sym <- unique(de_mapped$SYMBOL)
common <- intersect(ref_sym, new_sym)
cat("Overlap (SYMBOL):\n")
cat("  Common:", length(common), "/", length(ref_sym), "ref =",
    round(100 * length(common) / length(ref_sym), 1), "%\n")
cat("  New only:", length(setdiff(new_sym, ref_sym)), "\n")
cat("  Ref only:", length(setdiff(ref_sym, new_sym)), "\n\n")

ref_cmp <- ref[!duplicated(ref$SYMBOL), c("SYMBOL", "logFC")]
new_cmp <- de_mapped[!duplicated(de_mapped$SYMBOL), c("SYMBOL", "logFC")]
colnames(ref_cmp)[2] <- "logFC_ref"
colnames(new_cmp)[2] <- "logFC_new"
m <- merge(ref_cmp, new_cmp, by = "SYMBOL")
cat("logFC on", nrow(m), "common genes:\n")
cat("  Pearson r:", round(cor(m$logFC_ref, m$logFC_new), 4), "\n")
cat("  Spearman r:",
    round(cor(m$logFC_ref, m$logFC_new, method = "spearman"), 4), "\n")
cat("  Direction concordance:",
    round(mean(sign(m$logFC_ref) == sign(m$logFC_new)) * 100, 1), "%\n\n")

write.csv(m, file.path(out_dir, "common_genes_logfc.csv"), row.names = FALSE)
writeLines(setdiff(new_sym, ref_sym),
           file.path(out_dir, "new_only_symbols.txt"))
writeLines(setdiff(ref_sym, new_sym),
           file.path(out_dir, "ref_only_symbols.txt"))

cat("Done. Output:", out_dir, "\n")
