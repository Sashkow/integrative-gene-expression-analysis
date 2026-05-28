#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

out_dir <- "scripts/integrative_analysis/phase2b_direct_merge/config_cvs_vs_termination_compare"

# ── Eisenberg-Levanon housekeeping genes ──────────────────────────────────────

el_file <- "data/reference/integration_methods_references/ruv_housekeeping_genes_search/eisenberg_levanon_HK_genes.txt"
el_raw  <- read.delim(el_file, header = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)
el_symbols <- unique(trimws(el_raw$V1))
cat(sprintf("Eisenberg-Levanon list: %d symbols\n", length(el_symbols)))

el_map <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db, keys = el_symbols, columns = "ENTREZID", keytype = "SYMBOL")
)
el_map <- el_map[!is.na(el_map$ENTREZID), ]
el_map <- el_map[!duplicated(el_map$ENTREZID), ]
cat(sprintf("Mapped to Entrez: %d\n", nrow(el_map)))

# ── GO term exclusions (same as GA-matched build) ────────────────────────────

exclude_go <- c(
  "GO:0006915", "GO:0012501", "GO:0097190", "GO:0043065", "GO:0043066", "GO:0070265",
  "GO:0006955", "GO:0006954", "GO:0045087", "GO:0002376", "GO:0006952",
  "GO:0034097", "GO:0071345", "GO:0019221", "GO:0002250", "GO:0050776", "GO:0002684",
  "GO:0001890", "GO:0060136", "GO:0007565", "GO:0001824", "GO:0001829",
  "GO:0090074", "GO:0007566", "GO:0001837", "GO:0016477", "GO:0030335", "GO:0042060",
  "GO:0009755", "GO:0032870", "GO:0071383", "GO:0030539", "GO:0030520",
  "GO:0043627", "GO:0032355", "GO:0048545",
  "GO:0006950", "GO:0001666", "GO:0071456", "GO:0034599", "GO:0006979",
  "GO:0009615", "GO:0006974",
  "GO:0001525", "GO:0001568", "GO:0048514"
)

all_exclude_go <- unique(exclude_go)
for (go_id in exclude_go) {
  offspring <- tryCatch(
    as.character(GOBPOFFSPRING[[go_id]]),
    error = function(e) character(0)
  )
  if (length(offspring) > 0)
    all_exclude_go <- unique(c(all_exclude_go, offspring[!is.na(offspring)]))
}
cat(sprintf("GO terms (with offspring): %d\n", length(all_exclude_go)))

go2entrez <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db, keys = all_exclude_go,
                        columns = "ENTREZID", keytype = "GOALL")
)
exclude_entrez <- unique(go2entrez$ENTREZID[!is.na(go2entrez$ENTREZID)])

unstable_symbols <- c("ACTB", "GAPDH", "RPLP0")
unstable_map <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db, keys = unstable_symbols,
                        columns = "ENTREZID", keytype = "SYMBOL")
)
unstable_entrez <- unstable_map$ENTREZID[!is.na(unstable_map$ENTREZID)]

all_exclude <- unique(c(exclude_entrez, unstable_entrez))
keep <- el_map[!el_map$ENTREZID %in% all_exclude, ]
cat(sprintf("After GO + unstable exclusion: %d control genes\n\n", nrow(keep)))

# ── Dataset gene lists ────────────────────────────────────────────────────────

genes_12767  <- as.character(rownames(read.delim(
  "data/mapped/cvs/GSE12767_entrez_protein_coding.tsv", row.names = 1, check.names = FALSE)))
genes_100051 <- as.character(rownames(read.delim(
  "data/mapped/GSE100051.tsv", row.names = 1, check.names = FALSE)))
genes_93520  <- as.character(rownames(read.delim(
  "data/mapped/GSE93520.tsv", row.names = 1, check.names = FALSE)))
genes_28551  <- as.character(rownames(read.delim(
  "data/mapped/GSE28551.tsv", row.names = 1, check.names = FALSE)))

cat(sprintf("Gene counts: GSE12767=%d, GSE100051=%d, GSE93520=%d, GSE28551=%d\n",
            length(genes_12767), length(genes_100051), length(genes_93520), length(genes_28551)))

# ── Build per-pair control gene lists ─────────────────────────────────────────

pairs <- list(
  list(name = "GSE100051", genes = genes_100051),
  list(name = "GSE93520",  genes = genes_93520),
  list(name = "GSE28551",  genes = genes_28551)
)

for (p in pairs) {
  shared <- intersect(genes_12767, p$genes)
  ctl <- keep[keep$ENTREZID %in% shared, ]
  out_file <- file.path(out_dir, sprintf("control_genes_cvs_%s.txt", p$name))
  writeLines(ctl$ENTREZID, out_file)
  cat(sprintf("CVS + %s: %d shared genes, %d control genes -> %s\n",
              p$name, length(shared), nrow(ctl), out_file))
}
