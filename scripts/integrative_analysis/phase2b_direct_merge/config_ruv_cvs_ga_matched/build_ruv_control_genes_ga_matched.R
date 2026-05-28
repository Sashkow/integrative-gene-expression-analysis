#!/usr/bin/env Rscript
#' Build RUV negative control gene list for GA-matched CVS-vs-abortion
#' comparison using 3 datasets: GSE12767, GSE93520, GSE100051.
#'
#' Same logic as config_ruv_cvs/build_ruv_control_genes_cvs.R but
#' filters to genes present in ALL THREE datasets.

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

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

cat(sprintf("\nExcluding genes in %d GO BP terms\n", length(exclude_go)))

all_exclude_go <- unique(exclude_go)
for (go_id in exclude_go) {
  offspring <- tryCatch(
    as.character(GOBPOFFSPRING[[go_id]]),
    error = function(e) character(0)
  )
  if (length(offspring) > 0) {
    all_exclude_go <- unique(c(all_exclude_go, offspring[!is.na(offspring)]))
  }
}
cat(sprintf("Including offspring: %d total GO terms\n", length(all_exclude_go)))

go2entrez <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db,
    keys = all_exclude_go,
    columns = "ENTREZID",
    keytype = "GOALL"
  )
)
exclude_entrez <- unique(go2entrez$ENTREZID[!is.na(go2entrez$ENTREZID)])
cat(sprintf("Genes in excluded GO terms: %d unique Entrez IDs\n", length(exclude_entrez)))

unstable_symbols <- c("ACTB", "GAPDH", "RPLP0")
unstable_map <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db, keys = unstable_symbols,
                        columns = "ENTREZID", keytype = "SYMBOL")
)
unstable_entrez <- unstable_map$ENTREZID[!is.na(unstable_map$ENTREZID)]

all_exclude <- unique(c(exclude_entrez, unstable_entrez))
keep <- el_map[!el_map$ENTREZID %in% all_exclude, ]

cat(sprintf("\nEisenberg-Levanon: %d\n", nrow(el_map)))
cat(sprintf("Excluded by GO terms: %d\n", sum(el_map$ENTREZID %in% exclude_entrez)))
cat(sprintf("Excluded by literature (unstable): %d\n", sum(el_map$ENTREZID %in% unstable_entrez)))
cat(sprintf("Final control gene set: %d\n", nrow(keep)))

# Filter to genes present in ALL THREE datasets
genes_12767  <- rownames(read.delim("data/mapped/cvs/GSE12767_entrez_protein_coding.tsv",
                                     row.names = 1, check.names = FALSE))
genes_93520  <- rownames(read.delim("data/mapped/GSE93520.tsv",
                                     row.names = 1, check.names = FALSE))
genes_100051 <- rownames(read.delim("data/mapped/GSE100051.tsv",
                                     row.names = 1, check.names = FALSE))

shared_genes <- Reduce(intersect, list(
  as.character(genes_12767),
  as.character(genes_93520),
  as.character(genes_100051)
))
cat(sprintf("\nShared genes across 3 datasets: %d\n", length(shared_genes)))

keep_in_data <- keep[keep$ENTREZID %in% shared_genes, ]
cat(sprintf("Control genes in all 3 datasets: %d / %d\n", nrow(keep_in_data), nrow(keep)))

out_path <- "scripts/integrative_analysis/phase2b_direct_merge/config_ruv_cvs_ga_matched/ruv_control_genes_ga_matched.txt"
writeLines(keep_in_data$ENTREZID, out_path)
cat(sprintf("\nWrote: %s (%d genes)\n", out_path, nrow(keep_in_data)))
