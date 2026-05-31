#!/usr/bin/env Rscript
#' Build RUV negative control gene list for CVS-vs-abortion comparison.
#'
#' Starting set: Eisenberg-Levanon housekeeping genes (same as main pipeline).
#'
#' Exclusion criteria — remove genes in GO categories likely affected by
#' CVS (ongoing pregnancy) vs abortion (terminated pregnancy):
#'   - Apoptosis / programmed cell death
#'   - Immune / inflammatory response
#'   - Trophoblast invasion / placental development
#'   - Hormone signaling (hCG, progesterone, estrogen)
#'   - Stress response / hypoxia
#'
#' Also excludes ACTB, GAPDH, RPLP0 (known unstable in placenta) per
#' the same literature verdicts used in the main pipeline.
#'
#' Output: ruv_control_genes_cvs.txt (Entrez IDs, one per line)

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

# GO terms to exclude: biological processes affected by CVS vs abortion
exclude_go <- c(
  # Apoptosis / cell death
  "GO:0006915",   # apoptotic process
  "GO:0012501",   # programmed cell death
  "GO:0097190",   # apoptotic signaling pathway
  "GO:0043065",   # positive regulation of apoptotic process
  "GO:0043066",   # negative regulation of apoptotic process
  "GO:0070265",   # necrotic cell death

  # Immune / inflammatory response
  "GO:0006955",   # immune response
  "GO:0006954",   # inflammatory response
  "GO:0045087",   # innate immune response
  "GO:0002376",   # immune system process
  "GO:0006952",   # defense response
  "GO:0034097",   # response to cytokine
  "GO:0071345",   # cellular response to cytokine stimulus
  "GO:0019221",   # cytokine-mediated signaling pathway
  "GO:0002250",   # adaptive immune response
  "GO:0050776",   # regulation of immune response
  "GO:0002684",   # positive regulation of immune system process

  # Trophoblast / placental development
  "GO:0001890",   # placenta development
  "GO:0060136",   # embryonic process involved in female pregnancy
  "GO:0007565",   # female pregnancy
  "GO:0001824",   # blastocyst development
  "GO:0001829",   # trophectodermal cell differentiation
  "GO:0090074",   # negative regulation of protein homodimerization activity
  "GO:0007566",   # embryo implantation
  "GO:0001837",   # epithelial to mesenchymal transition
  "GO:0016477",   # cell migration
  "GO:0030335",   # positive regulation of cell migration
  "GO:0042060",   # wound healing

  # Hormone signaling
  "GO:0009755",   # hormone-mediated signaling pathway
  "GO:0032870",   # cellular response to hormone stimulus
  "GO:0071383",   # cellular response to steroid hormone stimulus
  "GO:0030539",   # male sex determination
  "GO:0030520",   # intracellular estrogen receptor signaling pathway
  "GO:0043627",   # response to estrogen
  "GO:0032355",   # response to estradiol
  "GO:0048545",   # response to steroid hormone

  # Stress / hypoxia
  "GO:0006950",   # response to stress
  "GO:0001666",   # response to hypoxia
  "GO:0071456",   # cellular response to hypoxia
  "GO:0034599",   # cellular response to oxidative stress
  "GO:0006979",   # response to oxidative stress
  "GO:0009615",   # response to virus
  "GO:0006974",   # cellular response to DNA damage stimulus

  # Angiogenesis (active in ongoing pregnancy)
  "GO:0001525",   # angiogenesis
  "GO:0001568",   # blood vessel development
  "GO:0048514"    # blood vessel morphogenesis
)

cat(sprintf("\nExcluding genes in %d GO BP terms\n", length(exclude_go)))

# Get all offspring terms (children) for each GO term
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

# Map GO terms to Entrez IDs
go2entrez <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db,
    keys = all_exclude_go,
    columns = "ENTREZID",
    keytype = "GOALL"
  )
)
exclude_entrez <- unique(go2entrez$ENTREZID[!is.na(go2entrez$ENTREZID)])
cat(sprintf("Genes in excluded GO terms: %d unique Entrez IDs\n", length(exclude_entrez)))

# Literature-unstable genes (same as main pipeline)
unstable_symbols <- c("ACTB", "GAPDH", "RPLP0")
unstable_map <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db, keys = unstable_symbols,
                        columns = "ENTREZID", keytype = "SYMBOL")
)
unstable_entrez <- unstable_map$ENTREZID[!is.na(unstable_map$ENTREZID)]

# Filter
all_exclude <- unique(c(exclude_entrez, unstable_entrez))
keep <- el_map[!el_map$ENTREZID %in% all_exclude, ]

cat(sprintf("\nEisenberg-Levanon: %d\n", nrow(el_map)))
cat(sprintf("Excluded by GO terms: %d\n", sum(el_map$ENTREZID %in% exclude_entrez)))
cat(sprintf("Excluded by literature (unstable): %d\n", sum(el_map$ENTREZID %in% unstable_entrez)))
cat(sprintf("Final control gene set: %d\n", nrow(keep)))

# Also filter to genes present in both datasets
cvs_genes <- rownames(read.delim("data/mapped/cvs/GSE12767_entrez_protein_coding.tsv",
                                  row.names = 1, nrows = 1, check.names = FALSE))
cvs_genes <- c(cvs_genes,
               rownames(read.delim("data/mapped/cvs/GSE12767_entrez_protein_coding.tsv",
                                    row.names = 1, check.names = FALSE)))
abort_genes <- rownames(read.delim("data/mapped/GSE93520.tsv",
                                    row.names = 1, check.names = FALSE))
shared_genes <- intersect(as.character(cvs_genes), as.character(abort_genes))

keep_in_data <- keep[keep$ENTREZID %in% shared_genes, ]
cat(sprintf("Present in both datasets: %d / %d\n", nrow(keep_in_data), nrow(keep)))

out_path <- "scripts/integrative_analysis/phase2b_direct_merge/config_ruv_cvs/ruv_control_genes_cvs.txt"
writeLines(keep_in_data$ENTREZID, out_path)
cat(sprintf("\nWrote: %s (%d genes)\n", out_path, nrow(keep_in_data)))

# Summary by excluded category
cat("\n=== Exclusion breakdown (HK genes only) ===\n")
for (go_id in exclude_go[1:20]) {
  go_name <- tryCatch(Term(GOTERM[[go_id]]), error = function(e) go_id)
  n_hk <- sum(el_map$ENTREZID %in% go2entrez$ENTREZID[go2entrez$GOALL == go_id])
  if (n_hk > 0) cat(sprintf("  %s (%s): %d HK genes\n", go_id, go_name, n_hk))
}
