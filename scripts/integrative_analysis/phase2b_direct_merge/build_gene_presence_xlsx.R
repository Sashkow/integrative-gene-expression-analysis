#!/usr/bin/env Rscript
#' Build per-gene presence/DEG table across all phase2b runs.
#'
#' Takes the union of protein-coding genes from the run with the largest
#' gene universe and produces an XLSX with columns:
#'   entrez_id, symbol, is_in_eisenberg_levanon,
#'   is_in_{run}, is_in_{run}_degs   (for each run)

suppressPackageStartupMessages({
  library(openxlsx)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

fdr_cutoff <- 0.05

runs <- list(
  list(id = "1_2_restoration_blockmask_imputed_0",
       dir = "output/phase2b_combat/phase2b_1_2_restoration_blockmask_imputed_0"),
  list(id = "2_3_restoration_blockmask_imputed_0",
       dir = "output/phase2b_combat/phase2b_2_3_restoration_blockmask_imputed_0"),
  list(id = "1_2_2nd_trim_only",
       dir = "output/phase2b_combat/phase2b_1_2_2nd_trim_only"),
  list(id = "2_3_2nd_trim_only",
       dir = "output/phase2b_combat/phase2b_2_3_2nd_trim_only"),
  list(id = "1_2_all_datasets",
       dir = "output/phase2b_combat/phase2b_1_2_all_datasets"),
  list(id = "2_3_all_datasets",
       dir = "output/phase2b_combat/phase2b_2_3_all_datasets"),
  list(id = "1_2_balanced",
       dir = "output/phase2b_combat/phase2b_1_2_balanced"),
  list(id = "2_3_balanced",
       dir = "output/phase2b_combat/phase2b_2_3_balanced"),
  list(id = "1_2_all_datasets_batch_in_limma",
       dir = "output/phase2b_batch_in_limma/phase2b_1_2_all_datasets_batch_in_limma"),
  list(id = "2_3_all_datasets_batch_in_limma",
       dir = "output/phase2b_batch_in_limma/phase2b_2_3_all_datasets_batch_in_limma"),
  list(id = "1_2_balanced_ruv",
       dir = "output/phase2b_ruv/phase2b_1_2_balanced_ruv"),
  list(id = "2_3_balanced_ruv",
       dir = "output/phase2b_ruv/phase2b_2_3_balanced_ruv"),
  list(id = "1_2_all_datasets_ruv",
       dir = "output/phase2b_ruv/phase2b_1_2_all_datasets_ruv"),
  list(id = "2_3_all_datasets_ruv",
       dir = "output/phase2b_ruv/phase2b_2_3_all_datasets_ruv")
)

runs <- runs[sapply(runs, function(r) dir.exists(r$dir))]
cat(sprintf("Found %d runs with output\n", length(runs)))

find_de_file <- function(run_dir) {
  candidates <- list.files(run_dir, pattern = "^difexp_(none|softimpute)_[^_].*\\.tsv$")
  candidates <- candidates[!grepl("^difexp_significant_", candidates)]
  if (length(candidates) == 0) return(NULL)
  # Prefer none (no imputation) over softimpute
  none_f <- grep("^difexp_none_", candidates, value = TRUE)
  if (length(none_f) > 0) return(none_f[1])
  candidates[1]
}

run_genes <- list()
run_degs  <- list()

for (r in runs) {
  exprs_file <- file.path(r$dir, "exprs_imputed_none.tsv")
  if (!file.exists(exprs_file)) next
  lines <- readLines(exprs_file)
  all_genes <- sub("\t.*", "", lines[-1])
  run_genes[[r$id]] <- unique(all_genes)

  de_file <- find_de_file(r$dir)
  if (!is.null(de_file)) {
    de <- read.table(file.path(r$dir, de_file), header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, check.names = FALSE, quote = "")
    sig <- de$gene[!is.na(de$adj.P.Val) & de$adj.P.Val < fdr_cutoff]
    run_degs[[r$id]] <- as.character(sig)
    cat(sprintf("  %s: %d genes, %d DEGs (FDR<%.2f) from %s\n",
                r$id, length(all_genes), length(sig), fdr_cutoff, de_file))
  } else {
    run_degs[[r$id]] <- character(0)
    cat(sprintf("  %s: %d genes, no DE file found\n", r$id, length(all_genes)))
  }
}

gene_counts <- sapply(run_genes, length)
largest_run <- names(which.max(gene_counts))
cat(sprintf("\nLargest gene universe: %s (%d genes)\n", largest_run, max(gene_counts)))

all_genes <- unique(unlist(run_genes))
cat(sprintf("Union of all genes across runs: %d\n", length(all_genes)))

# Map Entrez IDs to symbols
sym_map <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = all_genes,
                                  columns = "SYMBOL",
                                  keytype = "ENTREZID")
entrez_to_symbol <- setNames(sym_map$SYMBOL, sym_map$ENTREZID)

# Load Eisenberg-Levanon housekeeping gene list (gene symbols)
el_file <- "data/reference/integration_methods_references/ruv_housekeeping_genes_search/eisenberg_levanon_HK_genes.txt"
el_raw <- read.delim(el_file, header = FALSE, stringsAsFactors = FALSE,
                     strip.white = TRUE)
el_symbols <- trimws(el_raw$V1)
cat(sprintf("Eisenberg-Levanon list: %d symbols\n", length(el_symbols)))

# Map EL symbols to Entrez for matching
el_entrez_map <- AnnotationDbi::select(org.Hs.eg.db,
                                        keys = el_symbols,
                                        columns = "ENTREZID",
                                        keytype = "SYMBOL")
el_entrez <- unique(na.omit(el_entrez_map$ENTREZID))
cat(sprintf("Eisenberg-Levanon mapped to Entrez: %d IDs\n", length(el_entrez)))

# Build output data frame
df <- data.frame(
  entrez_id = all_genes,
  symbol = entrez_to_symbol[all_genes],
  is_in_eisenberg_levanon = as.integer(all_genes %in% el_entrez),
  stringsAsFactors = FALSE
)

for (r in runs) {
  col_present <- paste0("is_in_", r$id)
  col_deg     <- paste0("is_deg_", r$id)
  df[[col_present]] <- as.integer(all_genes %in% run_genes[[r$id]])
  df[[col_deg]]     <- as.integer(all_genes %in% run_degs[[r$id]])
}

df <- df[order(df$symbol), ]
rownames(df) <- NULL

cat(sprintf("\nFinal table: %d genes x %d columns\n", nrow(df), ncol(df)))
cat(sprintf("Eisenberg-Levanon overlap: %d / %d genes\n",
            sum(df$is_in_eisenberg_levanon), nrow(df)))

wb <- createWorkbook()
addWorksheet(wb, "gene_presence")
writeData(wb, "gene_presence", df, withFilter = TRUE)
freezePane(wb, "gene_presence", firstActiveRow = 2, firstActiveCol = 4)
setColWidths(wb, "gene_presence", cols = seq_len(ncol(df)), widths = "auto")

out_path <- "articles/imputation_article/phase2b_gene_presence.xlsx"
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
saveWorkbook(wb, out_path, overwrite = TRUE)
cat(sprintf("Wrote: %s\n", out_path))
