suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")
options(timeout = 600)

dest <- "data/raw_geo"
dir.create(dest, showWarnings = FALSE, recursive = TRUE)

sm_file <- file.path(dest, "GSE70102_series_matrix.txt.gz")
if (!file.exists(sm_file)) {
  url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70102/matrix/GSE70102_series_matrix.txt.gz"
  download.file(url, sm_file, mode = "wb")
}
cat("Loading series matrix from", sm_file, "\n")
gset <- getGEO(filename = sm_file, getGPL = FALSE)

annot_file <- file.path(dest, "GPL570.annot.gz")
if (!file.exists(annot_file)) {
  download.file(
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
    annot_file, mode = "wb"
  )
}
cat("Parsing", annot_file, "\n")
raw <- readLines(gzfile(annot_file))
start <- which(raw == "!platform_table_begin") + 1L
end_idx <- which(raw == "!platform_table_end") - 1L
gpl_tbl <- read.table(
  text = raw[start:end_idx], header = TRUE, sep = "\t",
  quote = "", comment.char = "", stringsAsFactors = FALSE,
  check.names = FALSE, fill = TRUE
)
cat("GPL570 annot rows:", nrow(gpl_tbl), " cols:", paste(colnames(gpl_tbl)[1:5], collapse = ", "), "...\n")
fData(gset) <- gpl_tbl[match(rownames(exprs(gset)), gpl_tbl$ID), , drop = FALSE]

exprs_data <- exprs(gset)
cat("Probes:", nrow(exprs_data), "  Samples:", ncol(exprs_data), "\n")
cat("Sample columns:", paste(head(colnames(exprs_data), 4), collapse = ", "), "...\n")

qx <- as.numeric(quantile(exprs_data, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
needs_log <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (needs_log) {
  cat("Applying log2(x+1) transform — values look raw.\n")
  exprs_data[exprs_data < 0] <- 0
  exprs_data <- log2(exprs_data + 1)
} else {
  cat("Data already on log scale.\n")
}

fd <- fData(gset)
cat("fData columns:", paste(colnames(fd), collapse = " | "), "\n")

entrez_col <- intersect(c("ENTREZ_GENE_ID", "Gene ID", "GENE", "Entrez_Gene_ID"), colnames(fd))[1]
if (is.na(entrez_col)) stop("No Entrez ID column found in platform annotation")
cat("Using fData column:", entrez_col, "\n")

probeid_col <- intersect(c("ID", "PROBEID", "Probe Set ID"), colnames(fd))[1]
if (is.na(probeid_col)) stop("No probe ID column found in platform annotation")
ann <- data.frame(
  PROBEID  = as.character(fd[[probeid_col]]),
  ENTREZID = as.character(fd[[entrez_col]]),
  stringsAsFactors = FALSE
)
ann <- ann[!is.na(ann$ENTREZID) & ann$ENTREZID != "" & ann$ENTREZID != "null", , drop = FALSE]
ann <- ann[!grepl("///", ann$ENTREZID, fixed = TRUE), , drop = FALSE]
ann <- ann[ann$PROBEID %in% rownames(exprs_data), , drop = FALSE]
cat("Probes with single ENTREZID:", nrow(ann), "\n")

exprs_mapped <- exprs_data[ann$PROBEID, , drop = FALSE]
ann$MEAN_EXPR <- rowMeans(exprs_mapped, na.rm = TRUE)
ann <- ann[order(ann$ENTREZID, -ann$MEAN_EXPR), ]
ann <- ann[!duplicated(ann$ENTREZID), ]
cat("Unique genes after collapse:", nrow(ann), "\n")

exprs_final <- exprs_mapped[ann$PROBEID, , drop = FALSE]
rownames(exprs_final) <- ann$ENTREZID

out <- "data/mapped/GSE70102.tsv"
write.table(exprs_final, out, sep = "\t", quote = FALSE, col.names = NA)
cat("Wrote:", out, "\n")
cat("Dims:", nrow(exprs_final), "genes x", ncol(exprs_final), "samples\n")
