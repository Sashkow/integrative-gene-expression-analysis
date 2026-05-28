suppressPackageStartupMessages({
  library(limma)
  library(illuminaHumanv4.db)
  library(AnnotationDbi)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")
options(timeout = 1800)

ds        <- "GSE55439"
raw_dir   <- file.path("data/raw_geo", ds)
nn_path   <- file.path(raw_dir, paste0(ds, "_non_normalized.txt.gz"))
sm_path   <- file.path(raw_dir, paste0(ds, "_series_matrix.txt.gz"))
out_path  <- "data/mapped/GSE55439_remapped.tsv"
cmp_path  <- "data/mapped/GSE55439.tsv"

stopifnot(file.exists(nn_path), file.exists(sm_path), file.exists(cmp_path))

# ---------- 1. SAMPLE label -> GSM accession ----------
sm_lines <- readLines(gzfile(sm_path), n = 200)
parse_sm <- function(tag) {
  line <- grep(paste0("^!Sample_", tag, "\\b"), sm_lines, value = TRUE)[1]
  toks <- strsplit(line, "\t")[[1]][-1]
  gsub("^\"|\"$", "", toks)
}
sample_desc <- parse_sm("description")
sample_gsm  <- parse_sm("geo_accession")
stopifnot(length(sample_desc) == length(sample_gsm))
label_to_gsm <- setNames(sample_gsm, sample_desc)
cat("Series-matrix samples:", length(sample_gsm), " first GSM:", sample_gsm[1], "\n")

# ---------- 2. Read non-normalized probe-level intensities ----------
cat("Reading", nn_path, "with limma::read.ilmn()\n")
raw <- read.ilmn(
  files     = nn_path,
  probeid   = "ID_REF",
  expr      = "SAMPLE",
  other.columns = "Detection"
)
cat("Probes:", nrow(raw), " samples:", ncol(raw), "\n")

# read.ilmn drops the "SAMPLE " prefix and keeps just the suffix
col_labels <- paste0("SAMPLE", colnames(raw))  # colnames already start with " 1", " 2", …
gsm_cols   <- label_to_gsm[col_labels]
stopifnot(!any(is.na(gsm_cols)))
colnames(raw) <- gsm_cols
cat("Mapped columns to GSM; first 4:", paste(head(colnames(raw), 4), collapse=", "), "\n")

# ---------- 3. neqc (normexp BG correct + quantile norm + log2) ----------
cat("Running limma::neqc() on all", ncol(raw), "samples (full dataset for best normalization)\n")
eset <- neqc(raw)
expr <- as.matrix(eset$E)
cat("After neqc:", nrow(expr), "probes x", ncol(expr), "samples  range:", round(range(expr), 2), "\n")

# ---------- 4. Filter by detection p-value (keep probes detected in >=1 sample) ----------
if (!is.null(eset$other$Detection)) {
  det <- as.matrix(eset$other$Detection)
  detected <- rowSums(det < 0.05) >= 1
  cat("Probes detected (p<0.05) in >=1 sample:", sum(detected), "/", nrow(expr), "\n")
  expr <- expr[detected, , drop = FALSE]
}

# ---------- 5. Map ILMN probes -> ENTREZID via illuminaHumanv4.db ----------
probes <- rownames(expr)
cat("Mapping", length(probes), "probes to ENTREZID via illuminaHumanv4.db\n")
map <- AnnotationDbi::select(
  illuminaHumanv4.db,
  keys = probes,
  columns = "ENTREZID",
  keytype = "PROBEID"
)
map <- map[!is.na(map$ENTREZID) & nzchar(map$ENTREZID), ]
map <- map[!duplicated(map$PROBEID), ]
cat("Probes with ENTREZID:", nrow(map), "\n")

expr <- expr[map$PROBEID, , drop = FALSE]
map$mean_expr <- rowMeans(expr, na.rm = TRUE)

# ---------- 6. Collapse probes->gene: pick highest-mean probe per gene ----------
ord <- order(map$ENTREZID, -map$mean_expr)
map <- map[ord, ]
map <- map[!duplicated(map$ENTREZID), ]
expr <- expr[map$PROBEID, , drop = FALSE]
rownames(expr) <- map$ENTREZID
cat("Unique ENTREZIDs after collapse:", nrow(expr), "\n")

# ---------- 7. Write remapped TSV (NEW path — does not overwrite) ----------
write.table(expr, out_path, sep = "\t", quote = FALSE, col.names = NA)
cat("Wrote:", out_path, "\n")
cat("Dims:", nrow(expr), "genes x", ncol(expr), "samples\n")

# ---------- 8. Compare against existing mapped file ----------
cat("\n=== Comparison vs", cmp_path, "===\n")
old <- read.table(cmp_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
cat("Old: ", nrow(old), "genes x", ncol(old), "samples\n")
cat("New: ", nrow(expr), "genes x", ncol(expr), "samples\n")
common_genes   <- intersect(rownames(old), rownames(expr))
only_in_old    <- setdiff(rownames(old), rownames(expr))
only_in_new    <- setdiff(rownames(expr), rownames(old))
cat("Genes in both:   ", length(common_genes), "\n")
cat("Only in old:     ", length(only_in_old), "\n")
cat("Only in new:     ", length(only_in_new), "\n")

# Correlation on shared 24 samples + shared genes
shared_samples <- intersect(colnames(old), colnames(expr))
cat("Shared samples:  ", length(shared_samples), "\n")
if (length(shared_samples) && length(common_genes)) {
  o <- as.matrix(old[common_genes, shared_samples])
  n <- expr[common_genes, shared_samples]
  cors <- sapply(seq_len(ncol(o)), function(i) cor(o[, i], n[, i]))
  cat("Per-sample Pearson correlation (old vs new) on shared genes:\n")
  cat("  mean:", round(mean(cors), 4), " min:", round(min(cors), 4), " max:", round(max(cors), 4), "\n")
}
