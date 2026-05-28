suppressPackageStartupMessages({
  library(affy)
  library(makecdfenv)
  library(Biobase)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")
options(timeout = 1800)

dataset_id   <- "GSE12767"
raw_dir      <- file.path("data/raws", dataset_id)
cel_dir      <- file.path(raw_dir, "cel")
tar_path     <- file.path(raw_dir, paste0(dataset_id, "_RAW.tar"))
tar_url      <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE12nnn/GSE12767/suppl/GSE12767_RAW.tar"
cdf_pkg_name <- "hgu133plus2hsentrezgcdf"
cdf_internal <- "hgu133plus2hsentrezg"
cdf_zip      <- "/home/shivers/y/microarray-data-pipeline/common/HGU133Plus2_Hs_ENTREZG_25.0.0.zip"
cdf_cache    <- "data/raws/brainarray_cdf_cache"

dir.create(raw_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(cel_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(cdf_cache, showWarnings = FALSE, recursive = TRUE)

# 1. Download RAW.tar
if (!file.exists(tar_path)) {
  cat("Downloading", tar_url, "\n")
  download.file(tar_url, tar_path, mode = "wb")
}
cat("RAW tar:", tar_path, "(", file.info(tar_path)$size, "bytes)\n")

# 2. Extract CEL files
existing_cels <- list.files(cel_dir, pattern = "GSM.*\\.[cC][eE][lL]$", full.names = TRUE)
if (length(existing_cels) == 0) {
  cat("Extracting tar to", cel_dir, "\n")
  utils::untar(tar_path, exdir = cel_dir)
  gz_files <- list.files(cel_dir, pattern = "\\.gz$", full.names = TRUE)
  for (gz in gz_files) {
    out <- sub("\\.gz$", "", gz)
    R.utils::gunzip(gz, destname = out, overwrite = TRUE, remove = TRUE)
  }
}
cel_files <- sort(list.files(cel_dir, pattern = "GSM.*\\.[cC][eE][lL]$", full.names = TRUE))
cat("CEL files:", length(cel_files), "\n")

# 3. Ensure Brainarray CDF is installed
if (!requireNamespace(cdf_pkg_name, quietly = TRUE)) {
  if (!file.exists(cdf_zip)) stop("Brainarray ZIP not found: ", cdf_zip)
  extract_dir <- file.path(cdf_cache, "unzipped")
  pkg_dir     <- file.path(cdf_cache, cdf_pkg_name)
  dir.create(extract_dir, showWarnings = FALSE, recursive = TRUE)
  utils::unzip(cdf_zip, exdir = extract_dir, overwrite = TRUE)
  cdf_file <- list.files(extract_dir, pattern = "\\.cdf$", recursive = TRUE, full.names = TRUE)[1]
  if (is.na(cdf_file)) stop("No .cdf inside ", cdf_zip)
  if (!dir.exists(pkg_dir)) {
    makecdfenv::make.cdf.package(
      filename    = basename(cdf_file),
      packagename = cdf_pkg_name,
      cdf.path    = dirname(cdf_file),
      package.path= cdf_cache,
      species     = "Homo_sapiens",
      unlink      = TRUE
    )
  }
  utils::install.packages(pkg_dir, repos = NULL, type = "source")
}
stopifnot(requireNamespace(cdf_pkg_name, quietly = TRUE))
cat("Brainarray CDF ready\n")

# 4. Read CELs with Brainarray CDF and run RMA
cat("Reading CELs with cdfname =", cdf_internal, "\n")
ab <- affy::ReadAffy(filenames = cel_files, cdfname = cdf_internal)
cat("AffyBatch:", nrow(exprs(ab)), "probes x", ncol(exprs(ab)), "samples\n")

cat("Running RMA...\n")
eset <- affy::rma(ab)
expr <- exprs(eset)
cat("After RMA:", nrow(expr), "probesets x", ncol(expr), "samples\n")

# 5. Set sample names to GSM IDs with .CEL suffix (matching samples.csv)
gsm_ids <- sub("^(GSM[0-9]+).*$", "\\1.CEL", basename(cel_files))
colnames(expr) <- gsm_ids

# 6. Map Brainarray probesets to Entrez IDs
probesets <- rownames(expr)
entrez    <- sub("_at$", "", probesets)
keep      <- grepl("^[0-9]+$", entrez)
cat("Numeric ENTREZID probesets:", sum(keep), "/", length(keep), "\n")
expr      <- expr[keep, , drop = FALSE]
rownames(expr) <- entrez[keep]
expr      <- expr[!duplicated(rownames(expr)), , drop = FALSE]

# 7. Filter to protein-coding genes
if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  library(AnnotationDbi)
  gene_types <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
    keys = rownames(expr), columns = "GENETYPE", keytype = "ENTREZID")
  pc_genes <- gene_types$ENTREZID[gene_types$GENETYPE == "protein-coding" & !is.na(gene_types$GENETYPE)]
  expr_pc <- expr[rownames(expr) %in% pc_genes, , drop = FALSE]
  cat("Protein-coding genes:", nrow(expr_pc), "/", nrow(expr), "\n")
} else {
  expr_pc <- expr
  cat("org.Hs.eg.db not available, keeping all genes\n")
}

# 8. Save
out_all <- "data/mapped/cvs/GSE12767_entrez.tsv"
out_pc  <- "data/mapped/cvs/GSE12767_entrez_protein_coding.tsv"
write.table(expr, out_all, sep = "\t", quote = FALSE, col.names = NA)
write.table(expr_pc, out_pc, sep = "\t", quote = FALSE, col.names = NA)
cat("Wrote:", out_all, "(", nrow(expr), "genes x", ncol(expr), "samples)\n")
cat("Wrote:", out_pc, "(", nrow(expr_pc), "genes x", ncol(expr_pc), "samples)\n")
