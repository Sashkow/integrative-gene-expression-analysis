suppressPackageStartupMessages({
  library(affy)
  library(makecdfenv)
  library(Biobase)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")
options(timeout = 1800)

dataset_id   <- "GSE70102"
raw_dir      <- file.path("data/raw_geo", dataset_id)
cel_dir      <- file.path(raw_dir, "cel")
tar_path     <- file.path(raw_dir, paste0(dataset_id, "_RAW.tar"))
tar_url      <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70102/suppl/GSE70102_RAW.tar"
cdf_zip      <- "/home/shivers/y/microarray-data-pipeline/common/HGU133Plus2_Hs_ENTREZG_25.0.0.zip"
cdf_pkg_name <- "hgu133plus2hsentrezgcdf"
cdf_internal <- "hgu133plus2hsentrezg"
cdf_cache    <- "data/raw_geo/brainarray_cdf_cache"

dir.create(raw_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(cel_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(cdf_cache, showWarnings = FALSE, recursive = TRUE)

# 1. Download GSE70102_RAW.tar if missing
if (!file.exists(tar_path)) {
  cat("Downloading", tar_url, "\n")
  download.file(tar_url, tar_path, mode = "wb")
}
cat("RAW tar:", tar_path, "(", file.info(tar_path)$size, "bytes )\n")

# 2. Extract CEL.gz files and decompress
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
print(basename(cel_files))

# 3. Ensure Brainarray ENTREZG CDF package is installed
if (!requireNamespace(cdf_pkg_name, quietly = TRUE)) {
  if (!file.exists(cdf_zip)) stop("Brainarray ZIP not found: ", cdf_zip)
  extract_dir <- file.path(cdf_cache, "unzipped")
  pkg_dir     <- file.path(cdf_cache, cdf_pkg_name)
  dir.create(extract_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Unzipping Brainarray archive:", cdf_zip, "\n")
  utils::unzip(cdf_zip, exdir = extract_dir, overwrite = TRUE)
  cdf_file <- list.files(extract_dir, pattern = "\\.cdf$", recursive = TRUE, full.names = TRUE)[1]
  if (is.na(cdf_file)) stop("No .cdf inside ", cdf_zip)
  cat("Building CDF package from", cdf_file, "\n")
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
  cat("Installing CDF package from source\n")
  utils::install.packages(pkg_dir, repos = NULL, type = "source")
}
stopifnot(requireNamespace(cdf_pkg_name, quietly = TRUE))
cat("Brainarray CDF ready:", cdf_pkg_name, "\n")

# 4. Read CELs with Brainarray CDF and run RMA
cat("Reading CELs with cdfname =", cdf_internal, "\n")
ab <- affy::ReadAffy(filenames = cel_files, cdfname = cdf_internal)
cat("AffyBatch:", nrow(exprs(ab)), "probes x", ncol(exprs(ab)), "samples\n")

cat("Running RMA (background correction + quantile normalization + log2 + median polish summarization)...\n")
eset <- affy::rma(ab)
expr <- exprs(eset)
cat("After RMA:", nrow(expr), "probesets x", ncol(expr), "samples\n")

# 5. Sample names: extract GSM accession from filename
gsm_ids <- sub("^(GSM[0-9]+).*$", "\\1", basename(cel_files))
colnames(expr) <- gsm_ids
cat("Sample columns:", paste(head(gsm_ids, 4), collapse=", "), "...\n")

# 6. Brainarray ENTREZG probesets are named like "780_at" → strip suffix to get ENTREZID
probesets <- rownames(expr)
entrez    <- sub("_at$", "", probesets)
keep      <- grepl("^[0-9]+$", entrez)
cat("Probesets keeping (numeric ENTREZID):", sum(keep), "/", length(keep), "\n")
expr      <- expr[keep, , drop = FALSE]
rownames(expr) <- entrez[keep]
expr      <- expr[!duplicated(rownames(expr)), , drop = FALSE]
cat("Final unique ENTREZIDs:", nrow(expr), "\n")

# 7. Save mapped TSV
out <- "data/mapped/GSE70102.tsv"
write.table(expr, out, sep = "\t", quote = FALSE, col.names = NA)
cat("Wrote:", out, "\n")
cat("Dims:", nrow(expr), "genes x", ncol(expr), "samples\n")
