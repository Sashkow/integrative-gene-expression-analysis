library(GEOquery)

raws_dir <- file.path(
  "data", "raws"
)

affy_datasets <- c(
  "GSE122214",  # GPL570, HG-U133 Plus 2.0, 4 samples
  "GSE22490",   # GPL570, HG-U133 Plus 2.0, 10 samples
  "GSE6573",    # GPL570, HG-U133 Plus 2.0, 6 samples
  "GSE37901",   # GPL570, HG-U133 Plus 2.0, 4 samples
  "GSE9984",    # GPL570, HG-U133 Plus 2.0, 12 samples
  "GSE73685",   # GPL6244, HuGene-1.0-st, 183 samples
  "GSE73374"    # GPL16686, HuGene-2.0-st, 36 samples
)

for (gse_id in affy_datasets) {
  dest_dir <- file.path(raws_dir, gse_id)
  cel_dir <- file.path(dest_dir, "cel")

  if (dir.exists(cel_dir) && length(list.files(cel_dir, pattern = "\\.[Cc][Ee][Ll]")) > 0) {
    cat(gse_id, ": CEL files already present, skipping\n")
    next
  }

  cat("\n=== Downloading", gse_id, "===\n")
  dir.create(cel_dir, recursive = TRUE, showWarnings = FALSE)

  tar_file <- file.path(dest_dir, paste0(gse_id, "_RAW.tar"))

  if (!file.exists(tar_file)) {
    tryCatch({
      getGEOSuppFiles(gse_id, baseDir = raws_dir, makeDirectory = TRUE)
    }, error = function(e) {
      cat("  ERROR downloading", gse_id, ":", conditionMessage(e), "\n")
      return(NULL)
    })
  }

  tar_file <- list.files(dest_dir, pattern = "_RAW\\.tar$", full.names = TRUE)
  if (length(tar_file) == 0) {
    cat("  No RAW.tar found for", gse_id, "\n")
    next
  }
  tar_file <- tar_file[1]

  cat("  Extracting", basename(tar_file), "to", cel_dir, "\n")
  untar(tar_file, exdir = cel_dir)

  gz_files <- list.files(cel_dir, pattern = "\\.gz$", full.names = TRUE)
  if (length(gz_files) > 0) {
    cat("  Decompressing", length(gz_files), "files...\n")
    for (f in gz_files) {
      system2("gunzip", args = c("-f", f))
    }
  }

  cel_count <- length(list.files(cel_dir, pattern = "\\.[Cc][Ee][Ll]$"))
  cat("  Done:", cel_count, "CEL files extracted\n")
}

cat("\n=== Summary ===\n")
for (gse_id in affy_datasets) {
  cel_dir <- file.path(raws_dir, gse_id, "cel")
  n <- length(list.files(cel_dir, pattern = "\\.[Cc][Ee][Ll]$"))
  cat(gse_id, ":", n, "CEL files\n")
}
