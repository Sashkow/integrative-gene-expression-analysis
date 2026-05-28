#!/usr/bin/env Rscript
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript build_qc_images_agilent.R <GSE_ID> [raw_dir] [output_dir]")

gse_id <- args[1]
raw_dir <- if (length(args) >= 2) args[2] else file.path("data", "raws", gse_id, "raw")
output_dir <- if (length(args) >= 3) args[3] else file.path("output", "qc", gse_id)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

txt_files <- list.files(raw_dir, pattern = "\\.txt(\\.gz)?$", full.names = TRUE)
cat("Found", length(txt_files), "Agilent FE files in", raw_dir, "\n")

read_agilent_fe <- function(filepath) {
  lines <- readLines(filepath, n = 500)
  features_line <- grep("^FEATURES", lines)
  if (length(features_line) == 0) stop("No FEATURES header found in ", filepath)
  header <- strsplit(lines[features_line], "\t")[[1]]
  data_start <- features_line + 1
  df <- read.delim(filepath, skip = data_start - 1, header = FALSE,
                   stringsAsFactors = FALSE, check.names = FALSE)
  df <- df[df[, 1] == "DATA", , drop = FALSE]
  colnames(df) <- header
  df$Row <- as.integer(df$Row)
  df$Col <- as.integer(df$Col)
  df$gMedianSignal <- as.numeric(df$gMedianSignal)
  df
}

for (filepath in txt_files) {
  fname <- basename(filepath)
  gsm_id <- sub("_.*", "", fname)
  cat("Processing", gsm_id, "...\n")

  fe <- tryCatch(read_agilent_fe(filepath), error = function(e) {
    cat("  ERROR:", e$message, "\n")
    NULL
  })
  if (is.null(fe)) next

  max_row <- max(fe$Row)
  max_col <- max(fe$Col)

  mat <- matrix(NA, nrow = max_row, ncol = max_col)
  for (i in seq_len(nrow(fe))) {
    mat[fe$Row[i], fe$Col[i]] <- fe$gMedianSignal[i]
  }
  mat_log <- log2(pmax(mat, 1))

  outfile <- file.path(output_dir, paste0("chip_image_", gsm_id, ".png"))
  png(outfile, width = max_col * 6, height = max_row * 2)
  par(mar = c(0, 0, 2, 0))
  image(t(mat_log[nrow(mat_log):1, ]), col = colorRampPalette(c("blue", "yellow", "red"))(256),
        axes = FALSE, main = paste0(gse_id, " — ", gsm_id))
  dev.off()
  cat("  Saved", basename(outfile), "\n")
}

n <- length(txt_files)
ncol_grid <- min(n, 4)
nrow_grid <- ceiling(n / ncol_grid)
grid_file <- file.path(output_dir, paste0("chip_images_all_", gse_id, ".png"))

png(grid_file, width = 500 * ncol_grid, height = 250 * nrow_grid)
par(mfrow = c(nrow_grid, ncol_grid), mar = c(0, 0, 2, 0))
for (filepath in txt_files) {
  fname <- basename(filepath)
  gsm_id <- sub("_.*", "", fname)
  fe <- tryCatch(read_agilent_fe(filepath), error = function(e) NULL)
  if (is.null(fe)) {
    plot.new()
    next
  }
  max_row <- max(fe$Row)
  max_col <- max(fe$Col)
  mat <- matrix(NA, nrow = max_row, ncol = max_col)
  for (i in seq_len(nrow(fe))) {
    mat[fe$Row[i], fe$Col[i]] <- fe$gMedianSignal[i]
  }
  mat_log <- log2(pmax(mat, 1))
  image(t(mat_log[nrow(mat_log):1, ]), col = colorRampPalette(c("blue", "yellow", "red"))(256),
        axes = FALSE, main = gsm_id)
}
dev.off()
cat("Saved combined grid:", basename(grid_file), "\n")
cat("Done.\n")
