#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript build_qc_images_nimblegen.R <GSE_ID> [raw_dir] [output_dir]")

gse_id <- args[1]
raw_dir <- if (length(args) >= 2) args[2] else file.path("data", "raws", gse_id, "raw")
output_dir <- if (length(args) >= 3) args[3] else file.path("output", "qc", gse_id)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pair_files <- list.files(raw_dir, pattern = "\\.pair(\\.gz)?$", full.names = TRUE)
cat("Found", length(pair_files), "NimbleGen PAIR files in", raw_dir, "\n")

read_pair <- function(filepath) {
  df <- read.delim(filepath, comment.char = "#", header = TRUE,
                   stringsAsFactors = FALSE, check.names = FALSE)
  df$X <- as.integer(df$X)
  df$Y <- as.integer(df$Y)
  df$PM <- as.numeric(df$PM)
  df
}

for (filepath in pair_files) {
  fname <- basename(filepath)
  gsm_id <- sub("_.*", "", fname)
  cat("Processing", gsm_id, "...\n")

  pair <- tryCatch(read_pair(filepath), error = function(e) {
    cat("  ERROR:", e$message, "\n")
    NULL
  })
  if (is.null(pair)) next

  max_x <- max(pair$X, na.rm = TRUE)
  max_y <- max(pair$Y, na.rm = TRUE)

  mat <- matrix(NA, nrow = max_y + 1, ncol = max_x + 1)
  for (i in seq_len(nrow(pair))) {
    mat[pair$Y[i] + 1, pair$X[i] + 1] <- pair$PM[i]
  }
  mat_log <- log2(pmax(mat, 1))

  outfile <- file.path(output_dir, paste0("chip_image_", gsm_id, ".png"))
  png(outfile, width = (max_x + 1) * 2, height = (max_y + 1) * 2)
  par(mar = c(0, 0, 2, 0))
  image(t(mat_log[nrow(mat_log):1, ]), col = colorRampPalette(c("blue", "yellow", "red"))(256),
        axes = FALSE, main = paste0(gse_id, " - ", gsm_id))
  dev.off()
  cat("  Saved", basename(outfile), "\n")
}

n <- length(pair_files)
ncol_grid <- min(n, 5)
nrow_grid <- ceiling(n / ncol_grid)
grid_file <- file.path(output_dir, paste0("chip_images_all_", gse_id, ".png"))

png(grid_file, width = 400 * ncol_grid, height = 550 * nrow_grid)
par(mfrow = c(nrow_grid, ncol_grid), mar = c(0, 0, 2, 0))
for (filepath in pair_files) {
  fname <- basename(filepath)
  gsm_id <- sub("_.*", "", fname)
  pair <- tryCatch(read_pair(filepath), error = function(e) NULL)
  if (is.null(pair)) {
    plot.new()
    next
  }
  max_x <- max(pair$X, na.rm = TRUE)
  max_y <- max(pair$Y, na.rm = TRUE)
  mat <- matrix(NA, nrow = max_y + 1, ncol = max_x + 1)
  for (i in seq_len(nrow(pair))) {
    mat[pair$Y[i] + 1, pair$X[i] + 1] <- pair$PM[i]
  }
  mat_log <- log2(pmax(mat, 1))
  image(t(mat_log[nrow(mat_log):1, ]), col = colorRampPalette(c("blue", "yellow", "red"))(256),
        axes = FALSE, main = gsm_id)
}
dev.off()
cat("Saved combined grid:", basename(grid_file), "\n")
cat("Done.\n")
