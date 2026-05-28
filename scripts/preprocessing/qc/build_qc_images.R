library(affy)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript build_qc_images.R <GSE_ID> [cel_dir] [output_dir]")

gse_id <- args[1]
cel_dir <- if (length(args) >= 2) args[2] else file.path("data", "raws", gse_id, "cel")
output_dir <- if (length(args) >= 3) args[3] else file.path("output", "qc", gse_id)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading CEL files from", cel_dir, "\n")
cel_files <- list.files(cel_dir, pattern = "\\.[Cc][Ee][Ll]$", full.names = TRUE)
cat("Found", length(cel_files), "CEL files\n")

abatch <- ReadAffy(filenames = cel_files)
sample_names <- sampleNames(abatch)

nrow_array <- nrow(abatch)
ncol_array <- ncol(abatch)
cat("Array dimensions:", nrow_array, "x", ncol_array, "\n")

# --- Pseudo-images (chip scan reproductions) ---
cat("\nGenerating chip pseudo-images...\n")
for (i in seq_along(sample_names)) {
  sname <- sample_names[i]
  clean_name <- gsub("\\.[Cc][Ee][Ll]$", "", sname)
  outfile <- file.path(output_dir, paste0("chip_image_", clean_name, ".png"))

  png(outfile, width = ncol_array, height = nrow_array)
  par(mar = c(0, 0, 2, 0))
  image(abatch[, i], main = paste0(gse_id, " — ", clean_name))
  dev.off()
  cat("  Saved", basename(outfile), "\n")
}

# --- Combined pseudo-image grid ---
n <- length(sample_names)
ncol_grid <- min(n, 3)
nrow_grid <- ceiling(n / ncol_grid)
grid_file <- file.path(output_dir, paste0("chip_images_all_", gse_id, ".png"))

panel_w <- ncol_array
panel_h <- nrow_array
png(grid_file, width = panel_w * ncol_grid, height = panel_h * nrow_grid)
par(mfrow = c(nrow_grid, ncol_grid), mar = c(0, 0, 2, 0))
for (i in seq_along(sample_names)) {
  clean_name <- gsub("\\.[Cc][Ee][Ll]$", "", sample_names[i])
  image(abatch[, i], main = clean_name)
}
dev.off()
cat("Saved combined grid:", basename(grid_file), "\n")

# --- RNA degradation plot ---
cat("\nComputing RNA degradation...\n")
rna_deg <- AffyRNAdeg(abatch)

deg_file <- file.path(output_dir, paste0("rna_degradation_", gse_id, ".png"))
png(deg_file, width = 800, height = 600)
plotAffyRNAdeg(rna_deg, cols = rainbow(n))
title(main = paste0("RNA Degradation Plot — ", gse_id))
legend("topleft",
  legend = gsub("\\.[Cc][Ee][Ll]$", "", sample_names),
  col = rainbow(n), lty = 1, cex = 0.7
)
dev.off()
cat("Saved", basename(deg_file), "\n")

# --- RNA degradation plot with ggplot ---
deg_data <- data.frame(
  probe_position = rep(seq_len(ncol(rna_deg$means.by.number)), each = n),
  mean_intensity = as.vector(t(rna_deg$means.by.number)),
  sample = rep(gsub("\\.[Cc][Ee][Ll]$", "", sample_names), ncol(rna_deg$means.by.number))
)

deg_gg_file <- file.path(output_dir, paste0("rna_degradation_ggplot_", gse_id, ".png"))
p <- ggplot(deg_data, aes(x = probe_position, y = mean_intensity, colour = sample)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  labs(
    title = paste0("RNA Degradation — ", gse_id),
    subtitle = "Mean intensity by probe position (5' to 3')",
    x = "Probe Position (5' → 3')",
    y = "Mean Intensity",
    colour = "Sample"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right")
png(deg_gg_file, width = 1500, height = 900, res = 150)
print(p)
dev.off()
cat("Saved", basename(deg_gg_file), "\n")

# --- Degradation slope summary ---
cat("\nRNA Degradation Slopes (3'/5' ratio proxy):\n")
slopes <- rna_deg$slope
names(slopes) <- gsub("\\.[Cc][Ee][Ll]$", "", sample_names)
for (i in seq_along(slopes)) {
  cat(sprintf("  %s: %.4f\n", names(slopes)[i], slopes[i]))
}

cat("\nDone. All QC images saved to", output_dir, "\n")
