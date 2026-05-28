#!/usr/bin/env Rscript

#' PCA 2x2 Multiplot for 1st vs 2nd Trimester
#' 5 datasets only: GSE122214, GSE22490, GSE37901, GSE9984 (Affymetrix), GSE93520 (Agilent)
#'
#' Uses pipeline output but filters to only the 5 requested datasets

cat("\n=== PCA Multiplot: 5 Datasets, 1st vs 2nd Trimester ===\n\n")

library(ggplot2)
library(gridExtra)
library(grid)
library(sva)

# Source pipeline modules
source("scripts/integrative_analysis/integrative_analysis_disser_pipeline/config.R")
source("scripts/integrative_analysis/integrative_analysis_disser_pipeline/utils.R")
source("scripts/integrative_analysis/integrative_analysis_disser_pipeline/01_data_merging.R")
source("scripts/integrative_analysis/integrative_analysis_disser_pipeline/02_batch_correction.R")

# Config
config <- load_config("config/volodymyr/config_volodymyr_1_2_trim.yaml")
target_datasets <- c("GSE122214", "GSE22490", "GSE37901", "GSE9984", "GSE93520")
output_dir <- "output/volodymyr/volodymyr_1_2_trim/5_datasets_only"
dir.create(file.path(output_dir, "pcas"), showWarnings = FALSE, recursive = TRUE)

# Load and filter phenodata to only target datasets
cat("Loading phenodata...\n")
pdata_full <- read.csv(config$paths$phenodata, stringsAsFactors = FALSE)

# Apply standard filters
for (filter_col in names(config$filtering)) {
  allowed <- config$filtering[[filter_col]]
  pdata_full <- pdata_full[pdata_full[[filter_col]] %in% allowed, ]
}

# Restrict to target datasets only
pdata <- pdata_full[pdata_full$secondaryaccession %in% target_datasets, ]
cat("Samples after filtering:", nrow(pdata), "\n")
cat("Datasets:", paste(unique(pdata$secondaryaccession), collapse = ", "), "\n")
cat("Trimesters:", table(pdata$Gestational.Age.Category), "\n\n")

# Load and merge expression data for target datasets only
cat("Loading expression data...\n")
mapped_path <- config$paths$mapped_data
sample_col <- config$merging$sample_column

# Use pipeline's merge function for selected datasets only
cat("Merging expression data for target datasets...\n")
ds_files <- paste0(target_datasets, ".tsv")
first_file <- ds_files[1]
merged_exprs <- read.table(file.path(mapped_path, first_file),
                           header = TRUE, sep = "\t", row.names = 1)
cat("  Starting with", nrow(merged_exprs), "genes from", first_file, "\n")

for (f in ds_files[-1]) {
  current <- read.table(file.path(mapped_path, f),
                        header = TRUE, sep = "\t", row.names = 1)
  merged_exprs <- merge(merged_exprs, current, by = "row.names", all = FALSE)
  rownames(merged_exprs) <- merged_exprs$Row.names
  merged_exprs <- merged_exprs[, !(colnames(merged_exprs) == "Row.names")]
  cat("  Merged", f, "-", nrow(merged_exprs), "common genes\n")
}

# Align phenodata with expression data (using make.names like the pipeline does)
pdata <- pdata[make.names(pdata[[sample_col]]) %in% colnames(merged_exprs), ]
merged_exprs <- merged_exprs[, make.names(pdata[[sample_col]])]

cat("Final matrix:", nrow(merged_exprs), "genes x", ncol(merged_exprs), "samples\n\n")

# Create trim_term
pdata$trim_term <- pdata$Gestational.Age.Category

# ComBat batch correction
cat("Running ComBat...\n")
batch <- pdata$secondaryaccession
mod <- model.matrix(~trim_term + Combined.Fetus.Sex, data = pdata)

exprs_before <- merged_exprs
exprs_after <- ComBat(dat = as.matrix(merged_exprs), batch = batch, mod = mod,
                       par.prior = TRUE, prior.plots = FALSE)

# Save data
write.table(exprs_before, file.path(output_dir, "merged_exprs_before_combat.tsv"), sep = "\t", quote = FALSE)
write.table(exprs_after, file.path(output_dir, "merged_exprs_after_combat.tsv"), sep = "\t", quote = FALSE)
write.csv(pdata, file.path(output_dir, "phenodata.csv"), row.names = FALSE)

cat("\n=== Creating PCA Plots ===\n\n")

create_pca_plot <- function(exprs, pdata, color_var, var_label, title_prefix) {
  exprs_clean <- na.omit(exprs)
  pca <- prcomp(t(exprs_clean), center = TRUE, scale. = FALSE)

  pc_data <- as.data.frame(pca$x[, c(1, 2)])
  colnames(pc_data) <- c("PC1", "PC2")
  pc_data$group <- as.factor(pdata[[color_var]])

  var_explained <- summary(pca)$importance[2, 1:2] * 100

  ggplot(pc_data, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(aes(group = group), type = "norm", level = 0.95, linetype = 2) +
    labs(
      title = paste(title_prefix, "- by", var_label),
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", var_explained[2]),
      color = var_label
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11)
    )
}

plot_a <- create_pca_plot(exprs_before, pdata, "secondaryaccession", "Dataset", "Before ComBat")
plot_b <- create_pca_plot(exprs_after, pdata, "secondaryaccession", "Dataset", "After ComBat")
plot_c <- create_pca_plot(exprs_before, pdata, "Gestational.Age.Category", "Trimester", "Before ComBat")
plot_d <- create_pca_plot(exprs_after, pdata, "Gestational.Age.Category", "Trimester", "After ComBat")

add_panel_label <- function(plot, label) {
  plot + annotation_custom(
    grob = textGrob(label, x = 0.05, y = 0.95,
                    gp = gpar(fontsize = 16, fontface = "bold")),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  )
}

plot_a <- add_panel_label(plot_a, "A")
plot_b <- add_panel_label(plot_b, "B")
plot_c <- add_panel_label(plot_c, "C")
plot_d <- add_panel_label(plot_d, "D")

width_inches <- 8.27
height_inches <- 5.85

output_png <- file.path(output_dir, "pcas", "combined_pca_figure.png")
png(output_png, width = width_inches, height = height_inches, units = "in", res = 300)
grid.arrange(plot_a, plot_b, plot_c, plot_d, ncol = 2, nrow = 2)
dev.off()

output_pdf <- file.path(output_dir, "pcas", "combined_pca_figure.pdf")
pdf(output_pdf, width = width_inches, height = height_inches)
grid.arrange(plot_a, plot_b, plot_c, plot_d, ncol = 2, nrow = 2)
dev.off()

output_svg <- file.path(output_dir, "pcas", "combined_pca_figure.svg")
svg(output_svg, width = width_inches, height = height_inches)
grid.arrange(plot_a, plot_b, plot_c, plot_d, ncol = 2, nrow = 2)
dev.off()

cat("Layout:\n")
cat("  A (top-left):     Before ComBat - Dataset\n")
cat("  B (top-right):    After ComBat - Dataset\n")
cat("  C (bottom-left):  Before ComBat - Trimester\n")
cat("  D (bottom-right): After ComBat - Trimester\n\n")
cat("Datasets:", paste(target_datasets, collapse = ", "), "\n")
cat("Output:", file.path(output_dir, "pcas"), "\n")
cat("\n=== DONE ===\n")
