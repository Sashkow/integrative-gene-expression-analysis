#!/usr/bin/env Rscript
# Frequency polygon of gestational age distribution by tissue type.
# All healthy, non-excluded samples. No dataset distinction.

library(ggplot2)
source("scripts/visualization/ga_helpers.R")

plot_ga_density_by_tissue <- function(
    output_dir = "output/visualization/gestational_age_estimation",
    width = 12, height = 7, dpi = 300) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  samples_csv <- read.csv("data/phenodata/samples.csv", stringsAsFactors = FALSE)

  pheno <- samples_csv[
    tolower(samples_csv$Diagnosis) == "healthy" &
    !(samples_csv$X.excluded. %in% TRUE), ]

  if (nrow(pheno) == 0) {
    cat("No matching samples found.\n")
    return(invisible(NULL))
  }

  df <- parse_ga_rows(pheno, extra_cols = c(tissue = "Biological.Specimen"))
  df <- collapse_rep_blocks(df)
  df <- fill_range_from_category(df)
  df <- estimate_weeks(df)

  df$tissue[df$tissue == "_"] <- NA
  df <- df[!is.na(df$tissue), ]

  tissue_counts <- sort(table(df$tissue), decreasing = TRUE)
  cat("\nSamples per tissue:\n")
  print(tissue_counts)

  tissue_labels <- sprintf("%s (n=%d)", names(tissue_counts), as.integer(tissue_counts))
  df$tissue_label <- sprintf("%s (n=%d)", df$tissue, as.integer(tissue_counts[df$tissue]))
  df$tissue_label <- factor(df$tissue_label, levels = tissue_labels)

  tissue_colors <- c(
    "#E69F00", "#56B4E9", "#009E73", "#D55E00",
    "#CC79A7", "#0072B2", "#F0E442", "#999999",
    "#882255", "#44AA99", "#332288")

  all_weeks <- seq(4, 42)
  grid <- expand.grid(tissue_label = levels(df$tissue_label),
                      estimated_week = all_weeks,
                      stringsAsFactors = FALSE)
  counts <- as.data.frame(table(tissue_label = df$tissue_label,
                                estimated_week = df$estimated_week),
                          stringsAsFactors = FALSE)
  counts$estimated_week <- as.integer(counts$estimated_week)
  grid <- merge(grid, counts, by = c("tissue_label", "estimated_week"), all.x = TRUE)
  grid$Freq[is.na(grid$Freq)] <- 0
  grid$tissue_label <- factor(grid$tissue_label, levels = levels(df$tissue_label))

  grid <- grid[order(grid$tissue_label, grid$estimated_week), ]
  grid$smooth <- ave(grid$Freq, grid$tissue_label, FUN = function(x) {
    stats::filter(x, rep(1/3, 3), sides = 2)
  })
  grid$smooth[is.na(grid$smooth)] <- grid$Freq[is.na(grid$smooth)]

  p <- ggplot(grid, aes(x = estimated_week, y = smooth,
                         color = tissue_label, group = tissue_label)) +
    geom_line(linewidth = 0.8) +
    scale_x_continuous(breaks = seq(4, 42, by = 1)) +
    scale_color_manual(values = tissue_colors[seq_along(tissue_counts)]) +
    geom_vline(xintercept = 12.5, lty = 2, color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 27.5, lty = 2, color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 33.5, lty = 2, color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 36.5, lty = 2, color = "grey60", linewidth = 0.4) +
    annotate("text", x = 8, y = Inf, label = "T1", vjust = 1.5, size = 3, color = "grey50") +
    annotate("text", x = 20, y = Inf, label = "T2", vjust = 1.5, size = 3, color = "grey50") +
    annotate("text", x = 30.5, y = Inf, label = "Early\nPreterm", vjust = 1.3, size = 2.5, color = "grey50") +
    annotate("text", x = 35, y = Inf, label = "Late\nPreterm", vjust = 1.3, size = 2.5, color = "grey50") +
    annotate("text", x = 39, y = Inf, label = "Term", vjust = 1.5, size = 3, color = "grey50") +
    labs(
      title = "Gestational age distribution by tissue type",
      subtitle = sprintf("Healthy samples, all datasets (n=%d)", nrow(df)),
      x = "Gestational week",
      y = "Number of samples",
      color = "Tissue"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      legend.position = "right",
      legend.text = element_text(size = 9),
      axis.text.x = element_text(size = 7)
    )

  plot_path <- file.path(output_dir, "gestational_age_density_by_tissue.png")
  png(plot_path, width = width, height = height, units = "in", res = dpi)
  print(p)
  dev.off()
  cat("\nSaved:", plot_path, "\n")
}

if (!interactive()) {
  plot_ga_density_by_tissue()
}
