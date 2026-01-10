#!/usr/bin/env Rscript

# Analyze PCA positions to describe dataset locations

library(stats)

# Load data
data_dir <- "output/test_first_datasets/GSE55439_GSE93520_GSE28551_GSE100051"
pdata <- read.csv(file.path(data_dir, "phenodata.csv"), stringsAsFactors = FALSE)

# Load expression matrices
exprs_before <- read.table(
  file.path(data_dir, "merged_exprs_before_combat.tsv"),
  header = TRUE, sep = "\t", check.names = FALSE
)
exprs_after <- read.table(
  file.path(data_dir, "merged_exprs_after_combat.tsv"),
  header = TRUE, sep = "\t", check.names = FALSE
)

# Perform PCA
pca_before <- prcomp(t(na.omit(exprs_before)), center = TRUE, scale. = FALSE)
pca_after <- prcomp(t(na.omit(exprs_after)), center = TRUE, scale. = FALSE)

# Get PC scores
pc_before <- data.frame(
  PC1 = pca_before$x[, 1],
  PC2 = pca_before$x[, 2],
  dataset = pdata$secondaryaccession,
  trimester = pdata$Gestational.Age.Category
)

pc_after <- data.frame(
  PC1 = pca_after$x[, 1],
  PC2 = pca_after$x[, 2],
  dataset = pdata$secondaryaccession,
  trimester = pdata$Gestational.Age.Category
)

cat("\n=== BEFORE COMBAT ===\n\n")
cat("Dataset centroids in PC space:\n")
for (ds in unique(pdata$secondaryaccession)) {
  subset_data <- pc_before[pc_before$dataset == ds, ]
  cat(sprintf("%-12s: PC1 mean = %7.1f, PC2 mean = %7.1f, n = %2d\n",
              ds,
              mean(subset_data$PC1),
              mean(subset_data$PC2),
              nrow(subset_data)))
}

cat("\n=== AFTER COMBAT ===\n\n")
cat("Dataset centroids in PC space:\n")
for (ds in unique(pdata$secondaryaccession)) {
  subset_data <- pc_after[pc_after$dataset == ds, ]
  cat(sprintf("%-12s: PC1 mean = %7.1f, PC2 mean = %7.1f, n = %2d\n",
              ds,
              mean(subset_data$PC1),
              mean(subset_data$PC2),
              nrow(subset_data)))
}

cat("\n=== TRIMESTER SEPARATION ===\n\n")
cat("BEFORE ComBat - Trimester centroids:\n")
for (trim in c("First Trimester", "Second Trimester")) {
  subset_data <- pc_before[pc_before$trimester == trim, ]
  cat(sprintf("%-16s: PC1 mean = %7.1f, PC2 mean = %7.1f, n = %2d\n",
              trim,
              mean(subset_data$PC1),
              mean(subset_data$PC2),
              nrow(subset_data)))
}

cat("\nAFTER ComBat - Trimester centroids:\n")
for (trim in c("First Trimester", "Second Trimester")) {
  subset_data <- pc_after[pc_after$trimester == trim, ]
  cat(sprintf("%-16s: PC1 mean = %7.1f, PC2 mean = %7.1f, n = %2d\n",
              trim,
              mean(subset_data$PC1),
              mean(subset_data$PC2),
              nrow(subset_data)))
}

cat("\n")
