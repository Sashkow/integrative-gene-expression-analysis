#!/usr/bin/env Rscript

#' Estimate Fetus Sex using massiR package
#'
#' This script estimates fetus sex for samples in a dataset using Y chromosome
#' gene expression patterns via the massiR package.
#'
#' @usage Rscript estimate_fetus_sex.R <dataset_id>
#' @example Rscript estimate_fetus_sex.R GSE93520
#'
#' @author Expression Integration Pipeline
#' @date 2025-10-26

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  cat("Usage: Rscript estimate_fetus_sex.R <dataset_id>\n")
  cat("Example: Rscript estimate_fetus_sex.R GSE93520\n\n")
  stop("Please provide a dataset ID")
}

dataset_id <- args[1]

cat("\n=== Estimating Fetus Sex for", dataset_id, "===\n\n")

# Load required libraries
suppressPackageStartupMessages({
  library(massiR)
  library(biomaRt)
})

# Define paths
phenodata_file <- "data/phenodata/samples.csv"
mapped_dir <- "data/mapped"
output_dir <- "output/sex_estimation"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ===================================================================
# STEP 1: LOAD PHENODATA
# ===================================================================

cat("Loading phenodata...\n")
phenodata <- read.csv(phenodata_file, stringsAsFactors = FALSE)

# Filter for the specified dataset
dataset_samples <- phenodata[phenodata$secondaryaccession == dataset_id, ]

if (nrow(dataset_samples) == 0) {
  stop("No samples found for dataset: ", dataset_id)
}

cat("  Found", nrow(dataset_samples), "samples for", dataset_id, "\n")
cat("  Platform:", unique(dataset_samples$microarrays)[1], "\n\n")

# ===================================================================
# STEP 2: LOAD EXPRESSION DATA
# ===================================================================

cat("Loading expression data...\n")

# Try to find expression file
expr_file_patterns <- c(
  paste0(dataset_id, "_matrix_no_filtering.tsv"),
  paste0(dataset_id, "_mapped.tsv"),
  paste0(dataset_id, ".tsv")
)

expr_file <- NULL
for (pattern in expr_file_patterns) {
  candidate <- file.path(mapped_dir, pattern)
  if (file.exists(candidate)) {
    expr_file <- candidate
    break
  }
}

if (is.null(expr_file)) {
  stop("Expression data file not found for ", dataset_id,
       "\nLooked for: ", paste(expr_file_patterns, collapse = ", "))
}

cat("  Loading:", expr_file, "\n")
exprs <- read.table(expr_file, header = TRUE, sep = "\t", row.names = 1,
                    stringsAsFactors = FALSE, check.names = FALSE)

cat("  Expression matrix dimensions:", nrow(exprs), "genes ×", ncol(exprs), "samples\n\n")

# ===================================================================
# STEP 3: SUBSET EXPRESSION DATA FOR DATASET SAMPLES
# ===================================================================

cat("Subsetting expression data for dataset samples...\n")

# Get proper column names (matching arraydatafile_exprscolumnnames)
sample_ids <- dataset_samples$arraydatafile_exprscolumnnames

# Check which samples are present in expression data
samples_in_exprs <- intersect(sample_ids, colnames(exprs))

if (length(samples_in_exprs) == 0) {
  # Try with make.names transformation
  sample_ids_clean <- make.names(sample_ids)
  samples_in_exprs <- intersect(sample_ids_clean, colnames(exprs))

  if (length(samples_in_exprs) == 0) {
    stop("No matching samples found in expression data.\n",
         "Sample IDs from phenodata: ", paste(head(sample_ids, 3), collapse = ", "), "\n",
         "Column names in exprs: ", paste(head(colnames(exprs), 3), collapse = ", "))
  }

  cat("  Using make.names() to match sample IDs\n")
  sample_ids <- sample_ids_clean
}

cat("  Found", length(samples_in_exprs), "matching samples\n")

# Subset expression data
exprs_subset <- exprs[, samples_in_exprs, drop = FALSE]
cat("  Subset expression matrix:", nrow(exprs_subset), "genes ×",
    ncol(exprs_subset), "samples\n\n")

# ===================================================================
# STEP 4: GET Y CHROMOSOME GENES
# ===================================================================

cat("Retrieving Y chromosome genes from Ensembl...\n")

# Get platform identifier for biomaRt
# Extract platform abbreviation from microarray name
platform_name <- unique(dataset_samples$microarrays)[1]
cat("  Platform:", platform_name, "\n")

# Determine massir pattern from platform
# This maps platform names to biomaRt filter names
massir_pattern <- if (grepl("agilent.*4x44k", tolower(platform_name))) {
  "agilent_wholegenome_4x44k_v1"
} else if (grepl("agilent.*4x44k.*v2", tolower(platform_name))) {
  "agilent_wholegenome_4x44k_v2"
} else if (grepl("illumina.*wg.*v2", tolower(platform_name))) {
  "illumina_humanwg_6_v2"
} else if (grepl("illumina.*wg.*v3", tolower(platform_name))) {
  "illumina_humanwg_6_v3"
} else if (grepl("illumina.*ht.*v4", tolower(platform_name))) {
  "illumina_humanht_12_v4"
} else if (grepl("hg.*u133.*plus.*2", tolower(platform_name))) {
  "affy_hg_u133_plus_2"
} else if (grepl("hugene.*1.*0.*st", tolower(platform_name))) {
  "affy_hugene_1_0_st_v1"
} else {
  # Default: try to use entrezgene mapping
  cat("  Warning: Could not determine exact platform filter, using entrezgene_id\n")
  "entrezgene_id"
}

cat("  Using biomaRt filter:", massir_pattern, "\n")

# Connect to Ensembl
tryCatch({
  mart <- useMart('ensembl', dataset = "hsapiens_gene_ensembl")

  # Get Y chromosome genes
  if (massir_pattern == "entrezgene_id") {
    # Use entrezgene directly if we couldn't determine platform
    gene.attributes <- getBM(
      mart = mart,
      values = TRUE,
      attributes = c("entrezgene_id", "chromosome_name")
    )

    y_genes <- gene.attributes[which(!is.na(gene.attributes$entrezgene_id) &
                                     gene.attributes$chromosome_name == "Y"), ]$entrezgene_id
  } else {
    # Use platform-specific filter
    filter_name <- paste0("with_", massir_pattern)

    gene.attributes <- getBM(
      mart = mart,
      values = TRUE,
      filters = c(filter_name),
      attributes = c(massir_pattern, "entrezgene_id", "chromosome_name")
    )

    y_genes <- gene.attributes[which(!is.na(gene.attributes$entrezgene_id) &
                                     gene.attributes$chromosome_name == "Y"), ]$entrezgene_id
  }

  y_genes <- unique(y_genes)
  cat("  Found", length(y_genes), "Y chromosome genes in Ensembl\n")

}, error = function(e) {
  cat("  Error connecting to Ensembl:", conditionMessage(e), "\n")
  cat("  Using fallback Y chromosome genes (Entrez IDs)\n")

  # Fallback: use Y chromosome gene Entrez IDs from literature
  # These are common Y chromosome genes used for sex determination
  # From massiR paper and Y chromosome literature
  y_genes <<- c(
    "8653",   # DDX3Y
    "9086",   # EIF1AY
    "8284",   # KDM5D
    "22829",  # NLGN4Y
    "5539",   # PRKY
    "6192",   # RPS4Y1
    "9087",   # TMSB4Y
    "246126", # TXLNGY
    "8287",   # USP9Y
    "7404",   # UTY
    "7543",   # ZFY
    "83869",  # TTTY14
    "83866",  # TTTY15
    "83867",  # TTTY19
    "266"     # AMELY
  )
  cat("  Using", length(y_genes), "known Y chromosome genes\n")
})

# Find Y genes present in expression data
y_genes_in_data <- intersect(as.character(y_genes), rownames(exprs_subset))

if (length(y_genes_in_data) == 0) {
  stop("No Y chromosome genes found in expression data.\n",
       "Genes searched: ", paste(head(y_genes, 10), collapse = ", "))
}

cat("  Found", length(y_genes_in_data), "Y chromosome genes in expression data\n\n")

# ===================================================================
# STEP 5: RUN MASSIR SEX ESTIMATION
# ===================================================================

cat("Running massiR sex estimation...\n")

# Create probe data frame for massiR
genes.df <- data.frame(matrix(nrow = length(y_genes_in_data), ncol = 0))
rownames(genes.df) <- y_genes_in_data

# Convert expression data to data frame (massiR requirement)
exprs_subset <- as.data.frame(exprs_subset)

# Step 1: Calculate coefficient of variation for Y genes
cat("  Step 1: Calculating Y gene coefficient of variation...\n")
massi.y.out <- massi_y(exprs_subset, genes.df)

# Step 2: Select informative probes (threshold = 4 means top 75% by CV)
cat("  Step 2: Selecting informative Y genes (threshold = 4)...\n")
massi.select.out <- massi_select(exprs_subset, genes.df, threshold = 4)
cat("    Selected", nrow(massi.select.out), "informative genes\n")

if (nrow(massi.select.out) < 2) {
  stop("Not enough informative Y genes (need at least 2, found ",
       nrow(massi.select.out), ")")
}

# Step 3: Test for bimodality (unimodality = bad, bimodal = good)
cat("  Step 3: Testing for bimodality...\n")
massi.dip.out <- massi_dip(massi.select.out)
dip_stat <- massi.dip.out$dip.statistics$statistic
dip_pval <- massi.dip.out$dip.statistics$p.value

cat("    Dip statistic:", round(dip_stat, 4), "(p =", round(dip_pval, 4), ")\n")
cat("    ", if (dip_stat > 0.05) "Good bimodal distribution" else
        "Warning: Distribution may not be clearly bimodal", "\n")

# Step 4: Cluster samples into male/female
cat("  Step 4: Clustering samples by sex...\n")
massi.cluster.out <- massi_cluster(massi.select.out)

# Extract results
sex_results <- massi.cluster.out$massi.results

cat("    Predicted sex distribution:\n")
sex_table <- table(sex_results$sex)
for (sex in names(sex_table)) {
  cat("      ", sex, ":", sex_table[sex], "\n")
}
cat("\n")

# ===================================================================
# STEP 6: UPDATE PHENODATA WITH SEX ESTIMATES
# ===================================================================

cat("Updating phenodata with sex estimates...\n")

# Create a mapping data frame
sex_mapping <- data.frame(
  arraydatafile_exprscolumnnames = sex_results$ID,
  Estimated.Fetus.Sex = ifelse(sex_results$sex == "male", "Male",
                               ifelse(sex_results$sex == "female", "Female", NA)),
  estimated_sex = ifelse(sex_results$sex == "male", "Male",
                        ifelse(sex_results$sex == "female", "Female", NA)),
  etimated_sex_unimodality = dip_stat,
  etimated_sex_unimodality_p = dip_pval,
  stringsAsFactors = FALSE
)

# If we used make.names, we need to map back to original sample IDs
if (all(sample_ids == make.names(dataset_samples$arraydatafile_exprscolumnnames))) {
  # Create reverse mapping
  original_to_clean <- data.frame(
    original = dataset_samples$arraydatafile_exprscolumnnames,
    clean = sample_ids,
    stringsAsFactors = FALSE
  )

  sex_mapping <- merge(
    data.frame(clean = sex_mapping$arraydatafile_exprscolumnnames,
               Estimated.Fetus.Sex = sex_mapping$Estimated.Fetus.Sex,
               estimated_sex = sex_mapping$estimated_sex,
               stringsAsFactors = FALSE),
    original_to_clean,
    by = "clean"
  )

  sex_mapping$arraydatafile_exprscolumnnames <- sex_mapping$original
  sex_mapping <- sex_mapping[, c("arraydatafile_exprscolumnnames",
                                  "Estimated.Fetus.Sex", "estimated_sex")]
  sex_mapping$etimated_sex_unimodality <- dip_stat
  sex_mapping$etimated_sex_unimodality_p <- dip_pval
}

# Read current phenodata
phenodata_current <- read.csv(phenodata_file, stringsAsFactors = FALSE)

# Update Estimated.Fetus.Sex for this dataset's samples
for (i in 1:nrow(sex_mapping)) {
  sample_id <- sex_mapping$arraydatafile_exprscolumnnames[i]
  idx <- which(phenodata_current$arraydatafile_exprscolumnnames == sample_id)

  if (length(idx) > 0) {
    phenodata_current$Estimated.Fetus.Sex[idx] <- sex_mapping$Estimated.Fetus.Sex[i]
    phenodata_current$estimated_sex[idx] <- sex_mapping$estimated_sex[i]
    phenodata_current$etimated_sex_unimodality[idx] <- sex_mapping$etimated_sex_unimodality[i]
    phenodata_current$etimated_sex_unimodality_p[idx] <- sex_mapping$etimated_sex_unimodality_p[i]
  }
}

# Save updated phenodata
phenodata_output <- phenodata_file
write.csv(phenodata_current, phenodata_output, row.names = FALSE, quote = TRUE)

cat("  ✓ Updated phenodata saved to:", phenodata_output, "\n\n")

# ===================================================================
# STEP 7: SAVE SEX ESTIMATION RESULTS
# ===================================================================

cat("Saving sex estimation results...\n")

# Save detailed results
results_file <- file.path(output_dir, paste0(dataset_id, "_sex_estimation.csv"))
write.csv(sex_mapping, results_file, row.names = FALSE, quote = TRUE)
cat("  ✓ Sex estimation results:", results_file, "\n")

# Save summary
summary_file <- file.path(output_dir, paste0(dataset_id, "_sex_summary.txt"))
sink(summary_file)
cat("Sex Estimation Summary for", dataset_id, "\n")
cat(rep("=", 60), "\n\n", sep = "")
cat("Dataset:", dataset_id, "\n")
cat("Platform:", platform_name, "\n")
cat("Total samples:", nrow(dataset_samples), "\n")
cat("Samples with sex estimated:", nrow(sex_mapping), "\n\n")
cat("Y chromosome genes found:", length(y_genes_in_data), "\n")
cat("Informative genes used:", nrow(massi.select.out), "\n\n")
cat("Bimodality test:\n")
cat("  Dip statistic:", round(dip_stat, 4), "\n")
cat("  P-value:", round(dip_pval, 4), "\n")
cat("  Interpretation:", if (dip_stat > 0.05) "Good bimodal distribution\n" else
                         "Warning: May not be clearly bimodal\n")
cat("\nSex distribution:\n")
print(table(sex_mapping$Estimated.Fetus.Sex))
cat("\n")
sink()
cat("  ✓ Summary saved to:", summary_file, "\n\n")

# Generate plots if there are both males and females
if (length(unique(sex_results$sex)) >= 2) {
  cat("Generating diagnostic plots...\n")

  # Save cluster plot
  plot_file <- file.path(output_dir, paste0(dataset_id, "_sex_cluster_plot.png"))
  png(plot_file, width = 10, height = 8, units = "in", res = 300)
  massi_cluster_plot(massi.select.out, massi.cluster.out)
  dev.off()
  cat("  ✓ Cluster plot saved to:", plot_file, "\n\n")
}

# ===================================================================
# SUMMARY
# ===================================================================

cat("=================================================================\n")
cat("                  SEX ESTIMATION COMPLETE                        \n")
cat("=================================================================\n\n")

cat("Dataset:", dataset_id, "\n")
cat("Samples processed:", nrow(sex_mapping), "\n")
cat("Sex distribution:\n")
print(table(sex_mapping$Estimated.Fetus.Sex))
cat("\nPhenodata updated:", phenodata_output, "\n")
cat("Results directory:", output_dir, "\n\n")
