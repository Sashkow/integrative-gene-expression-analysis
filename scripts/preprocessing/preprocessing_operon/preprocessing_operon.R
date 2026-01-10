# Preprocessing script for Operon arrays (GSE28551)
# This script downloads, processes, and maps GSE28551 data to generate
# expression matrix and phenodata files in the format used by IGEA pipeline

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
pkgs <- c("GEOquery", "Biobase", "org.Hs.eg.db", "AnnotationDbi")
BiocManager::install(pkgs, ask=FALSE, update=FALSE)

library(GEOquery)
library(Biobase)
library(org.Hs.eg.db)
library(AnnotationDbi)

setwd('/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis')

# Define paths
mappedpath = 'data/preprocessed_operon'
pdatapath = 'data/preprocessed_operon'

# Create directories if they don't exist
dir.create(mappedpath, showWarnings = FALSE, recursive = TRUE)
dir.create(pdatapath, showWarnings = FALSE, recursive = TRUE)

# Load GSE28551 data and GPL platform annotation
cat("Loading GSE28551 data from file...\n")
if (file.exists("GSE28551_series_matrix.txt.gz")) {
  # Load expression data without GPL first
  gse <- getGEO(filename = "GSE28551_series_matrix.txt.gz", getGPL = FALSE)
  gset <- gse

  # Load GPL platform annotation separately
  if (file.exists("GPL2986_family.soft.gz")) {
    cat("Loading GPL2986 platform annotation...\n")
    gpl <- getGEO(filename = "GPL2986_family.soft.gz")
    # Add GPL annotation to expression set
    fData(gset) <- Table(gpl)[match(rownames(exprs(gset)), Table(gpl)$ID), ]
  } else {
    cat("Downloading GPL2986 platform annotation...\n")
    gpl <- getGEO("GPL2986")
    fData(gset) <- Table(gpl)[match(rownames(exprs(gset)), Table(gpl)$ID), ]
  }
} else {
  cat("Downloading GSE28551 data from GEO...\n")
  options(timeout = 300)
  gse <- getGEO("GSE28551", GSEMatrix=TRUE, destdir = ".", getGPL = TRUE)
  gset <- gse[[1]]
}

# Extract expression data (already normalized according to GEO)
exprs_data <- exprs(gset)
cat("Downloaded expression data with", nrow(exprs_data), "probes and", ncol(exprs_data), "samples\n")

# Extract phenotype data
pheno <- pData(gset)
cat("Extracted phenotype data for", nrow(pheno), "samples\n")

# Check if data is log2 transformed
qx <- as.numeric(quantile(exprs_data, c(0., 0.25, 0.5, 0.75, 0.99, 1), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (!LogC) {
  cat("Data appears to be already log2 transformed\n")
} else {
  cat("Warning: Data may not be log2 transformed\n")
}

# Map probes to Entrez Gene IDs using GPL platform annotation
cat("\nMapping probes to Entrez Gene IDs...\n")

# Get annotation from the GEO object
feature_data <- fData(gset)
probe_ids <- rownames(exprs_data)

cat("Total probes:", length(probe_ids), "\n")
cat("Available annotation columns:\n")
print(colnames(feature_data))

# Use the GENE column from GPL platform for Entrez Gene IDs
if ("GENE" %in% colnames(feature_data)) {
  entrez_ids <- feature_data[["GENE"]]
  cat("Using 'GENE' column for Entrez ID mapping\n")
} else if ("Gene ID" %in% colnames(feature_data)) {
  entrez_ids <- feature_data[["Gene ID"]]
  cat("Using 'Gene ID' column for Entrez ID mapping\n")
} else {
  stop("No GENE or Gene ID column found in platform annotation")
}

# Create annotation data frame
annotation_df <- data.frame(
  PROBEID = probe_ids,
  ENTREZID = as.character(entrez_ids),
  stringsAsFactors = FALSE
)

cat("Before filtering:", nrow(annotation_df), "probes\n")

# Remove probes with missing, null, or empty Entrez IDs
annotation_df <- annotation_df[!is.na(annotation_df$ENTREZID) &
                                annotation_df$ENTREZID != "" &
                                annotation_df$ENTREZID != "null", ]

cat("After filtering null/empty:", nrow(annotation_df), "probes\n")

# Filter expression data to keep only mapped probes
exprs_mapped <- exprs_data[annotation_df$PROBEID, ]

# Handle multiple probes mapping to the same gene
# Take the probe with highest mean expression for each gene
cat("Handling multiple probes per gene...\n")
probe_means <- rowMeans(exprs_mapped, na.rm = TRUE)
annotation_df$MEAN_EXPR <- probe_means[annotation_df$PROBEID]

# Sort by ENTREZID and MEAN_EXPR (descending)
annotation_df <- annotation_df[order(annotation_df$ENTREZID, -annotation_df$MEAN_EXPR), ]

# Keep only the first (highest expression) probe for each gene
annotation_df <- annotation_df[!duplicated(annotation_df$ENTREZID), ]
cat("After collapsing to unique genes:", nrow(annotation_df), "genes\n")

# Create final expression matrix with Entrez IDs as rownames
exprs_final <- exprs_mapped[annotation_df$PROBEID, ]
rownames(exprs_final) <- annotation_df$ENTREZID

# Save mapped expression data
cat("\nSaving mapped expression data...\n")
write.table(exprs_final,
            paste(mappedpath, "/E-GEOD-28551_mapped_affymetrix.tsv", sep=""),
            sep="\t", quote=FALSE)

# Create and save phenodata file compatible with IGEA format (samples.csv)
cat("Creating phenodata file in IGEA samples.csv format...\n")

# Map trimester information
trimester_map <- pheno$`trimester:ch1`
gestational_category <- ifelse(trimester_map == "first", "First Trimester",
                               ifelse(trimester_map == "third", "Third Trimester", "_"))

# Create phenodata in the same format as samples.csv
# Note: Using dots in column names to match existing samples.csv format
pdata_final <- data.frame(
  arraydatafile_exprscolumnnames = rownames(pheno),
  accession = "E-GEOD-28551",
  secondaryaccession = "GSE28551",
  name = "Operon array expression profiling of first and third trimester placenta",
  releasedate = "_",
  lastupdatedate = "_",
  samples = nrow(pheno),
  microarrays = "Operon Human Genome Array v3.0",
  status = "_",
  Experiment = "E-GEOD-28551",
  Sample.Name = pheno$title,
  exprs_column_names = rownames(pheno),
  platform_batch = "_",
  medical_sample_name = "_",
  Diagnosis = "Healthy",
  Gestational.Age.Category = gestational_category,
  Gestational.Age = "_",
  Average.Gestational.Age = "_",
  Deviation.Gestational.Age = "_",
  Biological.Specimen = "Placenta",
  Maternal.Age = "_",
  Fetus.Sex = "_",
  Gestation = "_",
  Cells..Cultured = "_",
  Average.Maternal.Age = "_",
  Gravidity = "_",
  Pregnancy.Trimesters = trimester_map,
  Scan.Time = "_",
  Diastolic.Pressure = "_",
  Gestational.Age.Upper.Bound = "_",
  Average.Fetal.Weight = "_",
  Systolic.Pressure.at...20.week = "_",
  Sex = "_",
  Parity = "_",
  Diastolic.Pressure.at.delivery = "_",
  Diastolic.Pressure.at...20.week = "_",
  overexpressed.protein = "_",
  Fetus = "_",
  X.excluded. = "_",
  Ethnology = "_",
  Fetal.Weight = "_",
  Arterial.Pressure = "_",
  Labor.Induction = "_",
  Batch = "_",
  Caesarean.Section = "_",
  Organism = "Humans",
  Urine.Protein.Concentration = "_",
  Common = "_",
  Biosource.Provider = "_",
  Gestational.Age.Lower.Bound = "_",
  Pathology = "_",
  Continental.Population.Groups = "_",
  Pre.Eclampsia.Onset = "_",
  Gestational.Age.at.Time.of.Sampling = "_",
  Placenta.Weight = "_",
  Array.Data.File = paste(rownames(pheno), ".txt", sep=""),
  Estimated.Fetus.Sex = "_",
  Scan.Date = "_",
  Transfection = "_",
  culture.time = "_",
  Systolic.Pressure = "_",
  Smoking.Status = "_",
  massir_pattern = "_",
  estimated_sex = "_",
  etimated_sex_unimodality = "_",
  etimated_sex_unimodality_p = "_",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Save phenodata as CSV
write.csv(pdata_final,
          paste(pdatapath, "/GSE28551_pdata.csv", sep=""),
          row.names = FALSE, quote = TRUE)

cat("\n=== Processing complete ===\n")
cat("Expression matrix:", paste(mappedpath, "/E-GEOD-28551_mapped_affymetrix.tsv", sep=""), "\n")
cat("  - Dimensions:", nrow(exprs_final), "genes x", ncol(exprs_final), "samples\n")
cat("Phenodata:", paste(pdatapath, "/GSE28551_pdata.csv", sep=""), "\n")
cat("  - Number of samples:", nrow(pdata_final), "\n")
cat("\nFiles are ready for use in IGEA pipeline!\n")
