if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
for (pkg in c("quantro", "affy")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

library(quantro)
library(affy)

BASE <- "/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis"
pdata <- read.csv(file.path(BASE, "data/phenodata/samples.csv"))

output_dir <- file.path(BASE, "output/single_dataset_analysis/quantro")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

N_PERM <- 1000
results <- list()

run_quantro_safe <- function(mat, groups, dataset, platform, contrast) {
  groups <- factor(groups)
  if (nlevels(groups) < 2) {
    message("  Skipping: only one group level")
    return(NULL)
  }
  grp_tab <- table(groups)
  if (any(grp_tab < 2)) {
    message("  Skipping: group with fewer than 2 samples (",
            paste(names(grp_tab), grp_tab, sep="=", collapse=", "), ")")
    return(NULL)
  }
  message("  Groups: ", paste(names(grp_tab), grp_tab, sep="=", collapse=", "))
  qt <- quantro(mat, groupFactor = groups, B = N_PERM)
  message("  quantro stat: ", round(quantroStat(qt), 4),
          "  p-value: ", round(quantroPvalPerm(qt), 4))
  data.frame(
    Dataset = dataset, Platform = platform, Contrast = contrast,
    Groups = paste(names(grp_tab), collapse=" / "),
    N_per_group = paste(grp_tab, collapse="/"),
    quantro_stat = round(quantroStat(qt), 4),
    perm_pvalue = round(quantroPvalPerm(qt), 4)
  )
}

get_gsm_ids <- function(series_matrix_path) {
  con <- gzfile(series_matrix_path)
  lines <- readLines(con)
  close(con)
  geo_line <- lines[grepl("^!Sample_geo_accession", lines)]
  gsub("\"", "", strsplit(geo_line, "\t")[[1]][-1])
}

valid_cats <- c("First Trimester", "Second Trimester", "Term")

# ============================================================
# GSE100051 (Illumina): non-normalized bead-level signal
# 1st=42, 2nd=7, Term=5
# ============================================================
message("\n=== GSE100051 (Illumina HumanHT-12 V4.0) ===")
message("  Loading non-normalized signal intensities...")
raw100051 <- read.delim(
  file.path(BASE, "data/raws/GSE100051/GSE100051_non-normalized.txt"),
  skip = 1, row.names = 1
)
signal_cols <- grep("AVG_Signal", colnames(raw100051))
mat100051 <- log2(pmax(as.matrix(raw100051[, signal_cols]), 1))

gsm_order <- get_gsm_ids(file.path(BASE,
  "data/raws/GSE100051/GSE100051_series_matrix.txt.gz"))
colnames(mat100051) <- gsm_order

sub100051 <- pdata[grepl("GSE100051", pdata$secondaryaccession), ]
sub100051 <- sub100051[sub100051$Gestational.Age.Category %in% valid_cats, ]
idx <- match(sub100051$Sample.Name, colnames(mat100051))
mat100051_sub <- mat100051[, idx[!is.na(idx)]]
sub100051 <- sub100051[!is.na(idx), ]

message("  GA categories test:")
results[[length(results) + 1]] <- run_quantro_safe(
  mat100051_sub, sub100051$Gestational.Age.Category,
  "GSE100051", "Illumina HumanHT-12 V4.0", "GA categories"
)

message("  Male vs Female test:")
sex_keep <- sub100051$Combined.Fetus.Sex %in% c("Male", "Female")
results[[length(results) + 1]] <- run_quantro_safe(
  mat100051_sub[, sex_keep], sub100051$Combined.Fetus.Sex[sex_keep],
  "GSE100051", "Illumina HumanHT-12 V4.0", "Male vs Female"
)

# ============================================================
# Affy datasets: raw PM probe values via affy::pm()
# Following Hicks & Irizarry 2015 Methods: "we extracted the
# raw perfect match (PM) values from the CEL files using the
# affy R/Bioconductor package"
# ============================================================
load_affy_pm <- function(cel_dir) {
  raw_data <- ReadAffy(celfile.path = cel_dir)
  pm_mat <- pm(raw_data)
  log2(pm_mat)
}

# GSE9984 (Affy): 1st=4, 2nd=4, Term=4
message("\n=== GSE9984 (Affy HG-U133 Plus 2.0) ===")
message("  Extracting raw PM values from CEL files...")
mat9984 <- load_affy_pm(file.path(BASE, "data/raws/GSE9984/cel"))
colnames(mat9984) <- gsub("\\.CEL$", "", colnames(mat9984), ignore.case = TRUE)

sub9984 <- pdata[grepl("GSE9984", pdata$secondaryaccession), ]
sub9984 <- sub9984[sub9984$Gestational.Age.Category %in% valid_cats, ]
sample_ids <- gsub("\\.CEL$", "", sub9984$arraydatafile_exprscolumnnames,
                   ignore.case = TRUE)
idx <- match(sample_ids, colnames(mat9984))
mat9984_sub <- mat9984[, idx[!is.na(idx)]]
sub9984 <- sub9984[!is.na(idx), ]

message("  GA categories test:")
results[[length(results) + 1]] <- run_quantro_safe(
  mat9984_sub, sub9984$Gestational.Age.Category,
  "GSE9984", "Affy HG-U133 Plus 2.0", "GA categories"
)

message("  Male vs Female test:")
sex_keep <- sub9984$Combined.Fetus.Sex %in% c("Male", "Female")
results[[length(results) + 1]] <- run_quantro_safe(
  mat9984_sub[, sex_keep], sub9984$Combined.Fetus.Sex[sex_keep],
  "GSE9984", "Affy HG-U133 Plus 2.0", "Male vs Female"
)

# GSE22490 (Affy): 1st=8, 2nd=2
message("\n=== GSE22490 (Affy HG-U133 Plus 2.0) ===")
message("  Extracting raw PM values from CEL files...")
mat22490 <- load_affy_pm(file.path(BASE, "data/raws/GSE22490/cel"))
colnames(mat22490) <- gsub("\\.CEL$", "", colnames(mat22490), ignore.case = TRUE)

sub22490 <- pdata[grepl("GSE22490", pdata$secondaryaccession), ]
sub22490 <- sub22490[sub22490$Gestational.Age.Category %in% valid_cats, ]
sample_ids <- gsub("\\.CEL$", "", sub22490$arraydatafile_exprscolumnnames,
                   ignore.case = TRUE)
idx <- match(sample_ids, colnames(mat22490))
mat22490_sub <- mat22490[, idx[!is.na(idx)]]
sub22490 <- sub22490[!is.na(idx), ]

message("  GA categories test:")
results[[length(results) + 1]] <- run_quantro_safe(
  mat22490_sub, sub22490$Gestational.Age.Category,
  "GSE22490", "Affy HG-U133 Plus 2.0", "GA categories"
)

# ============================================================
# GSE28551 (ABI): raw Signal from individual probe-level files
# 1st=16, Term=21
# ============================================================
message("\n=== GSE28551 (ABI Human Genome Survey v2) ===")
message("  Loading raw probe-level signal from individual files...")
raw_dir_28551 <- file.path(BASE, "data/raws/GSE28551/raw")
raw_files_28551 <- list.files(raw_dir_28551, pattern = "\\.txt\\.gz$",
                               full.names = TRUE)

mat_list_28551 <- list()
for (f in raw_files_28551) {
  gsm <- gsub("\\.txt\\.gz$", "", basename(f))
  d <- read.delim(gzfile(f))
  mat_list_28551[[gsm]] <- setNames(d$Signal, d$Probe_ID)
}
all_probes <- Reduce(intersect, lapply(mat_list_28551, names))
mat28551 <- sapply(mat_list_28551, function(x) x[all_probes])
mat28551 <- log2(pmax(mat28551, 1))

sub28551 <- pdata[grepl("GSE28551", pdata$secondaryaccession), ]
sub28551_ga <- sub28551[sub28551$Gestational.Age.Category %in% valid_cats, ]
idx <- match(sub28551_ga$arraydatafile_exprscolumnnames, colnames(mat28551))
mat28551_ga <- mat28551[, idx[!is.na(idx)]]
sub28551_ga <- sub28551_ga[!is.na(idx), ]

message("  GA categories test (1st vs Term):")
results[[length(results) + 1]] <- run_quantro_safe(
  mat28551_ga, sub28551_ga$Gestational.Age.Category,
  "GSE28551", "ABI Human Genome Survey v2", "GA categories"
)

message("  Male vs Female test (1st trimester only):")
sub28551_sex <- sub28551_ga[sub28551_ga$Gestational.Age.Category == "First Trimester" &
                              sub28551_ga$Combined.Fetus.Sex %in% c("Male","Female"), ]
idx_sex <- match(sub28551_sex$arraydatafile_exprscolumnnames, colnames(mat28551_ga))
if (sum(!is.na(idx_sex)) >= 4) {
  results[[length(results) + 1]] <- run_quantro_safe(
    mat28551_ga[, idx_sex[!is.na(idx_sex)]],
    sub28551_sex$Combined.Fetus.Sex[!is.na(idx_sex)],
    "GSE28551", "ABI Human Genome Survey v2", "Male vs Female (1st tri)"
  )
}

# ============================================================
# GSE93520 (Agilent): raw gMedianSignal from feature extraction files
# 1st only, F:19 M:17
# ============================================================
message("\n=== GSE93520 (Agilent 4x44K) ===")
message("  Loading raw gMedianSignal from feature extraction files...")
raw_dir_93520 <- file.path(BASE, "data/raws/GSE93520/raw")
raw_files_93520 <- list.files(raw_dir_93520, pattern = "\\.txt\\.gz$",
                               full.names = TRUE)

mat_list_93520 <- list()
for (f in raw_files_93520) {
  gsm <- sub("_.*", "", basename(f))
  lines <- readLines(gzfile(f))
  feat_start <- which(grepl("^FEATURES", lines))
  header <- strsplit(lines[feat_start], "\t")[[1]]
  data_lines <- lines[(feat_start + 1):length(lines)]
  data_lines <- data_lines[grepl("^DATA", data_lines)]
  d <- read.delim(textConnection(data_lines), header = FALSE)
  colnames(d) <- header
  probe_col <- which(colnames(d) == "ProbeName")
  signal_col <- which(colnames(d) == "gMedianSignal")
  mat_list_93520[[gsm]] <- setNames(as.numeric(d[[signal_col]]),
                                     d[[probe_col]])
}
all_probes_93520 <- Reduce(intersect, lapply(mat_list_93520, names))
mat93520 <- sapply(mat_list_93520, function(x) x[all_probes_93520])
mat93520 <- log2(pmax(mat93520, 1))

sub93520 <- pdata[grepl("GSE93520", pdata$secondaryaccession), ]
sub93520 <- sub93520[sub93520$Combined.Fetus.Sex %in% c("Male", "Female"), ]
idx <- match(sub93520$Sample.Name, colnames(mat93520))
mat93520_sub <- mat93520[, idx[!is.na(idx)]]
sub93520 <- sub93520[!is.na(idx), ]

message("  Male vs Female test:")
results[[length(results) + 1]] <- run_quantro_safe(
  mat93520_sub, sub93520$Combined.Fetus.Sex,
  "GSE93520", "Agilent 4x44K", "Male vs Female"
)

# ============================================================
# Combined results
# ============================================================
results <- results[!sapply(results, is.null)]
results_df <- do.call(rbind, results)

message("\n=== All Results ===")
print(results_df, row.names = FALSE)
write.csv(results_df, file.path(output_dir, "quantro_results.csv"),
          row.names = FALSE)

message("\nInterpretation:")
message("  p < 0.05 => global distribution differences exist between groups")
message("           => quantile normalization may remove biological signal")
message("  p > 0.05 => no evidence of global shift")
message("           => quantile normalization assumption is reasonable")
message("\nMethod follows Hicks & Irizarry (2015) Genome Biology 16:117")
message("  Affy: raw PM values via affy::pm(), log2-transformed")
message("  Illumina: non-normalized AVG_Signal, log2-transformed")
message("  ABI: raw Signal from probe-level files, log2-transformed")
message("  Agilent: raw gMedianSignal from feature extraction, log2-transformed")
