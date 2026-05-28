if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
for (pkg in c("quantro", "oligo", "openxlsx")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

library(quantro)

BASE <- dirname(sys.frame(1)$ofile)
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

load_affy_no_qn <- function(cel_dir) {
  cel_files <- list.files(cel_dir, pattern = "\\.CEL$", full.names = TRUE,
                          ignore.case = TRUE)
  if (requireNamespace("oligo", quietly = TRUE)) {
    raw_data <- oligo::read.celfiles(cel_files)
    eset <- oligo::rma(raw_data, normalize = FALSE)
  } else {
    raw_data <- affy::ReadAffy(celfile.path = cel_dir)
    eset <- affy::rma(raw_data, normalize = FALSE)
  }
  exprs(eset)
}

get_gsm_ids <- function(series_matrix_path) {
  con <- gzfile(series_matrix_path)
  lines <- readLines(con)
  close(con)
  geo_line <- lines[grepl("^!Sample_geo_accession", lines)]
  gsm_ids <- gsub("\"", "", strsplit(geo_line, "\t")[[1]][-1])
  gsm_ids
}

valid_cats <- c("First Trimester", "Second Trimester", "Term")

# ============================================================
# GSE100051 (Illumina): 1st=42, 2nd=7, Term=5
# ============================================================
message("\n=== GSE100051 (Illumina HumanHT-12 V4.0) ===")
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
# GSE9984 (Affy): 1st=4, 2nd=4, Term=4
# ============================================================
message("\n=== GSE9984 (Affy HG-U133 Plus 2.0) ===")
mat9984 <- load_affy_no_qn(file.path(BASE, "data/raws/GSE9984_cel"))
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

# ============================================================
# GSE22490 (Affy): 1st=8, 2nd=2
# ============================================================
message("\n=== GSE22490 (Affy HG-U133 Plus 2.0) ===")
mat22490 <- load_affy_no_qn(file.path(BASE, "data/raws/GSE22490_cel"))
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
# GSE122214 (Affy): 1st only, 4 samples
# ============================================================
message("\n=== GSE122214 (Affy HG-U133 Plus 2.0) ===")
mat122214 <- load_affy_no_qn(file.path(BASE, "data/raws/GSE122214/cel"))
colnames(mat122214) <- gsub("\\.CEL$", "", colnames(mat122214), ignore.case = TRUE)

sub122214 <- pdata[grepl("GSE122214", pdata$secondaryaccession), ]
sub122214 <- sub122214[sub122214$Combined.Fetus.Sex %in% c("Male", "Female"), ]
sample_ids <- gsub("\\.CEL$", "", sub122214$arraydatafile_exprscolumnnames,
                   ignore.case = TRUE)
idx <- match(sample_ids, colnames(mat122214))
mat122214_sub <- mat122214[, idx[!is.na(idx)]]
sub122214 <- sub122214[!is.na(idx), ]

message("  Male vs Female test:")
results[[length(results) + 1]] <- run_quantro_safe(
  mat122214_sub, sub122214$Combined.Fetus.Sex,
  "GSE122214", "Affy HG-U133 Plus 2.0", "Male vs Female"
)

# ============================================================
# GSE28551 (ABI): 1st=16, Term=21
# ============================================================
message("\n=== GSE28551 (ABI Human Genome Survey v2) ===")
con28551 <- gzfile(file.path(BASE,
  "data/raws/GSE28551/GSE28551_series_matrix.txt.gz"))
lines28551 <- readLines(con28551)
close(con28551)
data_start <- which(grepl("^\"ID_REF\"", lines28551))
mat28551 <- as.matrix(read.delim(
  textConnection(lines28551[data_start:length(lines28551)]), row.names = 1))
mat28551 <- mat28551[complete.cases(mat28551), ]
if (max(mat28551, na.rm = TRUE) > 100)
  mat28551 <- log2(pmax(mat28551, 1))

sub28551 <- pdata[grepl("GSE28551", pdata$secondaryaccession), ]
sub28551_ga <- sub28551[sub28551$Gestational.Age.Category %in% valid_cats, ]
idx <- match(sub28551_ga$Sample.Name, colnames(mat28551))
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
idx_sex <- match(sub28551_sex$Sample.Name, colnames(mat28551_ga))
results[[length(results) + 1]] <- run_quantro_safe(
  mat28551_ga[, idx_sex[!is.na(idx_sex)]],
  sub28551_sex$Combined.Fetus.Sex[!is.na(idx_sex)],
  "GSE28551", "ABI Human Genome Survey v2", "Male vs Female (1st tri)"
)

# ============================================================
# GSE93520 (Agilent): 1st only, F:19 M:17
# ============================================================
message("\n=== GSE93520 (Agilent 4x44K) ===")
con93520 <- gzfile(file.path(BASE,
  "data/raws/GSE93520/GSE93520_series_matrix.txt.gz"))
lines93520 <- readLines(con93520)
close(con93520)
data_start <- which(grepl("^\"ID_REF\"", lines93520))
mat93520 <- as.matrix(read.delim(
  textConnection(lines93520[data_start:length(lines93520)]), row.names = 1))
mat93520 <- mat93520[complete.cases(mat93520), ]
if (max(mat93520, na.rm = TRUE) > 100)
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
