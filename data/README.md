# Data Directory Structure

This directory contains all input and intermediate data files for the expression integration pipeline.

## Directory Layout

```
data/
├── raw/                    # Raw expression data files (CEL, IDAT, etc.)
├── mapped/                 # Mapped expression matrices (PIPELINE INPUT)
├── phenodata/              # Sample metadata files
├── merged/                 # Merged expression data
└── temp/                   # Temporary files
```

## Required Input Files

### 1. Mapped Expression Files (`mapped/`)

**Format**: Tab-separated values (TSV)

**Structure**:
- **Rows**: Genes (ENTREZID as row names)
- **Columns**: Samples (sample IDs as column names)
- **Values**: Normalized expression values

**Example**:
```
ENTREZID    Sample1    Sample2    Sample3
1           7.234      7.891      6.543
10          5.678      5.234      5.890
100         8.901      9.123      8.765
...
```

**File naming**: `[STUDY_ID]_mapped_[PLATFORM].tsv`
- Examples:
  - `E-GEOD-73374_mapped_affymetrix.tsv`
  - `GSE12767_mapped_illumina.tsv`

### 2. Phenodata File (`phenodata/samples.csv`)

**Format**: CSV (comma-separated values)

**Required columns**:
- `arraydatafile_exprscolumnnames` - Sample IDs matching expression file column names
- `secondaryaccession` - Study/batch identifier (or custom batch variable)
- Biological variables for your analysis (e.g., `trim_term`, `Diagnosis`, `estimated_sex`)

**Recommended columns**:
- `Biological.Specimen` - Tissue type
- `Gestational.Age.Appr` - Gestational age (for trimester calculation)
- `Gestational.Age.Category` - Gestational age category
- `Diagnosis` - Clinical diagnosis
- `estimated_sex` - Estimated sex
- `Fetus.Sex` - Actual fetal sex (if known)

**Example**:
```csv
arraydatafile_exprscolumnnames,secondaryaccession,Biological.Specimen,Gestational.Age.Appr,trim_term,Diagnosis,estimated_sex
Sample1,GSE73374,Placenta,8,First Trimester,Healthy,Male
Sample2,GSE73374,Placenta,10,First Trimester,Healthy,Female
Sample3,GSE12767,Placenta,20,Second Trimester,Healthy,Male
Sample4,GSE12767,Placenta,22,Second Trimester,Healthy,Female
...
```

## Output Directories

### Merged Data (`merged/`)

Contains merged expression matrices and aligned phenodata after data integration.

**Files**:
- `mrgd.tsv` - Merged expression matrix (common genes across all studies)
- `pdata.tsv` - Aligned phenodata

### Temporary Files (`temp/`)

Contains intermediate files used during processing.

**Files**:
- `mrgd.tsv` - Temporary merged data
- `pdata.tsv` - Temporary phenodata
- Other processing files

## Data Preparation Guidelines

### 1. Expression Data Requirements

✅ **DO**:
- Use gene-level expression data (not probe-level)
- Map to common gene identifiers (ENTREZID recommended)
- Ensure consistent normalization across studies
- Remove genes with excessive missing values
- Use meaningful sample column names

❌ **DON'T**:
- Mix different normalization methods
- Include samples with missing metadata
- Use duplicate gene IDs
- Include completely missing genes

### 2. Phenodata Requirements

✅ **DO**:
- Use consistent column names across studies
- Include all relevant biological and technical variables
- Document factor levels clearly
- Handle missing values appropriately (use "" or NA)

❌ **DON'T**:
- Use special characters in column names
- Mix data types in same column
- Leave critical variables undefined
- Use inconsistent factor level naming

### 3. Sample ID Matching

**Critical**: Sample IDs in phenodata column `arraydatafile_exprscolumnnames` must match expression file column names.

R's `make.names()` function is used to standardize column names:
- Spaces → periods (`.`)
- Special characters → periods
- Leading numbers → prepended with `X`

**Example matching**:
```
Expression file column: "Sample 1" or "Sample.1"
Phenodata value:        "Sample 1"  (will be converted to "Sample.1")
```

## Example Data Structure

### Study 1: GSE73374 (Affymetrix)
```
data/mapped/GSE73374_mapped_affymetrix.tsv
- 18,000 genes × 30 samples
- First and Second Trimester placenta
```

### Study 2: GSE12767 (Illumina)
```
data/mapped/GSE12767_mapped_illumina.tsv
- 20,000 genes × 25 samples
- Term placenta
```

### Combined Phenodata
```
data/phenodata/samples.csv
- 55 samples total
- Columns: sample IDs, study, tissue, gestational age, diagnosis, sex
```

### After Merging
```
data/merged/mrgd.tsv
- 16,500 common genes × 55 samples
```

## Quality Control Checks

Before running the pipeline:

1. **Check file format**:
   ```r
   exprs <- read.table("data/mapped/study1.tsv", header=TRUE, sep="\t", row.names=1)
   head(exprs)
   ```

2. **Check phenodata**:
   ```r
   pdata <- read.csv("data/phenodata/samples.csv")
   head(pdata)
   table(pdata$secondaryaccession)
   ```

3. **Verify sample matching**:
   ```r
   make.names(pdata$arraydatafile_exprscolumnnames) %in% colnames(exprs)
   ```

4. **Check for missing data**:
   ```r
   sum(is.na(exprs))
   sum(is.na(pdata))
   ```

## File Size Estimates

Typical file sizes:
- Mapped expression file: 5-50 MB per study
- Phenodata: < 1 MB
- Merged data: 10-200 MB depending on number of studies

Large dataset tips:
- Use data.table for faster reading
- Consider subsetting genes before merging
- Increase memory limits in R configuration

---

**Note**: Always backup original data files before processing!
