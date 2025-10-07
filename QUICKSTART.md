# Quick Start Guide

Get started with the Expression Integration Pipeline in 5 minutes!

## Prerequisites

- R >= 4.0.0
- RStudio (recommended) or R command line

## Step 1: Install Dependencies (2 min)

```bash
cd scripts/exprs_integration_pipeline
Rscript setup.R
```

This will:
- ✓ Install all required Bioconductor and CRAN packages
- ✓ Verify all modules load correctly
- ✓ Create necessary directory structure

## Step 2: Prepare Your Data (5 min)

### A. Expression Files

Place your mapped expression files in `data/mapped/`:

```bash
# Example files (tab-separated, genes × samples)
data/mapped/GSE73374_mapped_affymetrix.tsv
data/mapped/GSE12767_mapped_illumina.tsv
data/mapped/GSE9984_mapped_affymetrix.tsv
```

**File format**:
```
ENTREZID    Sample1    Sample2    Sample3
1           7.234      7.891      6.543
10          5.678      5.234      5.890
...
```

### B. Phenodata

Create `data/phenodata/samples.csv`:

```csv
arraydatafile_exprscolumnnames,secondaryaccession,Biological.Specimen,Gestational.Age,trim_term,Diagnosis,estimated_sex
Sample1,GSE73374,Placenta,39,Term,Healthy,Male
Sample2,GSE73374,Placenta,38,Term,Healthy,Female
...
```

**Required columns**:
- `arraydatafile_exprscolumnnames` - Sample IDs (must match expression file columns)
- `secondaryaccession` - Study/batch identifier
- Biological variables for your analysis

**See**: `data/phenodata/example_phenodata.csv` for full example

## Step 3: Configure Pipeline (2 min)

Edit `config/config.yaml`:

```yaml
# Minimal configuration
paths:
  mapped_data: "data/mapped"
  phenodata: "data/phenodata/samples.csv"

batch_correction:
  batch_column: "secondaryaccession"
  model_formula: "~trim_term"

differential_expression:
  design_formula: "~0 + trim_term"
  contrasts:
    - "trim_termTerm - trim_termSecond Trimester"

stringdb:
  score_threshold: 400
```

**Key parameters to set**:
- `batch_column`: Column defining batches
- `model_formula`: Biological variables to preserve
- `design_formula`: Design matrix for DE analysis
- `contrasts`: Comparisons to test

## Step 4: Run Pipeline (15-60 min)

### Option A: RStudio (Recommended)

1. Open `notebooks/pipeline.Rmd` in RStudio
2. Click "Knit" button
3. Wait for analysis to complete
4. View HTML report

### Option B: Command Line

```bash
Rscript -e "rmarkdown::render('notebooks/pipeline.Rmd')"
```

## Step 5: View Results

### Main Output Directory: `output/`

```
output/
├── exprs_corrected.tsv         # Batch-corrected expression
├── pdata_aligned.tsv           # Aligned metadata
├── difexp_final.csv            # Final DE results with clusters
├── difexp/
│   ├── difexp_all.csv          # All DE results
│   └── difexp_filtered.csv     # Filtered DE results
├── clusters/
│   ├── 1/                      # Cluster 1
│   │   ├── cluster_genes.csv
│   │   ├── cluster_enrichment.csv
│   │   └── cluster_to_revigo.tsv
│   ├── 2/                      # Cluster 2
│   └── ...
├── plots/
│   ├── pca_trim_term_PC1_PC2.svg
│   ├── volcano_plot.png
│   ├── network_clusters.png
│   └── ...
└── reports/
    └── pipeline_summary.csv
```

### View HTML Report

Open `notebooks/pipeline.html` in your browser for interactive results with:
- Step-by-step analysis
- All figures embedded
- Summary statistics
- Session information

## What the Pipeline Does

### 1. Data Merging
- Combines expression matrices from multiple studies
- Keeps only common genes
- Aligns with phenodata

### 2. Batch Correction (ComBat)
- Removes technical batch effects
- Preserves biological variation
- Tests for confounding

### 3. PCA Analysis
- Visualizes data structure
- Colors by metadata variables
- Identifies outliers

### 4. Differential Expression (limma)
- Tests specified contrasts
- Multiple testing correction
- Annotates with gene symbols

### 5. Network Clustering
- Maps genes to STRING database
- Builds PPI network
- **Fastgreedy community detection**
- Recursive sub-clustering

### 6. Functional Enrichment
- Cluster-by-cluster enrichment
- GO Biological Process
- Coverage metrics
- REVIGO-ready output

## Common Use Cases

### Trimester Comparison
```yaml
differential_expression:
  design_formula: "~0 + trim_term"
  contrasts:
    - "trim_termSecond Trimester - trim_trimFirst Trimester"
    - "trim_termTerm - trim_termSecond Trimester"
```

### Disease vs. Control
```yaml
differential_expression:
  design_formula: "~0 + Diagnosis"
  contrasts:
    - "DiagnosisPre-Eclampsia - DiagnosisHealthy"
```

### Sex Differences
```yaml
differential_expression:
  design_formula: "~0 + estimated_sex"
  contrasts:
    - "estimated_sexMale - estimated_sexFemale"
```

## Troubleshooting

### Issue: Setup fails
**Solution**: Check R version (>= 4.0.0) and internet connection

### Issue: Sample mismatch error
**Solution**: Ensure `arraydatafile_exprscolumnnames` matches expression column names
```r
# In R, check matching:
pdata <- read.csv("data/phenodata/samples.csv")
exprs <- read.table("data/mapped/study1.tsv", header=TRUE, sep="\t", row.names=1)
make.names(pdata$arraydatafile_exprscolumnnames) %in% colnames(exprs)
```

### Issue: No genes on network
**Solution**: Lower STRING `score_threshold` in config (try 200 or 100)

### Issue: Memory error
**Solution**: Reduce dataset size or increase R memory:
```r
# Before running
memory.limit(16000)  # Windows
# Or use ulimit on Linux/Mac
```

### Issue: Empty enrichment
**Solution**: Increase `max_pvalue` or check gene IDs are ENTREZID

## Getting Help

1. **Documentation**: See `README.md` for detailed docs
2. **Data format**: See `data/README.md`
3. **Example data**: See `data/phenodata/example_phenodata.csv`
4. **Full summary**: See `PROJECT_SUMMARY.md`

## Next Steps

After successful run:

1. **Explore clusters**: Check `output/clusters/*/cluster_enrichment.csv`
2. **REVIGO analysis**: Upload `cluster_to_revigo.tsv` to http://revigo.irb.hr/
3. **Custom analysis**: Modify `notebooks/pipeline.Rmd`
4. **Iterate**: Adjust config and re-run with different parameters

## Example Session

```r
# In R/RStudio
library(knitr)
library(rmarkdown)

# Run pipeline
render("notebooks/pipeline.Rmd")

# Or run step-by-step interactively
source("R/config.R")
config <- load_config("config/config.yaml")

source("R/01_data_merging.R")
merged <- merge_expression_data(config$paths$mapped_data)

# ... continue with other steps
```

---

**Time to first results**: ~25 minutes (including setup)

**Need help?** Check the documentation or examine the example data!
