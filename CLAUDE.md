# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

R pipeline for integrating gene expression data from multiple GEO studies. Performs batch correction (ComBat), differential expression (limma), PPI network clustering (STRING + fastgreedy), and functional enrichment analysis. Primary use case: placental gene expression across pregnancy trimesters.

## Common Commands

```bash
# Run main pipeline notebook
Rscript -e "rmarkdown::render('notebooks/pipeline.Rmd')"

# Run separate trimester analyses (1_2 and 2_3)
Rscript scripts/analysis/run_separate_analyses.R

# Test dataset combinations for First Trimester comparison
Rscript scripts/analysis/test_first_datasets.R

# Test dataset combinations for Term comparison
Rscript scripts/analysis/test_term_datasets.R

# Compare results to dissertation reference data
Rscript scripts/validation/compare_all_to_disser.R

# Install required packages
Rscript install_packages.R
```

## Data Flow

1. Input: `data/mapped/*.tsv` (ENTREZID-keyed expression), `data/phenodata/samples.csv`
2. Merge → ComBat → PCA → limma DE → STRING mapping → fastgreedy clustering → enrichment
3. Output: `output/{analysis_name}/` with `difexp_filtered.tsv`, `exprs_corrected.tsv`, plots

## Key Patterns

- **Temporary scripts go in `one_off_scripts/`** - Put debug scripts, one-time analyses, and experimental code there
- Scripts archive previous output to `output/{name}/archive/{timestamp}/` before overwriting
- Phenodata sample column `arraydatafile_exprscolumnnames` must match expression file column names
- Output folders match script names (e.g., `test_first_datasets.R` → `output/test_first_datasets/`)
- DEG filtering defaults: `adj.P.Val < 0.05` and `|logFC| > 1`
