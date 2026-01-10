# Expression Data Integration Pipeline

A comprehensive R pipeline for integrating gene expression data from multiple studies with batch correction, PCA, differential expression analysis, protein-protein interaction network clustering, and functional enrichment analysis.

## Features

- **Multi-dataset Integration**: Merge expression matrices from different studies
- **Batch Effect Correction**: Remove technical batch effects using ComBat
- **PCA Analysis**: Visualize data structure and relationships
- **Differential Expression**: Identify DE genes using limma with multiple contrasts
- **Network Clustering**: Build PPI networks from STRING database and detect communities using fastgreedy algorithm
- **Recursive Clustering**: Hierarchical sub-clustering of communities
- **Functional Enrichment**: Pathway enrichment analysis cluster-by-cluster
- **Comprehensive Visualization**: Volcano plots, heatmaps, network plots, and more

## Directory Structure

```
integrative-gene-expression-analysis/
├── R/                                     # Core pipeline modules (source these)
│   ├── 01_data_merging.R                  # Merge expression matrices
│   ├── 02_batch_correction.R              # ComBat batch correction
│   ├── 03_pca_analysis.R                  # PCA and visualization
│   ├── 04_differential_expression.R       # Limma differential expression
│   ├── 05_network_clustering.R            # STRING mapping & fastgreedy clustering
│   ├── 06_enrichment_analysis.R           # Functional enrichment
│   ├── 07_visualization.R                 # Plotting functions
│   ├── config.R                           # Configuration management
│   └── utils.R                            # Utility functions
│
├── scripts/                               # Runnable analysis scripts
│   ├── analysis/                          # Main analysis scripts
│   │   ├── run_separate_analyses.R        # Run trim_1_2 and trim_2_3 analyses
│   │   ├── test_first_datasets.R          # Test dataset combinations (1_2)
│   │   └── test_term_datasets.R           # Test dataset combinations (2_3)
│   ├── validation/                        # Validation and comparison scripts
│   │   └── compare_all_to_disser.R        # Compare results to dissertation
│   ├── preprocessing/                     # Data preprocessing
│   │   ├── preprocessing_illumina/        # Illumina array preprocessing
│   │   └── preprocessing_operon/          # Operon array preprocessing
│   └── visualization/                     # Plotting scripts
│
├── config/                                # YAML configuration files
│   ├── config.yaml                        # Main pipeline config
│   ├── config_separate_analyses.yaml      # Trimester comparison config
│   ├── config_test_first_datasets.yaml    # First trimester testing config
│   └── config_test_term_datasets.yaml     # Term testing config
│
├── notebooks/
│   └── pipeline.Rmd                       # Main analysis notebook (RUN THIS)
│
├── data/
│   ├── mapped/                            # INPUT: Mapped expression files
│   ├── phenodata/                         # INPUT: Sample metadata
│   │   └── samples.csv                    # Main sample metadata
│   ├── raw/                               # Raw downloaded data
│   └── stringdb_cache/                    # STRING database cache
│
├── output/                                # OUTPUT: All results (gitignored)
│   ├── test_first_datasets/               # First trimester dataset testing
│   ├── test_term_datasets/                # Term dataset testing
│   ├── trim_1_2/                          # First vs Second trimester analysis
│   └── trim_2_3/                          # Second vs Third trimester analysis
│
├── one_off_scripts/                       # Temporary/debug scripts (gitignored)
├── tests/                                 # Test scripts
│
├── setup.R                                # Setup script (RUN FIRST)
├── install_packages.R                     # Install required packages
├── CLAUDE.md                              # Claude Code guidance
└── README.md                              # This file
```

## Installation

### Prerequisites

1. **R** (>= 4.0.0)
2. **Bioconductor** packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "limma",
    "sva",
    "org.Hs.eg.db",
    "STRINGdb"
))
```

3. **CRAN** packages:

```r
install.packages(c(
    "yaml",
    "ggplot2",
    "factoextra",
    "igraph",
    "pheatmap",
    "ggrepel",
    "enrichR"
))
```

## Quick Start

### 1. Setup (First Time Only)

```bash
Rscript setup.R
```

This installs all required packages and verifies the installation.

### 2. Prepare Your Data

**Expression Files** (`data/mapped/`):
- Tab-separated values
- Genes as rows (ENTREZID), samples as columns
- First column must be named "ENTREZID"
- Example: `GSE73374_mapped_affymetrix.tsv`

**Phenodata** (`data/phenodata/samples.csv`):
- CSV format with required columns:
  - `arraydatafile_exprscolumnnames` - Sample IDs matching expression file column names
  - `secondaryaccession` - Study/batch ID
  - Biological variables (e.g., `trim_term`, `Diagnosis`, `estimated_sex`)

See `data/phenodata/example_phenodata.csv` for format.

### Validate Your Data

Before running the pipeline, validate metadata coverage:

```bash
Rscript check_metadata_coverage_enhanced.R
```

This checks:
- Metadata completeness (Diagnosis, Gestational Category, Biological Specimen)
- Column name matching between expression files and metadata
- Sample count consistency

### 3. Configure Pipeline

Edit `config/config.yaml` to set:
- Input/output paths
- Batch correction parameters
- Differential expression contrasts
- STRING database settings
- Visualization preferences

Example minimal configuration:

```yaml
paths:
  mapped_data: "data/mapped"
  phenodata: "data/phenodata/samples.csv"
  output: "output"

batch_correction:
  batch_column: "secondaryaccession"
  model_formula: "~trim_term"

differential_expression:
  design_formula: "~0 + trim_term"
  contrasts:
    - "trim_termSecond Trimester - trim_termFirst Trimester"
```

### 4. Run Pipeline

Open and run the RMarkdown notebook:

```r
rmarkdown::render("notebooks/pipeline.Rmd")
```

Or run interactively in RStudio by opening `notebooks/pipeline.Rmd`.

## Usage Examples

### Load Individual Modules

```r
# Source specific modules
source("R/01_data_merging.R")
source("R/02_batch_correction.R")
# ... etc

# Load configuration
config <- load_config("config/config.yaml")

# Merge expression data
merged <- merge_expression_data(
    mapped_path = "data/mapped",
    pattern = "\\.tsv$"
)

# Apply batch correction
corrected <- apply_combat_correction(
    exprs = merged,
    pdata = pdata,
    batch_col = "secondaryaccession",
    mod_formula = ~trim_term
)
```

### Network Clustering

```r
# Initialize STRING database
string_db <- initialize_stringdb(
    version = "11",
    species = 9606,
    score_threshold = 400
)

# Build network
G <- build_ppi_network(difexp$STRING_id, string_db)

# Fastgreedy clustering
fgreedy <- fastgreedy_clustering(G)

# Add cluster assignments
difexp_clustered <- add_cluster_column(difexp, fgreedy, G)
```

### Enrichment Analysis

```r
# Cluster-by-cluster enrichment
enrichment_summary <- cluster_enrichment(
    difexp = difexp_clustered,
    string_db = string_db,
    output_dir = "output",
    max_pvalue = 0.05
)
```

## Output Files

### Main Results

- `exprs_corrected.tsv` - Batch-corrected expression matrix
- `pdata_aligned.tsv` - Aligned phenodata
- `difexp/difexp_filtered.csv` - Filtered differential expression results
- `difexp_final.csv` - Final results with cluster and enrichment info

### Cluster Results

Each cluster has a subdirectory in `output/clusters/[cluster_id]/`:
- `cluster_genes.csv` - Genes in the cluster with sub-cluster assignments
- `cluster_enrichment.csv` - Enrichment results
- `cluster_to_revigo.tsv` - Input file for REVIGO analysis

### Visualizations

- `plots/pca_*.svg` - PCA plots colored by different variables
- `plots/volcano_plot.png` - Volcano plot of differential expression
- `plots/difexp_heatmap.png` - Heatmap of top DE genes
- `plots/network_clusters.png` - Network visualization with clusters
- `plots/cluster_sizes.png` - Bar plot of cluster sizes
- `plots/enrichment_coverage.png` - Enrichment coverage by cluster

### Reports

- `reports/pipeline_summary.csv` - Summary statistics

## Core Algorithms

### ComBat Batch Correction
- Method: Empirical Bayes batch effect correction
- Preserves biological variation while removing technical batches
- Reference: Johnson et al. (2007) Biostatistics

### Differential Expression
- Method: Linear models with empirical Bayes (limma)
- Multiple contrast support
- Reference: Ritchie et al. (2015) NAR

### Network Clustering
- Database: STRING v11 protein-protein interactions
- Algorithm: Fastgreedy community detection
- Recursive sub-clustering of communities
- Reference: Clauset et al. (2004) Phys Rev E

### Functional Enrichment
- STRING enrichment API
- Cluster-by-cluster analysis
- Coverage metrics
- REVIGO-ready output

## Configuration Options

### Key Parameters

| Section | Parameter | Description | Default |
|---------|-----------|-------------|---------|
| **batch_correction** | `batch_column` | Column defining batches | `"secondaryaccession"` |
| | `model_formula` | Biological variables to preserve | `"~trim_term"` |
| **differential_expression** | `design_formula` | Design matrix formula | `"~0 + trim_term"` |
| | `logfc_threshold` | Log fold-change cutoff | `1` |
| | `pvalue_threshold` | Adjusted p-value cutoff | `0.05` |
| **stringdb** | `score_threshold` | Confidence score (0-1000) | `400` |
| **clustering** | `min_cluster_size` | Minimum genes per cluster | `3` |
| | `recursive` | Enable recursive clustering | `TRUE` |
| **enrichment** | `max_pvalue` | P-value threshold | `0.05` |

## Use Case Examples

### Example 1: Trimester Comparison
Compare placental gene expression across trimesters:
```yaml
differential_expression:
  design_formula: "~0 + trim_term"
  contrasts:
    - "trim_termSecond Trimester - trim_termFirst Trimester"
    - "trim_termThird Trimester - trim_termSecond Trimester"
```

### Example 2: Disease vs. Healthy
Compare disease to healthy samples:
```yaml
differential_expression:
  design_formula: "~0 + Diagnosis"
  contrasts:
    - "DiagnosisPre-Eclampsia - DiagnosisHealthy"
```

### Example 3: Sex Differences
Identify sex-specific genes:
```yaml
differential_expression:
  design_formula: "~0 + estimated_sex"
  contrasts:
    - "estimated_sexMale - estimated_sexFemale"
```

## Performance

### Memory Requirements
- Small dataset (<50 samples): 4-8 GB RAM
- Medium dataset (50-200 samples): 8-16 GB RAM
- Large dataset (>200 samples): 16-32 GB RAM

### Runtime Estimates
- Data merging: < 1 min
- Batch correction: 1-5 min
- PCA: < 1 min
- Differential expression: 1-3 min
- STRING mapping: 2-5 min
- Network clustering: 5-15 min
- Enrichment analysis: 5-30 min

**Total: ~15-60 minutes** depending on dataset size and cluster count

## Troubleshooting

### Common Issues

1. **Missing packages**: Install all required packages (see Installation)
2. **Memory errors**: Reduce dataset size or increase memory limit in config
3. **No genes on network**: Lower STRING score_threshold
4. **Empty enrichment**: Increase max_pvalue or check gene IDs
5. **Sample mismatch**: Verify `arraydatafile_exprscolumnnames` matches expression column names exactly

### Check Input Data

```r
source("R/utils.R")

# Validate input
validate_input(exprs, pdata)

# Check configuration
config <- load_config("config/config.yaml")
validate_config(config)
```

### Debug Mode

```r
# Load specific module
source("R/05_network_clustering.R")

# Test function
string_db <- initialize_stringdb(species = 9606, score_threshold = 400)
```

## Citation

If you use this pipeline or database, please cite:

### This Pipeline
- Lykhenko O., Frolova A., Obolenska M. (2021) **Consecutive integration of available microarray data for analysis of differential gene expression in human placenta.** Biotechnologia Acta, 14(1), pp. 38–45. DOI: [10.15407/biotech14.01.38](https://doi.org/10.15407/biotech14.01.38)
- Lykhenko O., Frolova A., Obolenska M. (2017) **Designing the database for microarray experiments metadata.** IEEE International Young Scientists Forum on Applied Physics and Engineering (YSF). DOI: [10.1109/YSF.2017.8126658](https://doi.org/10.1109/YSF.2017.8126658)

### Dependencies
- **limma**: Ritchie ME et al. (2015) Nucleic Acids Research
- **ComBat/sva**: Leek JT et al. (2012) Nature Reviews Genetics
- **STRING**: Szklarczyk D et al. (2021) Nucleic Acids Research
- **igraph**: Csardi G & Nepusz T (2006) InterJournal

## License

This pipeline is provided as-is for research purposes.

## Contact

For questions or issues, please open an issue in the repository or contact the maintainer.

## Recent Updates

### 2026-01-10: Project Reorganization
- **Reorganized project structure**
  - Moved scripts to `scripts/` subfolders (analysis, validation, preprocessing, visualization)
  - Added `one_off_scripts/` for temporary/debug scripts (gitignored)
  - Output folders now match script names
- **Dataset testing improvements**
  - Added `test_first_datasets.R` and `test_term_datasets.R` for testing dataset combinations
  - Scripts now archive previous output before overwriting
  - Added STRING DB excluded genes logging
- **Added CLAUDE.md** for Claude Code guidance

### 2025-10-16: Illumina Preprocessing & Metadata Validation
- **Added Illumina preprocessing pipeline** (`scripts/preprocessing/preprocessing_illumina/`)
  - Download raw data from GEO/ArrayExpress
  - Background correction and quantile normalization
  - Probe-to-gene mapping using max-mean strategy
- **Metadata validation improvements**
  - Enhanced check script (`check_metadata_coverage_enhanced.R`)
  - Column name matching validation

---

**Created**: 2025-10-04
**Last Updated**: 2026-01-10
**Version**: 1.2
