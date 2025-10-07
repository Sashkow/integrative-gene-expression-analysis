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
exprs_integration_pipeline/
├── R/                          # Core pipeline functions
│   ├── 01_data_merging.R       # Data merging functions
│   ├── 02_batch_correction.R   # Batch correction with ComBat
│   ├── 03_pca_analysis.R       # PCA and visualization
│   ├── 04_differential_expression.R  # Limma DE analysis
│   ├── 05_network_clustering.R # STRING mapping and fastgreedy clustering
│   ├── 06_enrichment_analysis.R # Functional enrichment
│   ├── 07_visualization.R      # Plotting functions
│   ├── config.R                # Configuration management
│   └── utils.R                 # Utility functions
├── config/
│   └── config.yaml             # Pipeline configuration file
├── notebooks/
│   └── pipeline.Rmd            # Main analysis notebook
├── data/
│   ├── raw/                    # Raw expression data
│   ├── mapped/                 # Mapped expression matrices (input)
│   ├── phenodata/              # Sample metadata
│   ├── merged/                 # Merged data
│   └── temp/                   # Temporary files
├── output/
│   ├── difexp/                 # Differential expression results
│   ├── clusters/               # Network clustering results
│   ├── enrichment/             # Enrichment analysis results
│   ├── plots/                  # Visualizations
│   └── reports/                # Summary reports
├── tests/                      # Unit tests (optional)
├── docs/                       # Additional documentation
└── README.md                   # This file
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

### 1. Prepare Your Data

Place your mapped expression files in `data/mapped/`:
- Files should be tab-separated with genes as rows, samples as columns
- First column should be gene IDs (ENTREZID)

Place your phenodata file in `data/phenodata/samples.csv`:
- Must include column `arraydatafile_exprscolumnnames` with sample IDs matching expression files
- Should include batch variable (e.g., `secondaryaccession`) and biological variables

### 2. Configure Pipeline

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

### 3. Run Pipeline

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

## Troubleshooting

### Common Issues

1. **Missing packages**: Install all required packages (see Installation)
2. **Memory errors**: Reduce dataset size or increase memory limit in config
3. **No genes on network**: Lower STRING score_threshold
4. **Empty enrichment**: Increase max_pvalue or check gene IDs

### Check Input Data

```r
source("R/utils.R")

# Validate input
validate_input(exprs, pdata)

# Check configuration
config <- load_config("config/config.yaml")
validate_config(config)
```

## Citation

If you use this pipeline, please cite:

- **limma**: Ritchie ME et al. (2015) Nucleic Acids Research
- **ComBat/sva**: Leek JT et al. (2012) Nature Reviews Genetics
- **STRING**: Szklarczyk D et al. (2021) Nucleic Acids Research
- **igraph**: Csardi G & Nepusz T (2006) InterJournal

## License

This pipeline is provided as-is for research purposes.

## Contact

For questions or issues, please open an issue in the repository or contact the maintainer.

---

**Last updated**: 2025-10-04
