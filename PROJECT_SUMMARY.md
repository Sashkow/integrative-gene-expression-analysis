# Expression Integration Pipeline - Project Summary

## Overview

This is a comprehensive, modular R pipeline for integrating gene expression data from multiple studies with:
- Batch effect correction using ComBat
- PCA analysis and visualization
- Differential expression analysis using limma
- Protein-protein interaction network construction from STRING database
- Community detection using fastgreedy algorithm
- Recursive hierarchical clustering
- Functional enrichment analysis cluster-by-cluster

## Project Structure

```
exprs_integration_pipeline/
├── R/                                  # Core pipeline modules (modular functions)
│   ├── 01_data_merging.R               # Merge expression matrices
│   ├── 02_batch_correction.R           # ComBat batch correction
│   ├── 03_pca_analysis.R               # PCA and visualization
│   ├── 04_differential_expression.R    # Limma differential expression
│   ├── 05_network_clustering.R         # STRING mapping & fastgreedy clustering
│   ├── 06_enrichment_analysis.R        # Functional enrichment
│   ├── 07_visualization.R              # Plotting functions
│   ├── config.R                        # Configuration management
│   └── utils.R                         # Utility functions
│
├── config/
│   └── config.yaml                     # Pipeline configuration (EDIT THIS)
│
├── notebooks/
│   └── pipeline.Rmd                    # Main analysis notebook (RUN THIS)
│
├── data/
│   ├── mapped/                         # INPUT: Mapped expression files
│   ├── phenodata/                      # INPUT: Sample metadata
│   │   └── example_phenodata.csv       # Example phenodata format
│   ├── merged/                         # Merged expression data
│   ├── temp/                           # Temporary files
│   └── README.md                       # Data documentation
│
├── output/                             # OUTPUT: All results
│   ├── difexp/                         # Differential expression results
│   ├── clusters/                       # Network clustering results
│   │   └── [cluster_id]/               # Per-cluster subdirectories
│   │       ├── cluster_genes.csv       # Genes with sub-clusters
│   │       ├── cluster_enrichment.csv  # Enrichment results
│   │       └── cluster_to_revigo.tsv   # REVIGO input
│   ├── enrichment/                     # Enrichment summaries
│   ├── plots/                          # All visualizations
│   └── reports/                        # Summary reports
│
├── setup.R                             # Setup script (RUN FIRST)
├── README.md                           # Main documentation
├── .gitignore                          # Git ignore rules
└── PROJECT_SUMMARY.md                  # This file
```

## Key Features

### 1. Modular Design
- Each step is in a separate R file with documented functions
- Easy to use individual modules independently
- Clean separation of concerns

### 2. Configuration-Driven
- All parameters in `config/config.yaml`
- No need to edit code for different analyses
- Easy to reproduce analyses with different settings

### 3. RMarkdown Notebook
- Step-by-step execution with explanations
- Generates HTML report with all results
- Interactive analysis in RStudio

### 4. Comprehensive Pipeline
- **Data Merging**: Combine multiple datasets by common genes
- **Batch Correction**: Remove technical variation with ComBat
- **PCA**: Visualize data structure
- **Differential Expression**: Limma with multiple contrasts
- **Network Clustering**: STRING PPI network + fastgreedy algorithm
- **Recursive Clustering**: Hierarchical sub-clustering of communities
- **Enrichment**: Cluster-by-cluster pathway analysis

### 5. Rich Visualizations
- PCA plots colored by metadata
- Volcano plots
- Heatmaps
- Network plots with cluster colors
- Cluster size and coverage plots

## Quick Start

### 1. Setup (First Time Only)

```bash
cd scripts/exprs_integration_pipeline
Rscript setup.R
```

This installs all required packages and verifies the installation.

### 2. Prepare Data

**Expression Files** (`data/mapped/`):
- Tab-separated values
- Genes as rows (ENTREZID), samples as columns
- Example: `GSE73374_mapped_affymetrix.tsv`

**Phenodata** (`data/phenodata/samples.csv`):
- CSV format with required columns:
  - `arraydatafile_exprscolumnnames` - Sample IDs
  - `secondaryaccession` - Study/batch ID
  - Biological variables (e.g., `trim_term`, `Diagnosis`, `estimated_sex`)

See `data/README.md` and `data/phenodata/example_phenodata.csv` for details.

### 3. Configure Pipeline

Edit `config/config.yaml`:

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
    - "trim_termTerm - trim_termSecond Trimester"

stringdb:
  score_threshold: 400  # Confidence score (0-1000)

clustering:
  recursive: TRUE
  min_cluster_size: 3

enrichment:
  max_pvalue: 0.05
```

### 4. Run Pipeline

**In R/RStudio**:
```r
rmarkdown::render("notebooks/pipeline.Rmd")
```

**From command line**:
```bash
Rscript -e "rmarkdown::render('notebooks/pipeline.Rmd')"
```

This generates `notebooks/pipeline.html` with complete analysis report.

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

## Output Files

### Main Results
- `exprs_corrected.tsv` - Batch-corrected expression matrix
- `pdata_aligned.tsv` - Aligned phenodata
- `difexp_final.csv` - Final DE results with clusters and enrichment

### Cluster Results
Each cluster in `output/clusters/[cluster_id]/`:
- `cluster_genes.csv` - Genes with sub-cluster assignments
- `cluster_enrichment.csv` - Enrichment results
- `cluster_to_revigo.tsv` - Input for REVIGO

### Visualizations
- `plots/pca_*.svg` - PCA plots
- `plots/volcano_plot.png` - Volcano plot
- `plots/network_clusters.png` - Network with clusters
- `plots/cluster_sizes.png` - Cluster sizes
- `plots/enrichment_coverage.png` - Coverage metrics

## Use Cases

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

## Customization

### Add New Modules
Create new R file in `R/` directory:
```r
#' My Custom Analysis Module
#'
#' @param data Input data
#' @return Analysis results
#' @export
my_custom_function <- function(data) {
  # Your code here
}
```

Source in notebook:
```r
source("../R/my_custom_module.R")
```

### Modify Pipeline
Edit `notebooks/pipeline.Rmd` to:
- Add/remove steps
- Change parameters
- Add custom visualizations
- Integrate external tools

## Dependencies

### Bioconductor
- limma - Differential expression
- sva - Batch correction (ComBat)
- org.Hs.eg.db - Gene annotation
- STRINGdb - Protein interactions
- AnnotationDbi - Annotation interface

### CRAN
- yaml - Configuration
- ggplot2 - Plotting
- factoextra - PCA visualization
- igraph - Network analysis
- pheatmap - Heatmaps
- ggrepel - Label repulsion
- enrichR - Enrichment (optional)

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

1. **Memory errors**: Reduce dataset size or increase R memory limit
2. **Package installation fails**: Check Bioconductor version compatibility
3. **No genes on network**: Lower STRING `score_threshold`
4. **Empty enrichment**: Check gene IDs and increase `max_pvalue`
5. **Sample mismatch**: Verify `arraydatafile_exprscolumnnames` matches expression column names

### Debug Mode

```r
# Load specific module
source("R/05_network_clustering.R")

# Test function
string_db <- initialize_stringdb(species = 9606, score_threshold = 400)
```

## Citation

If you use this pipeline, cite:
- **limma**: Ritchie ME et al. (2015) NAR
- **ComBat**: Johnson WE et al. (2007) Biostatistics
- **STRING**: Szklarczyk D et al. (2021) NAR
- **igraph/fastgreedy**: Clauset A et al. (2004) Phys Rev E

## License

Research use only. See individual package licenses for component tools.

## Support

- Documentation: See `README.md` and `data/README.md`
- Example data: `data/phenodata/example_phenodata.csv`
- Setup help: Run `Rscript setup.R`

---

**Created**: 2025-10-04
**Version**: 1.0
**Pipeline Type**: Expression Integration with Network Clustering
