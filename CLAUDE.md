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
- Output folders match script names (e.g., `test_first_datasets.R` → `output/dataset_testing/test_first_datasets/`)
- DEG filtering defaults: `adj.P.Val < 0.05` and `|logFC| > 1`
- Use unambigous date format

## Python Environment

Use `uv` for Python package management:

```bash
# Initialize new Python project
uv init

# Install dependencies
uv sync

# Run Python scripts
uv run python script.py
```

## MCP Vector Database

A local vector database MCP server is available for semantic search over project documentation.

**Location**: `mcp-vectordb/`

**Available tools** (via MCP):
- `search_docs` - Semantic search over all markdown documentation
- `ingest_docs` - Re-index documents after updates
- `list_sources` - List all indexed document sources

**Setup** (if needed):
```bash
cd mcp-vectordb
uv sync
```

**Manual test**:
```bash
cd mcp-vectordb
uv run python -c "from server import ingest_documents; print(f'Ingested {ingest_documents()} chunks')"
```

The MCP server is configured in `.mcp.json` and will be auto-started by Claude Code.
