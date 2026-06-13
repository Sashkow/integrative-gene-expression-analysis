## Why

The imputation article (BMC Bioinformatics target) identifies 447 DEGs between first- and second-trimester placenta using 6 datasets (softImpute + ComBat-ref), of which ~226 are "gained" by softImpute imputation (447 minus 221 intersection-only DEGs). The central reviewer concern will be: "Are these gained genes biologically real or imputation artefacts?" The article currently has no biological validation — no enrichment, no literature concordance, no cross-tissue test. This change adds layers of biological validation and restructures the article to present them.

A preliminary concordance check during exploration (on the earlier 7ds run) found that ~90% of gained DEGs show concordant logFC direction with the independent 2021 Lykhenko analysis (4 Affy datasets, 22 samples), and ~67% were already significant at FDR < 0.05 in that analysis — lost only because the intersection step discarded them when non-Affy platforms were added. These numbers need re-validation against the current 6ds DEG set, but the pattern is expected to hold. This is strong evidence against the artefact hypothesis.

## What Changes

### Layer 1: Article layout restructuring (main.tex)

Add new Results subsections after the existing Differential Expression section:

- **GO/KEGG enrichment of DEGs** — Functional enrichment analysis of the full 447, intersection-only 221, and gained ~226 DEGs. Table of top enriched terms per group + dot plots.

- **Prior-study concordance of gained DEGs** — Compare the gained DEGs against the full limma table from Lykhenko 2021a (20,162 genes, `output/dissertation/integrative_pipeline_disser/1_2/difexp/difexp_all.csv`). Add logFC scatter plot, volcano plot with gained DEGs highlighted, and categorization table.

- **(Optional) Literature concordance with published placental studies** — For each DEG with published placental qPCR/RT-PCR/IHC/Western data, report published direction vs. pipeline logFC. Separate manual effort.

Add corresponding Discussion paragraphs interpreting enrichment and concordance results.

### Layer 2: Coding tasks

**A. Prior-study concordance analysis script** (tasks 1.x — completed for 7ds, needs re-running for 6ds)
- Input: `output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko/` DEG tables
- Re-run `scripts/integrative_analysis/article_validation/prior_study_concordance.R` with updated 6ds paths
- Update output in `output/article_validation/prior_study_concordance/`

**B. GO/KEGG enrichment script** (tasks 2.x — completed for 7ds, needs re-running for 6ds)
- Input: `difexp_significant_softimpute_combat_ref.tsv` (447 DEGs) + `difexp_significant_none_combat_ref.tsv` (221 intersection DEGs)
- Re-run `scripts/integrative_analysis/article_validation/deg_enrichment.R` with updated 6ds paths
- Update output in `output/article_validation/deg_enrichment/`

**C. Mean-imputation benchmark**
- Add per-gene-mean and per-batch-mean imputation to the cross-validation framework
- Run alongside softImpute under both random-cell and block masking
- Update imputation accuracy table

**D. ComBat sensitivity analysis**
- Compare ComBat gamma/delta batch parameters between intersection-only (8,260 genes) and full imputed (17,531 genes) matrices
- Report correlation and max deviation

**E. Article text integration**
- Add Results subsections and Discussion paragraphs to main.tex
- Add Venn diagram figure and caption

## Capabilities

### New Capabilities
- `deg-enrichment-analysis`: GO/KEGG enrichment pipeline for DEG lists with comparison across pipeline variants
- `prior-study-concordance`: Script to compare gained DEGs against historical full limma output and produce concordance figures
- `literature-concordance`: Curated database of published placental gene expression results with automated comparison to pipeline output

### Modified Capabilities

## Impact

### Files to modify
- `articles/imputation_article/main.tex` — add Results subsections, Discussion paragraphs, new figure/table references
- `scripts/integrative_analysis/article_validation/prior_study_concordance.R` — update base_dir to 6ds output path
- `scripts/integrative_analysis/article_validation/deg_enrichment.R` — update base_dir to 6ds output path
- `scripts/integrative_analysis/article_validation/volcano_arrows*.R` — update base_dir to 6ds output path

### Key data (current 6ds run)
- 447 DEGs in `difexp_significant_softimpute_combat_ref.tsv` (6ds, softImpute + ComBat-ref)
- 221 DEGs in `difexp_significant_none_combat_ref.tsv` (6ds, intersection + ComBat-ref)
- ~226 gained DEGs (447 − 221)
- 20,162 genes in 2021 full limma: `output/dissertation/integrative_pipeline_disser/1_2/difexp/difexp_all.csv`
- Source data: `output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko/`
