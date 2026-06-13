## Context

The imputation article (target: BMC Bioinformatics) identifies 447 DEGs (6ds, softImpute + ComBat-ref) but has no biological validation. The autoreview flagged this as the #1 reviewer concern. During exploration we confirmed that the data for the strongest validation already exists: the 2021 Lykhenko full limma table (20,162 genes) overlaps most of the gained DEGs with strong direction concordance.

Current state:
- `output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko/` — all pipeline outputs (DEG tables, expression matrices)
- `output/dissertation/integrative_pipeline_disser/1_2/difexp/difexp_all.csv` — 2021 full limma (20,162 genes)
- `articles/imputation_article/main.tex` — updated for 6ds results, no biological validation sections yet
- Tasks 1 (prior-study concordance) and 2 (GO/KEGG enrichment) were completed for the earlier 7ds run — scripts exist but need re-running with 6ds data paths

## Goals / Non-Goals

**Goals:**
- Re-run completed validation scripts (concordance, enrichment) with 6ds data
- Add mean-imputation benchmark to cross-validation table
- Add ComBat sensitivity analysis
- Structure article text for new Results subsections and Discussion paragraphs

**Non-Goals:**
- Wet-lab experiments (out of scope)
- Literature mining (separate manual effort)
- Prostate cross-tissue validation (pending decision on sysbio article overlap)

## Decisions

### 1. Re-run tasks 1-2 with 6ds paths rather than rewriting scripts

The existing R scripts (`prior_study_concordance.R`, `deg_enrichment.R`, `volcano_arrows*.R`) just need their `base_dir` path updated from `phase2b_1_2_yehor_7ds_no_37653_enriched_sashko` to `phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko`. The analysis logic is unchanged.

### 2. Implementation language: R

The entire pipeline is R-based. All analysis scripts use existing R packages (clusterProfiler, limma, ggplot2, VennDiagram).

### 3. Enrichment uses clusterProfiler with org.Hs.eg.db

Standard Bioconductor workflow. ENTREZID-keyed input. enrichGO (BP) + enrichKEGG. compareCluster for side-by-side comparison of three DEG sets (full 447, intersection 221, gained ~226).

### 4. Mean-imputation benchmark modifies existing cross-validation code

The imputation cross-validation already exists in `scripts/integrative_analysis/phase2b_direct_merge/imputation.R`. Add per-gene-mean and per-batch-mean as additional methods.

### 5. Script location

Validation scripts are in `scripts/integrative_analysis/article_validation/`. Output goes to `output/article_validation/`.

## Risks / Trade-offs

- [Risk] Concordance numbers may differ slightly with 6ds DEGs vs 7ds → Mitigation: the concordance is against the independent 2021 analysis, so the direction of the result (strong concordance) should hold regardless of minor DEG set changes.

- [Risk] Mean-imputation benchmark may require re-running the full cross-validation → Mitigation: mean imputation is trivial to compute so wall time increase is minimal.
