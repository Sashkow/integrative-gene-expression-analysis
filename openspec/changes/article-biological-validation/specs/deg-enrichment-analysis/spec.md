## ADDED Requirements

### Requirement: GO/KEGG enrichment of DEG lists

The pipeline MUST run functional enrichment analysis on three DEG sets: the full 380 DEGs (softImpute + ComBat-ref), the 182 intersection-only DEGs, and the ~201 gained DEGs. It MUST use clusterProfiler (enrichGO with org.Hs.eg.db, enrichKEGG) and produce comparison dot plots and enrichment tables for the article.

#### Scenario: Enrichment of full, intersection, and gained DEG sets

Given the DEG lists from `output/yehor_sashko/phase2b_1_2_yehor_7ds_no_37653_enriched_sashko/`
When enrichGO and enrichKEGG are run on all three sets (full 380, intersection 182, gained ~201)
Then enrichment tables and dot plots are produced
And gained DEGs are checked for placental development pathway enrichment
And results are formatted as a new Results subsection and figure in main.tex

#### Scenario: Mean-imputation baseline benchmark

Given the existing cross-validation framework in imputation.R
When per-gene-mean and per-batch-mean imputation are added as baseline methods
Then Table 5 in main.tex is updated with Pearson r, RMSE, MAE for both baselines under random-cell and block masking

#### Scenario: ComBat sensitivity analysis

Given the intersection-only (8,260 genes) and full imputed (17,531 genes) expression matrices
When ComBat is run on both and batch parameters (gamma, delta) are extracted
Then parameter correlation and max deviation are reported as a supplementary table or figure
