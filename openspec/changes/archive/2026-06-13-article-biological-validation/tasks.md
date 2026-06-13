## 1. Prior-study concordance analysis (re-run for 6ds)

- [x] 1.1 Write R script `scripts/integrative_analysis/article_validation/prior_study_concordance.R` — EXISTS, needs base_dir update from 7ds to 6ds
- [x] 1.2 Archive 7ds outputs: move `output/article_validation/prior_study_concordance/` contents to `output/article_validation/prior_study_concordance/archive/7ds/` (preserve clustered arrows and transition matrix outputs)
- [x] 1.3 Update `base_dir` in `prior_study_concordance.R` and all `volcano_arrows*.R` scripts from `phase2b_1_2_yehor_7ds_no_37653_enriched_sashko` to `phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko`
- [x] 1.4 Re-run `prior_study_concordance.R` to regenerate concordance table, scatter plot, volcano plot, and Venn diagram with 6ds DEGs (447 softImpute, 221 intersection, 231 gained)
- [x] 1.5 Re-run `volcano_arrows_clustered.R` and `transition_matrix.R` to regenerate clustered arrows and transition matrix/heatmap with 6ds data
- [x] 1.6 Verify updated concordance numbers: 231 gained, 211/231 matched (91.3%), 188/211 same direction (89.1%), 141 FDR<0.05 (66.8%), Pearson r=0.606

## 1b. Internal concordance: 6ds softImpute vs 6ds intersection-only

- [x] 1b.1 Write R script `internal_concordance.R` comparing softImpute vs intersection-only on 8,260 shared genes — scatter, volcano, clustered arrows, transition matrix
- [x] 1b.2 logFC scatter plot: r=0.9982, 97.8% same direction across 8,260 shared genes
- [x] 1b.3 Volcano plot of intersection-only full limma with softImpute-only DEGs highlighted (23 genes)
- [x] 1b.4 Clustered arrows plot: 8 clusters from 216 shared DEGs
- [x] 1b.5 Transition matrix/heatmap: 216 DEGs in both, 23 softImpute-only, 5 intersection-only

## 2. GO/KEGG enrichment analysis (re-run for 6ds)

- [x] 2.1 Write R script `scripts/integrative_analysis/article_validation/deg_enrichment.R` — EXISTS, needs base_dir update
- [x] 2.2 Update labels in `deg_enrichment.R` from 380→447, 182→221, 201→226 (base_dir already pointed to 6ds)
- [x] 2.3 Re-run `deg_enrichment.R` — regenerated enrichment tables and dot plots with 6ds DEGs (447/221/231)
- [x] 2.4 Verified: Full 447: 635 GO BP terms, 31 KEGG; Intersection 221: 296 GO BP, 14 KEGG; Gained 226: 247 GO BP, 15 KEGG. Top terms: humoral immune response, cytokine production, leukocyte migration, complement/coagulation cascades

## 3. Mean-imputation benchmark

- [x] 3.1 Added `impute_gene_mean()` and `impute_batch_mean()` to `imputation.R` + registered in IMPUTERS
- [x] 3.2 Ran CV (5 repeats x 3 mask types x 3 methods) via `mean_imputation_benchmark.R`. softImpute dominates: block mask r=0.816 vs baselines r=0.356
- [x] 3.3 Results in `output/article_validation/mean_imputation_benchmark/cv_summary.csv`; baseline rows added to Table 3 in article

## 4. ComBat sensitivity analysis

- [x] 4.1 Write R script `scripts/integrative_analysis/article_validation/combat_sensitivity.R` that extracts ComBat batch parameters (gamma, delta) from intersection-only (8,260 genes) and full imputed (17,531 genes) runs using 6ds output
- [x] 4.2 Compute and plot parameter correlation + report max deviation

## 5. Article text integration

- [x] 5.1 Add Results subsection "Concordance with prior single-platform analysis" — external (Lykhenko 2021) and internal (softImpute vs intersection) concordance
- [x] 5.2 Add Results subsection "Functional enrichment of differentially expressed genes" — GO BP and KEGG for full/intersection/gained sets
- [x] 5.3 Add Discussion subsection "Biological plausibility of gained genes" — three lines of evidence (external concordance, internal concordance, functional coherence)
- [x] 5.4 Added baseline rows (gene mean, batch mean) to Table 3 + paragraph explaining results
- [x] 5.5 Add ComBat sensitivity results (blocked by task 4)
- [x] 5.6 Add Venn diagram figure (fig_venn_three_way.png) and caption to Results — regenerated with 6ds labels (447/221/328)
