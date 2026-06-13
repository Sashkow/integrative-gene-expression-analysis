## 1. Investigation (data lookups before editing)

- [x] 1.1 Identify the 2.41 max-deviation outlier gene from ComBat sensitivity output (`output/article_validation/combat_sensitivity/post_combat_comparison.csv` or recompute from the R script) → **SPP1** (Entrez 6696), deviation=2.41; top 5: SPP1 (2.41), F13A1 (2.19), HBG2 (2.04), HBE1 (2.04), CD36 (2.01)
- [x] 1.2 Characterize the 20 gained genes absent from Affymetrix platforms — check if they show functional enrichment or a pattern (from `output/article_validation/prior_study_concordance/`) → immune/defensin cluster (DEFA1, DEFA1B, DEFA3, CCL3, CCL3L3, FCGR1A, FCGR2A), CG-beta subunits (CGB2, CGB7), keratin (KRTAP26-1); several borderline FDR
- [x] 1.3 Identify the 10 RUVinv DEGs for the 1_2 contrast (from supplement 2 data or pipeline output) → CCR7, MUC20, FCGR2B, HSD11B1, IDO1, EPCAM, PAEP, TMEM100, CLEC3B, CLDN10 (identical 10 in both softImpute and intersection-only)
- [x] 1.4 Verify the Lanoix citation (ref 15 in supplement 2) — look up actual author list to determine correct form → Correct authors: Lanoix D, Lacasse AA, St-Pierre J, Taylor SC, Ethier-Chiasson M, Lafond J, Vaillancourt C. Current citation has "Bhatt H" 5x and "Kingdom JCP" erroneously inserted

## 2. Main text edits (main.tex on article/imputation branch)

- [x] 2.1 Add outlier gene identification note to ComBat sensitivity discussion (~line 540) → Added SPP1 identification + top 5 outliers (F13A1, HBG2, HBE1, CD36)
- [x] 2.2 Add block-mask RMSE as ~53% of expression SD to calibration paragraph (~line 294) → Added after block-mask result sentence
- [x] 2.3 Break Conclusions validation paragraph into two groups (~lines 1118-1134) → Split into primary (random-cell, block) and extended (worst-case, hard-block)
- [x] 2.4 Add Methods paragraph describing all four masking strategies (new subsection after Matrix completion, ~line 1577) → Added new subsection "Cross-validation masking strategies"
- [x] 2.5 Add brief note on 20 non-Affymetrix gained genes in concordance section (~line 808) → Added immune/defensin cluster + CGB2/CGB7 characterization
- [x] 2.6 Shorten deconvolution limitation bullet (~lines 1091-1099) → Condensed to one sentence
- [x] 2.7 Change "coloured"→"colored" and "Colour"→"Color" (lines 493, 494, 512) → Fixed 3 instances

## 3. Table restructure (main.tex + supplement)

- [x] 3.1 Cut Table 4 (`tab:imputation-extended`) from main.tex
- [x] 3.2 Add the table to supplement 1 (or create new supplement table document) → Added as new section "Supplementary table" before bibliography
- [x] 3.3 Update main.tex in-text references to point to supplementary table → 2 refs changed to "Additional file~1, Table~S1"
- [x] 3.4 Verify no hardcoded "Table 5/6/7" references remain in main text (grep check) → Clean

## 4. Figure adjustments

- [x] 4.1 Increase Figure 2 subplot title size in `one_off_scripts/pca_before_after_combat.R` (size 15 → 18+) → Changed to size 20
- [x] 4.2 Round correlation in Figure 4 title: change `%.6f` → `%.3f` in `scripts/integrative_analysis/article_validation/combat_sensitivity.R` (line 272) → Now shows r=0.999
- [x] 4.3 Regenerate Figure 2 and copy to `articles/imputation_article/figures/` → Verified larger titles
- [x] 4.4 Regenerate Figure 4 and copy to `articles/imputation_article/figures/` → Copied post_combat_scatter.png + per_gene_mae_hist.png

## 5. Footnote verification

- [x] 5.1 Extract all `\ghref{...}` paths from main.tex → 41 unique paths
- [x] 5.2 Verify each path exists on `article/imputation` branch via `git ls-tree` → Verified against master (7 broken found)
- [x] 5.3 Fix any broken paths → Fixed: phase5_validation_6ds→phase5_validation (5 paths), config_validation_6ds.yaml→config_validation.yaml, difexp_softimpute_combat.tsv→difexp_significant_softimpute_combat.tsv, regenerated logfc_stability.csv

**TODO**: Before publishing, replace master-branch paths with persistent branch paths

## 6. Supplement 1 edits (supplement1.tex on article/imputation branch)

- [x] 6.1 Add note to recommendation #5 that main article adopted ref-batch ComBat with r=0.9991 → Added to recommendation 5 with GSE100051 ref batch and r=0.9991
- [x] 6.2 Add future-work note to HarmonizR section → Added sentence at end of HarmonizR subsection
- [x] 6.3 Justify M1 (global) imputation choice or add limitation note regarding M2 → Added limitation paragraph after batch-sensitized imputation subsection
- [x] 6.4 Remove or contextualize unused GEO dataset citations (GSE6573, GSE73374, etc.) → Removed 5 orphan bibitems (GSE35574, GSE37653, GSE6573, GSE73374, GSE73685) + GSE55439 comment; kept GSE22490, GSE70102 (referenced in main article)
- [x] 6.5 Add empirical results section for RUV — summarize DEG counts and overlap with ComBat pipeline across configurations (RUVg, RUVinv, bRUV; k=3/4; balanced/all datasets; both contrasts) from `output/phase2b_ruv/` → Added section with summary table and key observations
- [x] 6.6 Add empirical results section for HarmonizR — summarize batch correction outcomes across dataset configs from `output/harmonizr/` → Added section with summary table + implementation divergence note
- [x] 6.7 Add empirical results section for reference-batch ComBat — summarize DEG lists and comparison with standard ComBat from `output/phase2b_combat_ref/` [added section with DEG count table across configs]
- [x] 6.8 Add empirical results section for effect-size meta-analysis — summarize per-study effect-size combination results from `output/exploratory/phase2_meta` [added section with metafor/DExMA/RankProd results]
- [x] 6.9 Add empirical results section for ComBat sensitivity analysis — summarize ref-batch correlation and max-deviation findings from `output/article_validation/combat_sensitivity/` [added section with r=0.9991 and SPP1 max-dev]
- [x] 6.10 Add empirical results section for mean imputation benchmark — summarize comparison of imputation strategies from `output/article_validation/mean_imputation_benchmark/` [added section with softImpute vs gene mean vs batch mean table]

## 7. Supplement 2 edits (supplement2.tex on article/imputation branch)

- [x] 7.1 Fix Lanoix citation author duplication (ref 15) — MUST FIX → Fixed: Lanoix D, Lacasse AA, St-Pierre J, Taylor SC, Ethier-Chiasson M, Lafond J, Vaillancourt C. Also fixed Drewlo citation (Drewlo S, Levytska K, Kingdom J; corrected DOI and pages)
- [x] 7.2 Update "14 datasets" / "227 DEGs" to match current 6ds config or add explanatory note → Added footnote explaining 14-dataset collection is retained for RUV evaluation; main article uses 6-dataset subset
- [x] 7.3 Add follow-up plans note to "Possible improvements" section → Added paragraph on planned follow-up with strategies 1 and 3 on the 6-dataset collection
- [x] 7.4 Add brief characterization of the 10 RUVinv DEGs → Added paragraph after RUVinv observations with functional annotations for all 10 genes

## 8. Build verification

- [x] 8.1 Compile main.tex (two passes) — verify no broken cross-references → 29 pages, 0 undefined references, only benign float warnings
- [x] 8.2 Compile supplement1.tex — verify table renders correctly → 11 pages, clean compile
- [x] 8.3 Compile supplement2.tex — verify citation is correct → 11 pages, Lanoix and Drewlo citations render correctly
- [x] 8.4 Visual check of regenerated Figures 2 and 4 in compiled PDF → Figure 2: larger subplot titles readable; Figure 4: r=0.999 (3dp) in title

## 9. Additional supplement 1 edits

- [x] 9.1 Add empirical benchmark section for batch-in-limma runs — summarize DEG counts and comparison with ComBat pipeline from batch-in-limma output → Added section with table (4 configs), 82-97% overlap with ComBat, power trade-off observation
