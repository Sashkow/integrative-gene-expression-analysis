## 1. Article text — dataset references

- [x] 1.1 Replace all "seven datasets"/"7 datasets"/"seven-dataset"/"7-dataset" with "six"/"6" equivalents in `main.tex`
- [x] 1.2 Remove GSE22490 from every dataset enumeration (Methods, Results, Data Availability)
- [x] 1.3 Remove GSE22490 row from Table 1 (dataset overview table)
- [x] 1.4 Remove GSE22490 from tissue description paragraphs in Methods
- [x] 1.5 Add GSE22490 exclusion rationale sentence (decidual contamination) in Methods dataset-selection paragraph

## 2. Article text — numerical results

- [x] 2.1 Update abstract: DEG count 380 → 447, sample count 123 → 117
- [x] 2.2 Update Results section DEG counts: softimpute 380→447 (up 346→369, down 70→78), intersection 182→221 (up 155→174, down 34→47)
- [x] 2.3 Update sample count references throughout (123 → 117)
- [x] 2.4 Update imputation validation Pearson r (0.823 → 0.816) and related accuracy metrics (RMSE, MAE if present)
- [x] 2.5 Update normalization method from "ComBat" to "reference-batch ComBat" (ref=GSE100051) where describing this study's method
- [x] 2.6 Review and update all remaining numerical values in Results and Discussion against 6ds `summary.txt` and `method_comparison.csv`

## 3. Figures

- [x] 3.1 Copy `pca_none_combat_ref.png` and `pca_softimpute_combat_ref.png` from 6ds output to `figures/` (renaming to match existing `\includegraphics` references or updating the references)
- [x] 3.2 Copy NA staircase plots from 6ds output to `figures/`
- [x] 3.3 Re-run phase5 subsampling validation for 6ds config (`config_validation_6ds.yaml`), generate plots, and copy to `figures/`
- [x] 3.4 Verify all `\includegraphics` paths in `main.tex` resolve to existing files
- [x] 3.5 Update figure captions that mention dataset/sample counts

## 4. Data files

- [x] 4.1 Gestational age data file — reusable from 7ds (per-dataset metadata is independent; removing GSE22490 does not affect other datasets' GA data)
- [x] 4.2 Update any `main.tex` references to the gestational age data filename — N/A, no tex references to this file

## 5. Supplements

- [x] 5.1 Review and update `supplement1/supplement1.tex` — no 7ds-specific references found (methodological content only)
- [x] 5.2 Review and update `supplement2/supplement2.tex` — no 7ds-specific references found (RUV control gene content only)

## 6. Subsampling validation text update

- [x] 6.1 Update subsampling validation section text with 6ds numbers: 95.3%→98.3% overlap, 87%→90% retention, Jaccard 0.49→0.54, CCC 0.83→0.88, FDR-only 7k→6.7k, split-half Jaccard 0.38→0.40, etc.
- [x] 6.2 Remove `%% TODO` comments after updating the validation numbers

## 7. Final verification

- [x] 7.1 Grep main.tex for any remaining "GSE22490", "seven", "123 sample", "380", "182" references — all clean (GSE22490 only in Lykhenko 2021 reference and exclusion rationale; "seven" only in pipeline steps count)
- [x] 7.2 Rebuild PDF with `pdflatex main.tex` (twice) and verify no warnings for missing figures (final rebuild after validation updates)
