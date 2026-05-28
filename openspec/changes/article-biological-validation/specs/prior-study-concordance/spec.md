## ADDED Requirements

### Requirement: Prior-study concordance of gained DEGs

The analysis MUST compare the ~201 gained DEGs (present in softImpute result but absent from intersection-only result) against the full limma table from Lykhenko 2021a (`output/dissertation/integrative_pipeline_disser/1_2/difexp/difexp_all.csv`, 20,162 genes, 4 Affy datasets, 22 samples). Exploration confirmed: 181/201 testable in 2021, 163/181 (90.1%) same logFC direction, 122 already FDR < 0.05 in 2021. This is the strongest rebuttal to the "imputation artefact" concern.

#### Scenario: Concordance categorization

Given the 201 gained DEGs and the 2021 full limma output
When each gained DEG is matched to the 2021 table
Then each gene is categorized as: same-direction-significant (FDR<0.05), same-direction-near-sig (FDR 0.05–0.20), same-direction-weak (FDR>0.20), opposite-direction, or absent (non-Affy platform genes)
And a summary table with counts per category is produced

#### Scenario: logFC scatter plot

Given matched gained DEGs with logFC values from both analyses
When logFC from the current 7-dataset study is plotted against logFC from the 2021 4-Affy study
Then a scatter plot is produced showing correlation with Pearson r annotation
And points are colour-coded by 2021 FDR significance

#### Scenario: Volcano plot with gained DEGs highlighted

Given the 2021 full limma table (all 20,162 genes)
When a volcano plot is drawn (-log10 FDR vs logFC)
Then the 181 gained DEGs found in 2021 are highlighted in colour
And the viewer can see that most gained DEGs already sat in the significant region of the 2021 analysis

#### Scenario: Three-way Venn diagram

Given the DEG lists from: (a) this study softImpute 380, (b) this study intersection-only 182, (c) Lykhenko 2021a 327
When a Venn diagram is drawn
Then it shows the overlap structure across all three analyses
And demonstrates that gained DEGs largely overlap with 2021 significant genes
