## ADDED Requirements

### Requirement: PCA plots replaced with 6ds versions
The PCA figures SHALL be replaced with plots from `output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko/`. Both intersection (`pca_none_combat_ref.png`) and softImpute (`pca_softimpute_combat_ref.png`) PCA plots SHALL be updated.

#### Scenario: PCA figure files updated
- **WHEN** checking `articles/imputation_article/figures/`
- **THEN** `fig_pca_intersection.png` and `fig_pca_softimpute.png` SHALL contain 6ds combat_ref PCA plots showing 6 datasets (no GSE22490 cluster)

### Requirement: NA staircase plots replaced with 6ds versions
The NA staircase figures SHALL be replaced with 6ds versions from the output directory.

#### Scenario: Staircase figure files updated
- **WHEN** checking `articles/imputation_article/figures/`
- **THEN** staircase plots SHALL show the 6-dataset coverage pattern

### Requirement: Validation plots replaced with 6ds versions
All validation figures (CCC convergence, retention, split-half) SHALL be replaced with 6ds equivalents if available in the output directory.

#### Scenario: Validation figures reflect 6ds
- **WHEN** viewing validation figures
- **THEN** they SHALL show metrics from the 6-dataset run

### Requirement: LaTeX figure references match filenames
All `\includegraphics` paths in `main.tex` SHALL resolve to existing files in the figures directory after the update.

#### Scenario: No broken figure references
- **WHEN** compiling `main.tex` with pdflatex
- **THEN** no "file not found" warnings SHALL appear for figure includes

### Requirement: Gestational age data file updated
The gestational age data file SHALL be updated from `data/gestational_age_7ds.xlsx` to a 6ds version reflecting only the 6 included datasets.

#### Scenario: GA data file reflects 6 datasets
- **WHEN** reading the gestational age data file
- **THEN** it SHALL contain data for exactly 6 datasets without GSE22490
