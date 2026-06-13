## ADDED Requirements

### Requirement: Outlier gene identification in ComBat sensitivity
The ComBat sensitivity discussion SHALL identify the gene with the maximum deviation of 2.41 log2 and briefly characterize it (e.g., coverage pattern, known volatility, or sparsity boundary status).

#### Scenario: Max-deviation gene described
- **WHEN** the reader reads the ComBat sensitivity results (near the r=0.9991 discussion)
- **THEN** they find a brief note identifying which gene produced the 2.41 max deviation and what characteristic explains it

### Requirement: Block-mask RMSE calibration as percentage of SD
Section 2.2 (Imputation accuracy) SHALL include the block-mask RMSE (1.52) expressed as a percentage of the expression standard deviation (~53%), alongside the existing random-cell RMSE calibration (~12%).

#### Scenario: Block-mask RMSE calibrated
- **WHEN** the reader reaches the RMSE calibration paragraph
- **THEN** they see both the random-cell RMSE (~12% of variability) and the block-mask RMSE (~53% of variability) presented for comparison

### Requirement: Conclusions paragraph restructured
The Conclusions paragraph covering validation results (approximately lines 1118-1134) SHALL be broken into two sentences or groups: one for primary results (random-cell, block) and one for extended validation (worst-case, hard-block thresholds).

#### Scenario: Conclusions paragraph is readable
- **WHEN** the reader reads the Conclusions section
- **THEN** the validation summary is structured as two clear groups rather than a single run-on passage

### Requirement: Masking strategies described in Methods
The Methods section SHALL include a paragraph describing all four cross-validation masking strategies: (1) random-cell, (2) gene-dataset block, (3) worst-case block, and (4) hard block mask. Currently only random-cell and gene-dataset block are described (implicitly, in Results).

#### Scenario: Methods section describes all masking types
- **WHEN** the reader reads the Methods section
- **THEN** they find explicit descriptions of all four masking types used in the imputation validation, located before the Results section references them

### Requirement: Brief note on 20 non-Affymetrix gained genes
The concordance section (Section 2.6, near line 808) SHALL include a brief note on the 20 gained genes absent from Affymetrix platforms — whether they show functional enrichment or any particular pattern.

#### Scenario: Non-Affymetrix genes addressed
- **WHEN** the reader reads the concordance section's discussion of gained genes
- **THEN** they find a brief characterization of the 20 genes absent from Affymetrix platforms

### Requirement: Shortened deconvolution limitation bullet
The deconvolution limitation bullet (Limitations, ~lines 1091-1099) SHALL be shortened to approximately one sentence, such as: "Exploratory dtangle deconvolution confirmed trophoblast dominance (>=90%) across all datasets, but could not resolve decidual contamination at bulk resolution."

#### Scenario: Deconvolution bullet is concise
- **WHEN** the reader reads the Limitations section
- **THEN** the deconvolution bullet is no longer than ~2 lines

### Requirement: Consistent American English spelling
The article SHALL use consistent American English spelling throughout. Specifically, "coloured" (lines 493-494) and "Colour" (line 512) SHALL be changed to "colored" and "Color".

#### Scenario: No British spelling variants remain
- **WHEN** the text is searched for "colour" (case-insensitive)
- **THEN** zero matches are found; all instances use "color"
