## ADDED Requirements

### Requirement: Dataset count references updated to six
All references to "seven datasets", "7 datasets", "seven-dataset", and "7-dataset" in `main.tex` SHALL be replaced with "six datasets", "6 datasets", "six-dataset", and "6-dataset" respectively.

#### Scenario: Abstract dataset count
- **WHEN** reading the abstract
- **THEN** it SHALL state "six placental/chorionic-villi datasets"

#### Scenario: Body text dataset counts
- **WHEN** searching main.tex for "seven" or "7 dataset"
- **THEN** zero matches SHALL remain (except in references to other studies)

### Requirement: GSE22490 removed from all dataset listings
All mentions of GSE22490 SHALL be removed from main.tex, including dataset enumerations, Table 1, tissue description paragraphs, and data availability statements. The dataset list SHALL be: GSE100051, GSE122214, GSE28551, GSE37901, GSE93520, GSE9984.

#### Scenario: Table 1 datasets
- **WHEN** reading Table 1 (dataset overview)
- **THEN** it SHALL list exactly 6 datasets with no GSE22490 row

#### Scenario: Methods dataset enumeration
- **WHEN** reading the Methods section dataset list
- **THEN** GSE22490 SHALL NOT appear

#### Scenario: Data availability
- **WHEN** reading the data availability statement
- **THEN** the GEO accession list SHALL contain exactly 6 accessions without GSE22490

### Requirement: GSE22490 exclusion rationale documented
The Methods section SHALL include a brief statement explaining that GSE22490 was excluded due to decidual tissue contamination in its collection protocol.

#### Scenario: Exclusion rationale present
- **WHEN** reading the Methods dataset selection paragraph
- **THEN** it SHALL mention GSE22490 exclusion with decidual contamination as the reason

### Requirement: DEG counts updated to 6ds results
All differential expression result numbers SHALL match the 6ds `softimpute_combat_ref` output: 447 total DEGs (369 up, 78 down) for the full imputed set, and 221 DEGs (174 up, 47 down) for intersection-only.

#### Scenario: Abstract DEG count
- **WHEN** reading the abstract
- **THEN** it SHALL report 447 DEGs (not 380)

#### Scenario: Results section DEG counts
- **WHEN** reading the Results section method comparison
- **THEN** counts SHALL match: softimpute_combat_ref 447 (369 up, 78 down), none_combat_ref 221 (174 up, 47 down)

### Requirement: Sample count updated
All references to sample counts SHALL reflect the 6ds run: 117 total samples (not 123).

#### Scenario: Sample count in Results
- **WHEN** reading sample count references
- **THEN** they SHALL state 117 samples

### Requirement: Imputation validation correlation updated
The imputation validation Pearson r SHALL be reported as 0.816 (from 6ds block-mask validation), not 0.823.

#### Scenario: Validation correlation in Results
- **WHEN** reading imputation accuracy results
- **THEN** Pearson r SHALL be ~0.816

### Requirement: Normalization method references updated
References to the normalization method SHALL reflect reference-batch ComBat (combat_ref with GSE100051 as reference), not standard ComBat.

#### Scenario: Methods normalization description
- **WHEN** reading the normalization description
- **THEN** it SHALL specify reference-batch ComBat with GSE100051 as the reference batch

### Requirement: Vendor and platform counts unchanged
The article SHALL continue to state four microarray platforms from four vendors (Affymetrix, Illumina, Applied Biosystems, Agilent), since dropping GSE22490 (Affymetrix) does not reduce vendor or platform count.

#### Scenario: Abstract platform statement
- **WHEN** reading the abstract
- **THEN** it SHALL state "four microarray platforms from four vendors"
