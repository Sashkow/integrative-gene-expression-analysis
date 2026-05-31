## ADDED Requirements

### Requirement: Supplement 1 updated to 6ds
Supplement 1 (`supplement1/supplement1.tex`) SHALL be reviewed and updated to reflect 6ds results. Any dataset listings, sample counts, or result numbers SHALL match the 6ds run.

#### Scenario: Supplement 1 dataset references
- **WHEN** reading supplement1.tex
- **THEN** all dataset references SHALL list 6 datasets without GSE22490

### Requirement: Supplement 2 updated to 6ds
Supplement 2 (`supplement2/supplement2.tex`) SHALL be reviewed and updated to reflect 6ds results.

#### Scenario: Supplement 2 consistency
- **WHEN** reading supplement2.tex
- **THEN** all numerical results and dataset references SHALL match the 6ds run

### Requirement: Data availability accessions updated
The data availability section in main.tex and any supplementary materials SHALL list exactly 6 GEO accessions: GSE100051, GSE122214, GSE28551, GSE37901, GSE93520, GSE9984.

#### Scenario: Accession list in data availability
- **WHEN** reading the data availability statement
- **THEN** it SHALL list exactly 6 accessions with no GSE22490
