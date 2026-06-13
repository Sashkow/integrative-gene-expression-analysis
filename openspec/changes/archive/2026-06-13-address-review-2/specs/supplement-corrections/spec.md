## ADDED Requirements

### Requirement: Lanoix citation author duplication fixed (S2.1, must-fix)
Supplement 2 (`supplement2.tex`) SHALL correct the Lanoix citation (ref 15) to remove duplicated "Bhatt H" author entries. The corrected citation SHALL match the actual publication's author list.

#### Scenario: Citation verified
- **WHEN** the Lanoix reference is checked in the supplement 2 bibliography
- **THEN** no author name appears more than once

### Requirement: Supplement 1 notes ref-batch ComBat adoption (S1.1)
Supplement 1 SHALL note that the main article adopted recommendation #5 (reference-batch ComBat) and reference the ComBat sensitivity result (r=0.9991) as empirical support.

#### Scenario: Recommendation adoption cross-referenced
- **WHEN** the reader reads supplement 1's recommendation ladder
- **THEN** recommendation #5 includes a note that the main article implemented this approach with measured impact r=0.9991

### Requirement: HarmonizR acknowledged as future work (S1.2)
Supplement 1 SHALL acknowledge that HarmonizR is described as "directly applicable" but was not benchmarked in this study, and note it as future work.

#### Scenario: HarmonizR limitation noted
- **WHEN** the reader reads the HarmonizR description in supplement 1
- **THEN** there is an explicit note that it was not benchmarked and is a candidate for future evaluation

### Requirement: M1 vs M2 imputation choice justified (S1.3)
Supplement 1 SHALL justify why the pipeline uses global imputation (M1) rather than batch-sensitized imputation (M2), or explicitly acknowledge this as a limitation.

#### Scenario: Imputation strategy choice addressed
- **WHEN** the reader reads supplement 1's discussion of M1 vs M2
- **THEN** they find either a justification for the M1 choice or an acknowledgment that M2 was not tested

### Requirement: Unused GEO citations removed from supplement 1 (S1.4)
Supplement 1 SHALL remove GEO dataset citations (e.g., GSE6573, GSE73374) that are not used in the main article, to reduce reader confusion.

#### Scenario: No dangling dataset citations
- **WHEN** GEO accession numbers in supplement 1 are checked against the main article
- **THEN** all cited accessions either appear in the main article or are clearly marked as illustrative examples

### Requirement: Supplement 2 dataset count matches main article (S2.2)
Supplement 2 references to "14 datasets" and "227 ComBat DEGs for 1_2 balanced" SHALL be updated to match the current 6-dataset configuration and the main article's numbers, or an explanatory note SHALL be added.

#### Scenario: Numbers consistent
- **WHEN** the reader compares supplement 2's dataset/DEG counts with the main article
- **THEN** the numbers match or the discrepancy is explicitly explained

### Requirement: RUV improvements noted as planned (S2.3)
The "Possible improvements" section in supplement 2 SHALL note which of the 6 listed strategies are planned for follow-up work.

#### Scenario: Follow-up plans noted
- **WHEN** the reader reads the possible improvements section
- **THEN** at least a brief note indicates which strategies are under consideration for follow-up

### Requirement: RUVinv 10 DEGs identified (S2.4)
Supplement 2 SHALL briefly identify or characterize the 10 DEGs produced by RUVinv for the 1_2 contrast (e.g., list them or note that they represent the highest-confidence DEGs).

#### Scenario: RUVinv genes described
- **WHEN** the reader reads the RUVinv results
- **THEN** they find a brief characterization of which 10 genes were identified
