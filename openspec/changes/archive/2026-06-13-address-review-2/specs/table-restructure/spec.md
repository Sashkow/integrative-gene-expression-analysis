## ADDED Requirements

### Requirement: Table 4 moved to supplementary material
The extended block-mask cross-validation table (currently Table 4, `\label{tab:imputation-extended}`) SHALL be moved from `main.tex` to a supplement. The main text SHALL retain a reference to the supplementary table. All in-text references to `tab:imputation-extended` SHALL be updated to point to the supplement.

#### Scenario: Table no longer in main text
- **WHEN** the main article PDF is compiled
- **THEN** the extended block-mask table does not appear in the main body; it appears in a supplement

#### Scenario: Table numbering updated
- **WHEN** the main article PDF is compiled
- **THEN** the remaining main-text tables (former Tables 5, 6, 7) are renumbered automatically via LaTeX (now Tables 4, 5, 6) and all `\ref{}` cross-references resolve correctly

#### Scenario: Main text retains reference
- **WHEN** the reader reaches the imputation results discussion
- **THEN** the text references the supplementary table for the extended block-mask results
