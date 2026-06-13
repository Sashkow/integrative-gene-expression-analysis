## ADDED Requirements

### Requirement: All source-code footnotes resolve on article branch
Every `\ghref{...}` path in `main.tex` SHALL resolve to an existing file or directory on the branch referenced by `\ghb`. Any paths that do not resolve SHALL be corrected to point to the correct location.

#### Scenario: Footnote paths verified
- **WHEN** all `\ghref{path}{label}` paths are extracted from `main.tex`
- **THEN** every path exists on the `article/imputation` branch (verified via `git ls-tree` or `git show`)

#### Scenario: Broken paths corrected
- **WHEN** a `\ghref` path does not resolve on the article branch
- **THEN** the path is updated to the correct location where the script or output file exists
