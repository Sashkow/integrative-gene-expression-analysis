## Context

The manuscript is on the `article/imputation` branch with source in `articles/imputation_article/`. The main article is `main.tex` (1945 lines). Supplements 1 and 2 are separate `.tex`/`.pdf` files in `supplement1/` and `supplement2/` subdirectories — their `.tex` sources exist on the `article/imputation` branch but not on `master`.

Key file locations:
- Main text: `articles/imputation_article/main.tex`
- Supplement 1: `articles/imputation_article/supplement1/supplement1.tex`
- Supplement 2: `articles/imputation_article/supplement2/supplement2.tex`
- Figure 2 script: `one_off_scripts/pca_before_after_combat.R` (generates `fig_pca_before_after_combat.png`)
- Figure 4 script: `scripts/integrative_analysis/article_validation/combat_sensitivity.R` (generates `post_combat_scatter.png` and `per_gene_mae_hist.png`)
- Article figures: `articles/imputation_article/figures/`

All work targets the `article/imputation` branch.

## Goals / Non-Goals

**Goals:**
- Address all 20 review findings (1 must-fix, 11 minor, 8 suggestions)
- Keep changes minimal and self-contained — each review item maps to a small, verifiable edit
- Regenerate only figures that need visual changes (Figures 2 and 4)

**Non-Goals:**
- Restructuring the article beyond what the reviewer requested
- Running the pipeline or regenerating results (all data-driven values are already correct)
- Adding new analyses or results
- Reformatting figures that the reviewer found acceptable

## Decisions

### 1. Work on the article/imputation branch directly
All edits target `article/imputation` branch files. The `.tex` sources for supplements exist there. Checkout or worktree to that branch before editing.

### 2. LaTeX text edits are direct — no programmatic generation
Review items 1, 2, 4, 5, 6, 7, 8, 10 are all direct LaTeX text edits. The outlier gene investigation (item 2) requires reading the combat sensitivity output data to identify the gene, then adding a brief note.

### 3. Table 4 → supplement via LaTeX reorganization
Move the `tab:imputation-extended` table environment from `main.tex` into `supplement1.tex` (or a new supplement). Update cross-references. This changes table numbering for Tables 5-7 — but since LaTeX auto-numbers, only `\ref{}` labels matter, not hardcoded numbers.

**Alternative considered:** Keep Table 4 in main text with reduced formatting. Rejected — reviewer specifically suggested supplementary placement, and the content is redundant with Table 3.

### 4. Figure 2 subplot titles — increase size in R script
The `pca_before_after_combat.R` script sets `plot.title` at `size = 15`. Increase to `size = 18` or `size = 20`. Regenerate the figure.

### 5. Figure 4 correlation rounding — format string change
The `combat_sensitivity.R` script uses `sprintf("...r=%.6f", overall_r)`. Change to `%.3f` to produce `r=0.999`. Regenerate.

### 6. Footnote verification — scripted check
Write a one-off script that extracts all `\ghref{...}` paths from `main.tex` and verifies each exists on the `article/imputation` branch via `git show`. Report any missing paths.

### 7. Supplement edits require checkout of article/imputation branch
The supplement `.tex` files are only on `article/imputation`. Either checkout that branch or use `git show` to read them, then write edits.

### 8. Outlier gene identification — read from existing output
The 2.41 max-deviation gene can be identified from `output/article_validation/combat_sensitivity/post_combat_comparison.csv` or by running a quick R command against the existing data. No new analysis needed.

## Risks / Trade-offs

- **Table renumbering**: Moving Table 4 to supplement changes table numbers in-text references from Table 5→4, 6→5, 7→6. LaTeX `\ref{}` handles this automatically, but any hardcoded "Table 5" in running text needs manual update. → Mitigation: grep for hardcoded table references.
- **Figure regeneration**: Changing R scripts and regenerating figures requires the R environment and packages to be installed. → Mitigation: changes are minimal (font size, format string), low risk of regression.
- **Branch complexity**: Edits span `master` (R scripts) and `article/imputation` (tex files, figures). → Mitigation: make R script edits on master, regenerate figures, copy to article branch.
- **Supplement 1 table placement**: If Table 4 moves to supplement 1, supplement 1 needs a table environment and possibly different formatting. → Mitigation: use the same tabular environment, just wrapped in the supplement's document style.
