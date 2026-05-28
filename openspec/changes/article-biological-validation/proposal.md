## Why

The imputation article (BMC Bioinformatics target) identifies 380 DEGs between first- and second-trimester placenta, of which 198 are "gained" by softImpute imputation. The central reviewer concern will be: "Are these 198 genes biologically real or imputation artefacts?" The article currently has no biological validation — no enrichment, no literature concordance, no cross-tissue test. This change adds four layers of biological validation and restructures the article to present them.

A preliminary concordance check during exploration found that 90.1% of gained DEGs show concordant logFC direction with the independent 2021 Lykhenko analysis (4 Affy datasets, 22 samples), and 67% were already significant at FDR < 0.05 in that analysis — lost only because the intersection step discarded them when non-Affy platforms were added. This is strong evidence against the artefact hypothesis.

## What Changes

### Layer 1: Article layout restructuring (main.tex)

Add three new Results subsections after the existing Differential Expression section (§4):

- **§4a GO/KEGG enrichment of DEGs** — Functional enrichment analysis of the full 380, intersection-only 182, and gained 198 DEGs. Table of top enriched terms per group + dot plots.

- **§4b Literature concordance with published placental studies** — For each DEG with published placental qPCR/RT-PCR/IHC/Western data, report published direction vs. pipeline logFC. Concordance summary table modeled on the sysbio_as_imbg_methodological_hub article's Table 4 (prostate qPCR comparison). Target: 30–80 genes with published data, expect ≥70% concordance.

- **§4c Prior-study concordance of gained DEGs** — The strongest validation: compare the 201 gained DEGs against the full limma table from Lykhenko 2021a (20,162 genes, `output/dissertation/integrative_pipeline_disser/1_2/difexp/difexp_all.csv`). Exploration already confirmed: 181/201 testable, 163/181 (90.1%) same direction, 122 already FDR < 0.05 in 2021. Add logFC scatter plot, volcano plot with gained DEGs highlighted, and categorization table.

- **§4d (optional) Cross-tissue pipeline validation: prostate cancer** — Brief subsection showing the pipeline applied to 3 prostate cancer GEO datasets (from sysbio article Dept 1), with qPCR concordance from 4 department publications (77 genes, 52% confirmed, 69% direction concordance). Demonstrates pipeline generalizability. Decision needed: whether this can appear here or is reserved for the sysbio paper.

Add corresponding Discussion paragraphs interpreting enrichment, literature concordance, and prior-study concordance results.

### Layer 2: Literature wet-results mining

Systematic PubMed search for published placental gene expression data matching our 380 DEGs:

- Search strategy per gene: `"{GENE_SYMBOL} placenta trimester expression qPCR OR RT-PCR OR immunohistochemistry OR Western"`
- Also check Human Protein Atlas (placenta tissue page) for protein-level expression by trimester
- Low-hanging fruit: well-known placental markers (LEP, FLT1, ENG, PAPPA, ADAM12, HSD3B1, CGA, CGB, GH2, PSG family)
- Special attention to the 198 gained DEGs — even 10–15 with literature confirmation directly addresses reviewer concern
- Output: curated CSV with gene, PMID, method, published direction, our logFC, concordant yes/no

### Layer 3: Coding tasks

**A. GO/KEGG enrichment script**
- Input: `difexp_significant_softimpute_combat_ref.tsv` (380 DEGs) + `difexp_significant_none_combat_ref.tsv` (182 intersection DEGs)
- Tool: `clusterProfiler` (enrichGO, enrichKEGG)
- Compute enrichment for: full 380, gained ~198, intersection-only 182
- Output: enrichment tables, dot plots, comparison figure
- Key question: do gained DEGs enrich for known placental development pathways?

**B. Prior-study concordance analysis script**
- Input: 201 gained DEG IDs, `output/dissertation/integrative_pipeline_disser/1_2/difexp/difexp_all.csv` (Lykhenko 2021a full limma)
- Analysis:
  1. Match gained DEGs to 2021 table (181/201 = 90% match confirmed)
  2. Categorize each: same-direction-significant (122), same-direction-near-sig (18), same-direction-weak (23), opposite (18), absent (20)
  3. logFC scatter plot: this study vs 2021
  4. Volcano plot of 2021 data with gained DEGs highlighted in colour
  5. Summary statistics and concordance rate
- Output: figures + concordance table for paper

**C. Literature concordance assembly script**
- Input: curated CSV from Layer 2 (gene + published direction from literature)
- Compare against: logFC from DEG table
- Output: concordance table (like sysbio article Table 4)

**D. Venn/upset diagram**
- Intersection-only DEGs (182) vs. softImpute DEGs (380) vs. Lykhenko 2021a DEGs (327)
- Three-way Venn showing shared/unique genes across analyses
- Highlights that the "gained" genes largely overlap with 2021 significant genes

**E. Mean-imputation benchmark**
- Add per-gene-mean and per-batch-mean imputation to the cross-validation
- Run alongside softImpute under both random-cell and block masking
- Update Table 5 with baseline comparators

**F. ComBat sensitivity analysis**
- Compare ComBat gamma/delta batch parameters between intersection-only (8,260 genes) and full imputed (17,531 genes) matrices
- Report correlation and max deviation
- Output: correlation plot + summary stats

### Layer 4: Wet-lab experiment suggestions (for future / if feasible)

- Targeted qPCR panel of ~20 top DEGs (mix of gained + intersection) on whatever 2nd-trimester placenta samples are accessible
- Existing biobanks or collaborator tissue as sample source
- Note: wet-lab validation is LOW priority for BMC Bioinformatics (methods journal) — enrichment + literature concordance + prior-study concordance should suffice for acceptance

## Capabilities

### New Capabilities
- `deg-enrichment-analysis`: GO/KEGG enrichment pipeline for DEG lists with comparison across pipeline variants
- `prior-study-concordance`: Script to compare gained DEGs against historical full limma output and produce concordance figures
- `literature-concordance`: Curated database of published placental gene expression results with automated comparison to pipeline output

### Modified Capabilities
- `article-results-section`: Three new subsections added to Results (enrichment, literature concordance, prior-study concordance)
- `article-discussion`: New paragraphs interpreting biological validation findings

## Impact

### Files to modify
- `articles/imputation_article/main.tex` — add Results subsections §4a–§4c (or §4d), Discussion paragraphs, new figure/table references
- `articles/imputation_article/supplement1/supplement1.tex` — may reference enrichment details
- Pipeline output used as input: `output/yehor_sashko/phase2b_1_2_yehor_7ds_no_37653_enriched_sashko/` (current DEGs), `output/dissertation/integrative_pipeline_disser/1_2/difexp/difexp_all.csv` (2021 full limma)

### New files
- R script(s) for enrichment analysis, concordance analysis, figures
- Curated CSV of literature-mined placental gene expression results
- New figures: enrichment dot plots, logFC scatter, volcano, Venn/upset diagram

### Key data confirmed during exploration
- 380 DEGs in `difexp_significant_softimpute_combat_ref.tsv` (7ds, softImpute + ComBat-ref)
- 182 DEGs in `difexp_significant_none_combat_ref.tsv` (7ds, intersection + ComBat-ref)
- 201 gained DEGs (380 − 182, accounting for a few genes significant only in intersection)
- 20,162 genes in 2021 full limma: `output/dissertation/integrative_pipeline_disser/1_2/difexp/difexp_all.csv`
- 181/201 gained DEGs found in 2021 table (90%)
- 163/181 same direction (90.1%), 122 already FDR < 0.05 in 2021 (67%)
- Sysbio prostate validation already written: 77 genes, 52% confirmed, 69% direction concordance
