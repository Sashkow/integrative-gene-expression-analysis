# Validation of the quantile normalization assumption for placental microarray data

## Rationale

Quantile normalization (QN) is a standard step in microarray preprocessing
pipelines such as RMA (Irizarry et al., 2003). QN forces the expression
distribution of every sample to the same shape, effectively removing
systematic differences in signal intensity between arrays. This is
appropriate when global distributional differences are technical in origin
(e.g., scanner gain, hybridisation efficiency), but it can erase genuine
biological signal if one group of samples has a systematically different
expression landscape from another.

Hicks and Irizarry (2015) formalised this concern and proposed
*quantro*, a data-driven test for whether the assumption underlying QN
holds in a given dataset. The quantro statistic compares the variability
of quantile distributions *between* groups to the variability *within*
groups. A significant result (permutation p < 0.05) indicates that the
groups have genuinely different global expression distributions, and
that applying QN may remove biologically meaningful variation.

Because our pipeline applies QN at the per-dataset preprocessing stage
(RMA for Affymetrix, analogous normalisation for Illumina and Agilent),
we tested whether the two biological contrasts present in our data ---
gestational age category (first trimester, second trimester, term) and
fetal sex (male vs female) --- exhibit global distributional differences
that QN might erase.

## Method

We ran the quantro test (Hicks and Irizarry, 2015; R/Bioconductor
package *quantro*) on pre-normalised expression data from each dataset
that contained at least two samples in each group:

- **Gestational age categories.** GSE100051 (Illumina, n = 54:
  42 first-trimester, 7 second-trimester, 5 term), GSE9984 (Affymetrix,
  n = 12: 4/4/4), GSE22490 (Affymetrix, n = 10: 8 first, 2 second).
- **Fetal sex (male vs female).** GSE100051 (34 F, 20 M), GSE9984
  (9 F, 3 M), GSE93520 (Agilent, 19 F, 17 M).

For Affymetrix datasets, CEL files were processed with `oligo::rma()`
with `normalize = FALSE` to obtain background-corrected,
median-polished expression values without quantile normalisation. For
the Illumina dataset (GSE100051), non-normalised bead-level signal
intensities were log2-transformed. For the Agilent dataset (GSE93520),
probe-level values from the series matrix were log2-transformed. Each
quantro test used 1 000 permutations to estimate the null distribution.

Datasets with only one gestational age category (GSE122214, GSE28551,
GSE37901, GSE93520) could not be tested for GA-category differences;
those with insufficient sex balance (< 2 per group after matching) were
excluded from the sex comparison.

## Results

| Dataset   | Platform                  | Contrast       | Groups                       | N per group | quantro stat | Perm. p-value |
|-----------|---------------------------|----------------|------------------------------|-------------|--------------|---------------|
| GSE100051 | Illumina HumanHT-12 V4.0 | GA categories  | 1st / 2nd / Term             | 42 / 7 / 5  | 2.79         | 0.058         |
| GSE100051 | Illumina HumanHT-12 V4.0 | Male vs Female | Female / Male                | 34 / 20     | 0.21         | 0.827         |
| GSE9984   | Affy HG-U133 Plus 2.0    | GA categories  | 1st / 2nd / Term             | 4 / 4 / 4   | 2.48         | 0.103         |
| GSE9984   | Affy HG-U133 Plus 2.0    | Male vs Female | Female / Male                | 9 / 3       | 2.66         | 0.110         |
| GSE22490  | Affy HG-U133 Plus 2.0    | GA categories  | 1st / 2nd                    | 8 / 2       | 0.64         | 0.575         |
| GSE93520  | Agilent 4x44K             | Male vs Female | Female / Male                | 19 / 17     | 0.34         | 0.705         |

No test reached statistical significance at the 0.05 level. The largest
quantro statistic was observed for GSE100051 under the GA-category
contrast (q = 2.79, p = 0.058), which is suggestive but does not cross
the significance threshold even without multiple-testing correction.
All male-vs-female comparisons yielded small quantro statistics and
large p-values, consistent with the absence of global transcriptomic
amplification or repression linked to fetal sex.

## Interpretation

These results indicate that, in the placental microarray datasets used
in our pipeline, the global shape of the expression distribution does
not differ significantly between gestational age groups or between fetal
sexes. The core assumption of quantile normalisation --- that observed
distributional differences between arrays are technical rather than
biological --- is therefore supported. Applying QN (as part of RMA or
analogous preprocessing) does not appear to erase a systematic
biological signal in these data.

The borderline result for GSE100051 (p = 0.058) warrants a note of
caution: in this dataset the first-trimester group (n = 42) vastly
outnumbers the second-trimester (n = 7) and term (n = 5) groups, and
the test may be picking up a subtle distributional shift associated
with gestational age. However, even in this case the effect is not
strong enough to reject the null at conventional thresholds, and the
practical impact on downstream differential expression is expected to
be minor given that our pipeline applies QN within each dataset before
cross-dataset integration.

## References

- Hicks, S. C. & Irizarry, R. A. (2015). quantro: a data-driven
  approach to guide the choice of an appropriate normalization method.
  *Genome Biology*, 16, 117. https://doi.org/10.1186/s13059-015-0679-0

- Irizarry, R. A., Hobbs, B., Collin, F., Beazer-Barclay, Y. D.,
  Antonellis, K. J., Scherf, U. & Speed, T. P. (2003). Exploration,
  normalization, and summaries of high density oligonucleotide array
  probe level data. *Biostatistics*, 4(2), 249--264.
  https://doi.org/10.1093/biostatistics/4.2.249

## Archive contents

```
quantro_archive/
  quantro_test.R              -- analysis script
  quantro_validation.md       -- this document
  data/
    samples.csv               -- phenodata
    raws/
      GSE100051/              -- Illumina non-normalised + series matrix
      GSE93520/               -- Agilent series matrix
      GSE9984_cel/            -- Affymetrix CEL files
      GSE22490_cel/           -- Affymetrix CEL files
  results/
    quantro_results.csv       -- raw results table
    quantro_results.xlsx      -- formatted results table
```
