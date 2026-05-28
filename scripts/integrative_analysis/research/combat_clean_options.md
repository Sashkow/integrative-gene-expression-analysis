# Batch Correction Strategy for Partial Genes in "Clean ComBat" Pipeline

## Problem

When using softImpute + ComBat + limma:
- ComBat runs on the full matrix including imputed values
- Imputed values contaminate ComBat's batch parameter estimation
- This shifts "real" values before limma sees them
- Setting weight=0 in limma for imputed values can't undo the upstream ComBat corruption

**Solution**: Run ComBat on inner-join genes only (clean estimation). But what about genes that are partially observed (present in some but not all studies)?

## Three Options for Partial Genes

### Option A: Per-batch median shift (simplest)

**How it works:**
1. Run ComBat on inner-join genes → get corrected matrix
2. For each batch, compute the median per-gene shift: `median(corrected[g,s] - original[g,s])` across all inner-join genes `g` for samples `s` in that batch
3. Apply that median shift to the real observations of partial genes in each batch
4. Leave imputed observations unchanged (weight=0 in limma)

**Pros:**
- Very simple to implement (~10 lines of R)
- No dependency on sva internals
- Robust (median, not mean)
- Captures the dominant batch correction effect (location shift)

**Cons:**
- Only corrects location (mean shift), not scale (variance)
- Assumes all genes in a batch experience the same shift — ignores gene-specific batch effects
- ComBat actually does per-gene corrections; this collapses that to one number per batch

**When this is sufficient:**
- When batch effects are primarily additive (location shifts)
- When the partial genes behave similarly to inner-join genes w.r.t. batch

### Option B: Full EB with reused priors (most principled)

**How it works:**
1. Run ComBat on inner-join genes
2. Extract the EB prior hyperparameters estimated across inner-join genes:
   - `gamma.bar[j]` — prior mean for batch location (per batch j)
   - `t2[j]` — prior variance for batch location
   - `a.prior[j]`, `b.prior[j]` — InvGamma params for batch scale
3. For each partial gene:
   - Compute gene-specific `var.pooled`, `stand.mean` using only real observations
   - Standardize real observations
   - Estimate per-batch `gamma.hat[j,g]` and `delta.hat[j,g]` from available data
   - Shrink using the saved priors: `gamma.star = postmean(gamma.hat, gamma.bar, n, delta, t2)`
   - Apply the correction formula to real observations only
4. Leave imputed observations unchanged

**Pros:**
- Statistically principled — same EB framework as ComBat
- Per-gene corrections (both location and scale)
- Priors are estimated from clean data (inner-join genes)
- For partial genes with few observations per batch, shrinkage toward the prior is strong (conservative)

**Cons:**
- Complex implementation (~50-80 lines of R)
- Requires accessing `sva:::` internals or reimplementing ComBat's EB formulas
- May break if sva package internals change across versions
- For partial genes with very few observations in a batch (e.g., 1-2 samples), the per-gene estimates are very noisy and the result is dominated by the prior anyway — collapsing to something similar to Option A

**Key ComBat formulas needed:**
```r
# Standardize
stand.mean[g,s] = grand.mean[g] + X[s] * beta[g]
s.data[g,s] = (dat[g,s] - stand.mean[g,s]) / sqrt(var.pooled[g])

# Naive batch estimates (on standardized data)
gamma.hat[j,g] = mean(s.data[g, samples_in_batch_j])
delta.hat[j,g] = var(s.data[g, samples_in_batch_j])

# EB priors (estimated from inner-join genes, REUSED for partial genes)
gamma.bar[j] = mean(gamma.hat[j,]) across inner-join genes
t2[j] = var(gamma.hat[j,]) across inner-join genes
a.prior[j], b.prior[j] = method-of-moments InvGamma fit to delta.hat[j,]

# EB shrinkage (iterative conditional modes)
gamma.star[j,g] = (t2[j]*n_j*gamma.hat[j,g] + delta.star*gamma.bar[j]) / (t2[j]*n_j + delta.star)
delta.star[j,g] = (0.5*sum2 + b.prior[j]) / (n_j/2 + a.prior[j] - 1)

# Correction
corrected[g,s] = (s.data[g,s] - gamma.star[j,g]) / sqrt(delta.star[j,g]) * sqrt(var.pooled[g]) + stand.mean[g,s]
```

### Option C: No correction for partial genes

**How it works:**
1. Run ComBat on inner-join genes only
2. Partial genes pass through to limma uncorrected
3. Imputed values get weight=0 in limma; real observations get weight=1
4. No batch correction on partial genes' real observations

**Pros:**
- Simplest possible implementation
- No risk of incorrect batch correction on partial genes
- Clean separation: well-corrected inner-join genes vs. uncorrected partial genes

**Cons:**
- Partial genes retain batch effects, which inflate variance and reduce power
- Creates a systematic difference between inner-join and partial genes (corrected vs. uncorrected scale)
- limma's eBayes variance shrinkage uses all genes — mixing corrected and uncorrected genes distorts the prior
- Partial genes may produce false positives if a batch effect aligns with the biological contrast

**When this might be OK:**
- When batch effects are small
- When you mainly care about inner-join gene results and partial genes are "bonus"

## Recommendation

Start with **Option A** (per-batch median shift) as the default. It's simple, captures the main correction, and avoids the complexity risks of Option B. If results look suspicious (e.g., partial-gene DEGs show batch-correlated patterns), upgrade to Option B.

Option C is too risky — uncorrected batch effects on partial genes will contaminate limma's global variance estimation.

## Implementation note

The infrastructure (na_mask tracking, weighted limma, ComBat on inner-join) is the same for all three options. The only difference is what happens to partial genes in the normalization step. This can be made into a parameter: `partial_gene_correction = c("median_shift", "eb_priors", "none")`.
