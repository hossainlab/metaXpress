# metaXpress Statistical Methods Reference

## Overview

metaXpress supports 6 meta-analysis methods, split into two families:
- **P-value combination** — combines per-study p-values into a single meta p-value
- **Effect size models** — pools log2FC estimates across studies with weighting

Choose based on the nature of your studies (see Decision Guide below).

---

## Decision Guide

```
How many studies?
├── 2 studies → Use "fisher" or "stouffer" (random effects unreliable with n=2)
└── 3+ studies →
    What is inter-study heterogeneity?
    ├── Low (I-squared < 25%) → Use "fixed_effects"
    ├── Moderate-High (I-squared >= 25%) → Use "random_effects"
    └── Unknown / Mixed → Use "awmeta" (adaptive, handles both)
```

---

## Method 1: Fisher's Combined P-value (`method = "fisher"`)

**Reference:** Rau, Marot & Jaffrézic (2013) BMC Bioinformatics

**Formula:**
```
T_Fisher = -2 * sum(log(p_i))   for i = 1..k studies

Under H0: T_Fisher ~ chi-squared(2k degrees of freedom)
```

**Implementation:**
```r
.fisher_combine <- function(pvalues) {
  # pvalues: numeric vector, one per study (NA allowed for missing)
  pvalues <- pvalues[!is.na(pvalues)]
  k <- length(pvalues)
  if (k < 2) return(NA_real_)
  T_stat <- -2 * sum(log(pvalues))
  meta_p <- pchisq(T_stat, df = 2 * k, lower.tail = FALSE)
  return(meta_p)
}
```

**When to use:** Many studies (k >= 3), independent p-values, no need for effect size estimate
**Limitation:** Does not produce a pooled effect size (log2FC). Cannot distinguish direction.

---

## Method 2: Stouffer's Z-score (`method = "stouffer"`)

**Reference:** Rau, Marot & Jaffrézic (2013) BMC Bioinformatics

**Formula:**
```
Z_i = qnorm(1 - p_i / 2) * sign(log2FC_i)   # signed Z-score per study
w_i = sqrt(n_i)                               # weight by sample size

Z_combined = sum(w_i * Z_i) / sqrt(sum(w_i^2))

Under H0: Z_combined ~ N(0, 1)
meta_p = 2 * pnorm(-|Z_combined|)
```

**Implementation:**
```r
.stouffer_combine <- function(pvalues, log2fc, n_samples) {
  valid <- !is.na(pvalues) & !is.na(log2fc) & !is.na(n_samples)
  pvalues <- pvalues[valid]; log2fc <- log2fc[valid]; n <- n_samples[valid]
  z_scores <- qnorm(1 - pvalues / 2) * sign(log2fc)
  weights  <- sqrt(n)
  z_comb   <- sum(weights * z_scores) / sqrt(sum(weights^2))
  meta_p   <- 2 * pnorm(-abs(z_comb))
  list(meta_p = meta_p, meta_z = z_comb)
}
```

**When to use:** Studies with different sample sizes; preserves direction information
**Advantage over Fisher:** Accounts for direction of effect and sample size

---

## Method 3: Fused Inverse-Normal (`method = "inverse_normal"`)

**Reference:** Prasad & Li (2021) BMC Bioinformatics

**Formula:**
```
# Transform each study's p-value to a Z-score via inverse normal
Z_i = qnorm(1 - p_i)   # one-sided, direction from sign(log2FC_i)

# Fused weight: combines sample size and effect size magnitude
w_i = sqrt(n_i) * |log2FC_i|

Z_fused = sum(w_i * Z_i) / sqrt(sum(w_i^2))
meta_p  = pnorm(-Z_fused)   # one-sided; convert to two-sided as needed
```

**When to use:** Small-n studies (n < 10 per group); outperforms Fisher/Stouffer in this regime
**Key advantage:** Upweights studies with larger effect sizes, reducing noise from underpowered studies

---

## Method 4: Fixed Effects Model (`method = "fixed_effects"`)

**Reference:** Keel & Lindholm-Perry (2022) Frontiers in Genetics

**Assumes:** True effect size is identical across all studies (homogeneous)

**Formula:**
```
# Inverse-variance weighting
w_i     = 1 / SE_i^2          # SE_i = standard error of log2FC in study i
theta_FE = sum(w_i * log2FC_i) / sum(w_i)   # pooled effect size
SE_FE    = 1 / sqrt(sum(w_i))

Z_FE    = theta_FE / SE_FE
meta_p  = 2 * pnorm(-|Z_FE|)
```

**SE estimation from DESeq2/edgeR output:**
```r
# DESeq2: lfcSE column directly available
# edgeR: SE = log2FC / qnorm(1 - pvalue/2)  [approximate]
# limma: SE = log2FC / t-statistic * sqrt(df/(df-2))  [approximate]
```

**When to use:** Low heterogeneity (I-squared < 25%); studies from very similar conditions
**Warning:** Produces anti-conservative results when heterogeneity is present

---

## Method 5: Random Effects Model (`method = "random_effects"`)

**Reference:** Keel & Lindholm-Perry (2022) Frontiers in Genetics; DerSimonian & Laird (1986)

**Assumes:** True effect sizes vary across studies (heterogeneous); drawn from a distribution

**Formula:**
```
# Step 1: Estimate between-study variance (tau-squared) via DerSimonian-Laird
Q       = sum(w_i * (log2FC_i - theta_FE)^2)   # Cochran's Q
df      = k - 1                                  # k = number of studies
C       = sum(w_i) - sum(w_i^2) / sum(w_i)
tau_sq  = max(0, (Q - df) / C)                  # between-study variance

# Step 2: Updated weights incorporating tau-squared
w_i_RE  = 1 / (SE_i^2 + tau_sq)

# Step 3: Pooled estimate
theta_RE = sum(w_i_RE * log2FC_i) / sum(w_i_RE)
SE_RE    = 1 / sqrt(sum(w_i_RE))

Z_RE    = theta_RE / SE_RE
meta_p  = 2 * pnorm(-|Z_RE|)
```

**Heterogeneity statistics:**
```
I_squared = max(0, (Q - df) / Q) * 100   # % of variance due to heterogeneity
# Interpretation: <25% low, 25-75% moderate, >75% high
```

**When to use:** Default choice for 3+ studies; robust to heterogeneity
**Minimum studies:** 3 (warn if k < 3)

---

## Method 6: Adaptive Weighting (`method = "awmeta"`)

**Reference:** Hu et al. (2025) bioRxiv — AWmeta

**Key innovation:** Adaptively combines p-value and effect size approaches based on observed heterogeneity per gene. Does not require pre-specifying fixed vs. random effects globally.

**Formula (simplified):**
```
# For each gene g:
# 1. Estimate heterogeneity (I-squared_g)
# 2. Compute adaptive weight lambda_g in [0, 1]:
#    lambda_g = 1 - I_squared_g / 100
#    (lambda=1 → pure effect size; lambda=0 → pure p-value combination)

# 3. Blend effect size Z-score and p-value Z-score:
Z_adaptive_g = lambda_g * Z_effect_g + (1 - lambda_g) * Z_pvalue_g

meta_p_g = 2 * pnorm(-|Z_adaptive_g|)
```

**When to use:** Mixed heterogeneity across genes (common in real datasets); most robust default
**Note:** This is the most recent method (2025 preprint); validate on your dataset

---

## Heterogeneity Assessment

Always compute and report heterogeneity alongside meta-analysis results.

```r
.compute_heterogeneity <- function(log2fc_vec, se_vec) {
  # log2fc_vec, se_vec: per-study vectors for one gene
  valid <- !is.na(log2fc_vec) & !is.na(se_vec)
  lfc <- log2fc_vec[valid]; se <- se_vec[valid]
  k   <- length(lfc)
  if (k < 2) return(list(Q = NA, I_sq = NA, tau_sq = NA))

  w       <- 1 / se^2
  theta_FE <- sum(w * lfc) / sum(w)
  Q       <- sum(w * (lfc - theta_FE)^2)
  df      <- k - 1
  I_sq    <- max(0, (Q - df) / Q) * 100
  C       <- sum(w) - sum(w^2) / sum(w)
  tau_sq  <- max(0, (Q - df) / C)

  list(Q = Q, df = df, I_sq = I_sq, tau_sq = tau_sq,
       p_heterogeneity = pchisq(Q, df = df, lower.tail = FALSE))
}
```

**Reporting thresholds:**
- I-squared < 25%: Low heterogeneity — fixed effects appropriate
- I-squared 25–75%: Moderate — random effects recommended
- I-squared > 75%: High — investigate sources; consider subgroup analysis

---

## Multiple Testing Correction

CRITICAL: Always apply BH (Benjamini-Hochberg) FDR correction to meta p-values.

```r
# After computing meta_pvalue for all genes:
meta_results$meta_padj <- p.adjust(meta_results$meta_pvalue, method = "BH")

# Default significance threshold:
sig_genes <- meta_results[meta_results$meta_padj <= 0.05 &
                           abs(meta_results$meta_log2FC) >= 1, ]
# Use <= and >= (inclusive inequalities)
```

---

## Method Comparison Summary

| Method | Effect Size | Direction | Handles Heterogeneity | Min Studies | Best For |
|--------|-------------|-----------|----------------------|-------------|----------|
| fisher | No | No | No | 2 | Many studies, no SE available |
| stouffer | No | Yes | No | 2 | Weighted by sample size |
| inverse_normal | No | Yes | Partial | 2 | Small-n studies |
| fixed_effects | Yes | Yes | No (assumes none) | 2 | Low heterogeneity |
| random_effects | Yes | Yes | Yes | 3 | Default for 3+ studies |
| awmeta | Yes | Yes | Yes (adaptive) | 3 | Mixed heterogeneity |
