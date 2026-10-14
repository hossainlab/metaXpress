# Run meta-analysis across per-study DE results

Integrates per-study differential expression results from
[`mx_de_all`](https://mdjubayerhossain.com/metaXpress/reference/mx_de_all.md)
into a single gene-level meta-analysis. Two families of methods are
available: *effect-size* models, which pool the per-study log2 fold
changes weighted by their sampling variance, and *p-value combination*
methods, which combine per-study evidence without pooling effects.

## Usage

``` r
mx_meta(
  de_results,
  method = c("random_effects", "fixed_effects", "fisher", "stouffer", "inverse_normal",
    "minp", "maxp", "rop", "awmeta"),
  min_studies = 2,
  tau_method = c("REML", "PM", "DL"),
  hksj = FALSE,
  directional = TRUE,
  weights = NULL,
  r = NULL,
  alpha = 0.05,
  n_samples = NULL,
  ...
)
```

## Arguments

- de_results:

  A named list of `data.frame` objects (the `de_result` slots produced
  by
  [`mx_de`](https://mdjubayerhossain.com/metaXpress/reference/mx_de.md)),
  or a named list of
  [`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
  objects after
  [`mx_de_all`](https://mdjubayerhossain.com/metaXpress/reference/mx_de_all.md).
  Required columns: `gene_id`, `log2FC`, `pvalue`, `padj`; effect-size
  methods also need `lfcSE`.

- method:

  Character scalar. One of:

  `"random_effects"`

  :   (default) Inverse-variance random-effects model; between-study
      variance \\\tau^2\\ estimated by `tau_method`. Estimates the
      *average* effect across populations and reports a prediction
      interval.

  `"fixed_effects"`

  :   Inverse-variance common-effect model. Appropriate only when
      studies are near-replicates (low \\I^2\\).

  `"fisher"`

  :   Fisher's sum of \\-2\log p\\. With `directional = TRUE` (default)
      the two one-sided Fisher tests are computed and combined as in
      Owen (2009), so evidence in opposite directions does not
      accumulate. Like all Fisher-type tests it can still be driven by
      one extreme study; check `direction_consistency`.

  `"stouffer"`

  :   Unweighted Stouffer Z on signed z-scores.

  `"inverse_normal"`

  :   Weighted inverse-normal (Stouffer) combination with weights
      \\\sqrt{n_s}\\, as in Marot et al. (2009) and Rau et al. (2014).

  `"minp"`

  :   Tippett's minimum p: DE in *at least one* study. Non-directional.

  `"maxp"`

  :   Wilkinson's maximum p: DE in *all* studies. Non-directional.

  `"rop"`

  :   r-th ordered p-value (Song & Tseng 2014): DE in at least `r`
      studies. Non-directional.

  `"awmeta"`

  :   AWmeta (Hu et al. 2026). The p-value is the adaptively weighted
      Fisher statistic (AW-Fisher; Li & Tseng 2011), calibrated with the
      AWFisher package (Huo et al. 2020), which also returns binary
      per-study weights indicating which studies contribute. The effect
      size is an REML random-effects estimate over the selected studies
      only (AW-REM). Requires AWFisher and `lfcSE`. Non-directional, as
      in the original method; check `direction_consistency`. The
      selected studies are reported in `aw_studies`.

- min_studies:

  Integer. Minimum number of studies in which a gene must be measured.
  Default: `2`.

- tau_method:

  Estimator of the between-study variance for `"random_effects"`:
  `"REML"` (default; restricted maximum likelihood), `"PM"`
  (Paule-Mandel) or `"DL"` (DerSimonian-Laird). REML and PM are
  preferred over DL, which underestimates \\\tau^2\\ (Veroniki et al.
  2016; Langan et al. 2019).

- hksj:

  Logical. Apply the Hartung-Knapp-Sidik-Jonkman adjustment to the
  random-effects standard error and use a \\t\_{k-1}\\ reference
  distribution. Recommended when few studies are pooled because the
  standard Wald test is anti-conservative (IntHout et al. 2014).
  Default: `FALSE`.

- directional:

  Logical. For `"fisher"`, combine one-sided p-values in each direction
  (default `TRUE`). `FALSE` gives the classical Fisher test on two-sided
  p-values, which ignores sign.

- weights:

  For `"stouffer"`-type methods and the effect estimate reported by
  p-value methods: `NULL` (default) uses per-study sample sizes from the
  `n_case`/`n_control` columns written by
  [`mx_de`](https://mdjubayerhossain.com/metaXpress/reference/mx_de.md)
  for `"inverse_normal"` and equal weights for `"stouffer"`. A numeric
  vector (one per study) overrides this.

- r:

  Integer. Order statistic for `"rop"`. Default: `ceiling(k / 2)` for a
  gene measured in `k` studies.

- alpha:

  Numeric. Confidence level is `1 - alpha` for confidence and prediction
  intervals. Default: `0.05`.

- n_samples:

  Deprecated alias for `weights` supplied as total sample sizes; its
  square root is used.

- ...:

  Unused.

## Value

A
[`metaXpressResult`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressResult-class.md).
Its `meta_table` has one row per gene with columns:

- gene_id:

  Gene identifier.

- meta_log2FC, meta_se:

  Pooled log2 fold change and its standard error. For p-value methods
  this is the inverse-variance common-effect estimate, reported as a
  descriptive effect size.

- ci_lower, ci_upper:

  `1 - alpha` confidence interval.

- pi_lower, pi_upper:

  Random-effects prediction interval for the effect in a new study
  (Higgins et al. 2009); requires \\k \ge 3\\.

- meta_stat:

  Test statistic (Z, t or \\\chi^2\\ as appropriate).

- meta_pvalue, meta_padj:

  Combined p-value and its Benjamini-Hochberg adjustment across genes.

- tau2, i_squared, q_stat, q_pvalue:

  Between-study variance, \\I^2\\ (%), Cochran's Q and its p-value.

- n_studies, n_up, n_down:

  Studies measuring the gene, and those with positive / negative log2FC.

- direction_consistency:

  Fraction of studies whose log2FC sign agrees with `meta_log2FC`.

## Details

Use `meta_padj` (BH-adjusted) for significance. For a typical
disease-signature meta-analysis we recommend `"random_effects"` with
`tau_method = "REML"`, and reporting \\I^2\\ and the prediction interval
alongside the pooled effect. P-value combination methods answer
different questions (any vs. all vs. most studies) and should be chosen
a priori, not by comparing numbers of significant genes.

A study contributes to a gene only when it reports `log2FC` and `pvalue`
(plus `lfcSE` for effect-size methods). DESeq2 sets `pvalue` to `NA` for
count outliers flagged by Cook's distance, so such genes are treated as
not measured in that study. (DESeq2's independent filtering affects only
`padj`, which is not used.)

All methods assume studies are independent. Fold-change direction must
be harmonised across studies (same case/control orientation); check with
[`mx_study_concordance`](https://mdjubayerhossain.com/metaXpress/reference/mx_study_concordance.md).

## References

DerSimonian, R. & Laird, N. (1986) Meta-analysis in clinical trials.
*Controlled Clinical Trials*, **7**(3), 177–188.
[doi:10.1016/0197-2456(86)90046-2](https://doi.org/10.1016/0197-2456%2886%2990046-2)

Viechtbauer, W. (2005) Bias and efficiency of meta-analytic variance
estimators in the random-effects model. *Journal of Educational and
Behavioral Statistics*, **30**(3), 261–293.
[doi:10.3102/10769986030003261](https://doi.org/10.3102/10769986030003261)

Paule, R.C. & Mandel, J. (1982) Consensus values and weighting factors.
*Journal of Research of the National Bureau of Standards*, **87**(5),
377–385.
[doi:10.6028/jres.087.022](https://doi.org/10.6028/jres.087.022)

Hartung, J. & Knapp, G. (2001) A refined method for the meta-analysis of
controlled clinical trials with binary outcome. *Statistics in
Medicine*, **20**(24), 3875–3889.
[doi:10.1002/sim.1009](https://doi.org/10.1002/sim.1009)

IntHout, J., Ioannidis, J.P.A. & Borm, G.F. (2014) The
Hartung-Knapp-Sidik-Jonkman method for random effects meta-analysis is
straightforward and considerably outperforms the standard
DerSimonian-Laird method. *BMC Medical Research Methodology*, **14**,
25.
[doi:10.1186/1471-2288-14-25](https://doi.org/10.1186/1471-2288-14-25)

Higgins, J.P.T., Thompson, S.G. & Spiegelhalter, D.J. (2009) A
re-evaluation of random-effects meta-analysis. *Journal of the Royal
Statistical Society: Series A*, **172**(1), 137–159.
[doi:10.1111/j.1467-985X.2008.00552.x](https://doi.org/10.1111/j.1467-985X.2008.00552.x)

Higgins, J.P.T. & Thompson, S.G. (2002) Quantifying heterogeneity in a
meta-analysis. *Statistics in Medicine*, **21**(11), 1539–1558.
[doi:10.1002/sim.1186](https://doi.org/10.1002/sim.1186)

Veroniki, A.A. et al. (2016) Methods to estimate the between-study
variance and its uncertainty in meta-analysis. *Research Synthesis
Methods*, **7**(1), 55–79.
[doi:10.1002/jrsm.1164](https://doi.org/10.1002/jrsm.1164)

Langan, D. et al. (2019) A comparison of heterogeneity variance
estimators in simulated random-effects meta-analyses. *Research
Synthesis Methods*, **10**(1), 83–98.
[doi:10.1002/jrsm.1316](https://doi.org/10.1002/jrsm.1316)

Owen, A.B. (2009) Karl Pearson's meta-analysis revisited. *Annals of
Statistics*, **37**(6B), 3867–3892.
[doi:10.1214/09-AOS697](https://doi.org/10.1214/09-AOS697)

Marot, G., Foulley, J.-L., Mayer, C.-D. & Jaffrezic, F. (2009) Moderated
effect size and P-value combinations for microarray meta-analyses.
*Bioinformatics*, **25**(20), 2692–2699.
[doi:10.1093/bioinformatics/btp444](https://doi.org/10.1093/bioinformatics/btp444)

Rau, A., Marot, G. & Jaffrezic, F. (2014) Differential meta-analysis of
RNA-seq data from multiple studies. *BMC Bioinformatics*, **15**, 91.
[doi:10.1186/1471-2105-15-91](https://doi.org/10.1186/1471-2105-15-91)

Song, C. & Tseng, G.C. (2014) Hypothesis setting and order statistic for
robust genomic meta-analysis. *Annals of Applied Statistics*, **8**(2),
777–800. [doi:10.1214/13-AOAS683](https://doi.org/10.1214/13-AOAS683)

Li, J. & Tseng, G.C. (2011) An adaptively weighted statistic for
detecting differential gene expression when combining multiple
transcriptomic studies. *Annals of Applied Statistics*, **5**(2A),
994–1019. [doi:10.1214/10-AOAS393](https://doi.org/10.1214/10-AOAS393)

Huo, Z., Tang, S., Park, Y. & Tseng, G. (2020) P-value evaluation,
variability index and biomarker categorization for adaptively weighted
Fisher's meta-analysis method in omics applications. *Bioinformatics*,
**36**(2), 524–532.
[doi:10.1093/bioinformatics/btz589](https://doi.org/10.1093/bioinformatics/btz589)

Hu, Y. et al. (2026) AWmeta empowers adaptively weighted transcriptomic
meta-analysis. *Current Issues in Molecular Biology*, **48**(5), 530.
[doi:10.3390/cimb48050530](https://doi.org/10.3390/cimb48050530)

## See also

[`mx_heterogeneity`](https://mdjubayerhossain.com/metaXpress/reference/mx_heterogeneity.md),
[`mx_sensitivity`](https://mdjubayerhossain.com/metaXpress/reference/mx_sensitivity.md),
[`mx_study_concordance`](https://mdjubayerhossain.com/metaXpress/reference/mx_study_concordance.md),
[`mx_forest`](https://mdjubayerhossain.com/metaXpress/reference/mx_forest.md)

## Examples

``` r
# Three per-study results in the format written by mx_de()
set.seed(1)
genes <- paste0("G", 1:200)
de <- lapply(1:3, function(s) {
  lfc <- rnorm(200, c(rep(1, 20), rep(0, 180)), 0.3)
  se  <- runif(200, 0.2, 0.4)
  p   <- 2 * pnorm(-abs(lfc / se))
  data.frame(gene_id = genes, log2FC = lfc, lfcSE = se, pvalue = p,
             padj = p.adjust(p, "BH"), n_case = 5, n_control = 5)
})
names(de) <- c("S1", "S2", "S3")
res <- mx_meta(de, method = "random_effects", tau_method = "REML")
#> Running meta-analysis (random_effects) on 3 studies...
head(res@meta_table[order(res@meta_table$meta_pvalue), ])
#>    gene_id meta_log2FC   meta_se  ci_lower ci_upper  pi_lower pi_upper
#> 5       G5   1.2844615 0.1879690 0.9160490 1.652874 -1.543210 4.112133
#> 12     G12   1.0837146 0.1735513 0.7435604 1.423869 -1.121463 3.288893
#> 10     G10   1.0446984 0.1680409 0.7153443 1.374053 -1.090464 3.179861
#> 8       G8   1.0759511 0.1752744 0.7324196 1.419483 -1.151122 3.303024
#> 1       G1   0.9789829 0.1600478 0.6652951 1.292671 -1.054617 3.012583
#> 11     G11   1.2198264 0.1994541 0.8289036 1.610749 -1.314478 3.754131
#>    meta_stat  meta_pvalue    meta_padj       tau2 i_squared     q_stat
#> 5   6.833369 8.294346e-12 1.658869e-09 0.01419286         0 1.89583209
#> 12  6.244349 4.255700e-10 3.077624e-08 0.00000000         0 0.04912692
#> 10  6.216929 5.069799e-10 3.077624e-08 0.00000000         0 0.15762653
#> 8   6.138665 8.321762e-10 3.077624e-08 0.00000000         0 0.55516997
#> 1   6.116817 9.546260e-10 3.077624e-08 0.00000000         0 1.27237086
#> 11  6.115825 9.605844e-10 3.077624e-08 0.00000000         0 0.75678390
#>     q_pvalue n_studies n_up n_down direction_consistency
#> 5  0.3875478         3    3      0                     1
#> 12 0.9757358         3    3      0                     1
#> 10 0.9242125         3    3      0                     1
#> 8  0.7576112         3    3      0                     1
#> 1  0.5293077         3    3      0                     1
#> 11 0.6849620         3    3      0                     1
```
