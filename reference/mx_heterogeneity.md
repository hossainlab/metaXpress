# Compute per-gene heterogeneity statistics

Returns Cochran's Q, degrees of freedom, \\I^2\\, the between-study
variance \\\tau^2\\ and the p-value of the Q test for each gene measured
in at least two studies.

## Usage

``` r
mx_heterogeneity(de_results, tau_method = c("REML", "PM", "DL"))
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

- tau_method:

  Estimator of the between-study variance for `"random_effects"`:
  `"REML"` (default; restricted maximum likelihood), `"PM"`
  (Paule-Mandel) or `"DL"` (DerSimonian-Laird). REML and PM are
  preferred over DL, which underestimates \\\tau^2\\ (Veroniki et al.
  2016; Langan et al. 2019).

## Value

A `data.frame` with columns `gene_id`, `Q`, `df`, `I_sq`, `tau_sq`,
`p_heterogeneity`.

## Details

\\I^2 = \max(0, (Q - df)/Q)\\ describes the proportion of variability
due to between-study differences rather than sampling error (Higgins &
Thompson 2002). With few studies the Q test has low power and \\I^2\\ is
imprecise, so heterogeneity should be interpreted alongside \\\tau^2\\
and the prediction interval from
[`mx_meta`](https://mdjubayerhossain.com/metaXpress/reference/mx_meta.md).

## References

Cochran, W.G. (1954) The combination of estimates from different
experiments. *Biometrics*, **10**(1), 101–129.
[doi:10.2307/3001666](https://doi.org/10.2307/3001666)

Higgins, J.P.T. & Thompson, S.G. (2002) Quantifying heterogeneity in a
meta-analysis. *Statistics in Medicine*, **21**(11), 1539–1558.
[doi:10.1002/sim.1186](https://doi.org/10.1002/sim.1186)

## See also

[`mx_meta`](https://mdjubayerhossain.com/metaXpress/reference/mx_meta.md),
[`mx_heterogeneity_plot`](https://mdjubayerhossain.com/metaXpress/reference/mx_heterogeneity_plot.md)

## Examples

``` r
set.seed(1)
de <- lapply(1:3, function(s) {
  lfc <- rnorm(50, s * 0.3, 0.2); se <- rep(0.2, 50)
  p <- 2 * pnorm(-abs(lfc / se))
  data.frame(gene_id = paste0("G", 1:50), log2FC = lfc, lfcSE = se,
             pvalue = p, padj = p.adjust(p, "BH"))
})
het <- mx_heterogeneity(de)
summary(het$I_sq)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    0.00   23.13   50.66   47.10   70.70   88.32 
```
