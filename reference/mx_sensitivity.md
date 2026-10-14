# Leave-one-out sensitivity analysis

Reruns the meta-analysis \\k\\ times, each time excluding one study, to
assess whether any single study drives the result (Cochrane Handbook,
Section 10.14).

## Usage

``` r
mx_sensitivity(
  de_results,
  method = c("random_effects", "fixed_effects", "fisher", "stouffer", "inverse_normal",
    "minp", "maxp", "rop", "awmeta"),
  padj_threshold = 0.05,
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

- padj_threshold:

  Numeric. FDR threshold used to call genes in the robustness summary.
  Default: `0.05`.

- ...:

  Further arguments passed to
  [`mx_meta`](https://mdjubayerhossain.com/metaXpress/reference/mx_meta.md)
  (e.g. `tau_method`, `hksj`).

## Value

A named list of
[`metaXpressResult`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressResult-class.md)
objects, one per excluded study. The attribute `"robustness"` holds a
`data.frame` with, for every gene, the number and fraction of
leave-one-out runs in which it remained significant
(`meta_padj <= padj_threshold`) with the same direction as in the full
analysis.

## References

Higgins, J.P.T. et al. (eds.) (2019) *Cochrane Handbook for Systematic
Reviews of Interventions*. 2nd edn. Chichester: John Wiley & Sons.
[doi:10.1002/9781119536604](https://doi.org/10.1002/9781119536604)

## See also

[`mx_meta`](https://mdjubayerhossain.com/metaXpress/reference/mx_meta.md)

## Examples

``` r
set.seed(1)
de <- lapply(1:4, function(s) {
  lfc <- rnorm(100, c(rep(1, 10), rep(0, 90)), 0.2); se <- rep(0.25, 100)
  p <- 2 * pnorm(-abs(lfc / se))
  data.frame(gene_id = paste0("G", 1:100), log2FC = lfc, lfcSE = se,
             pvalue = p, padj = p.adjust(p, "BH"))
})
names(de) <- paste0("S", 1:4)
loo <- mx_sensitivity(de, method = "random_effects")
#> Running leave-one-out sensitivity analysis...
#>   Excluding: S1
#>   Excluding: S2
#>   Excluding: S3
#>   Excluding: S4
head(attr(loo, "robustness"))
#>   gene_id    full_padj n_loo_sig frac_loo_sig
#> 1      G1 1.480181e-14         4            1
#> 2      G2 1.167633e-15         4            1
#> 3      G3 1.562896e-11         4            1
#> 4      G4 7.823927e-16         4            1
#> 5      G5 9.692259e-08         4            1
#> 6      G6 8.707619e-17         4            1
```
