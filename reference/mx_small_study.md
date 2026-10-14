# Small-study effects and publication bias tests

Tests, gene by gene, whether smaller (less precise) studies report
systematically larger effects, the pattern produced by publication bias
or by other small-study effects. Two standard tests are provided:
Egger's regression test and Begg and Mazumdar's rank correlation test.

## Usage

``` r
mx_small_study(
  de_results,
  genes = NULL,
  method = c("egger", "begg", "both"),
  min_studies = 3
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

- genes:

  Optional character vector of genes to test (e.g. the significant genes
  of
  [`mx_meta`](https://mdjubayerhossain.com/metaXpress/reference/mx_meta.md)).
  Default: all genes measured in at least `min_studies` studies.

- method:

  `"egger"` (default), `"begg"` or `"both"`.

- min_studies:

  Minimum number of studies with a standard error for a gene to be
  tested. Default `3` (the minimum for Egger's test).

## Value

A `data.frame` with `gene_id`, `n_studies` and, for Egger,
`egger_intercept`, `egger_se`, `egger_pvalue`; for Begg, `begg_tau`,
`begg_pvalue`. BH-adjusted columns (`*_padj`) are added across the genes
tested.

## Details

**Egger's test** regresses the standardized effect \\y_i / s_i\\ on
precision \\1/s_i\\ by ordinary least squares; an intercept different
from zero indicates funnel-plot asymmetry (t-test with \\k - 2\\ df).
**Begg's test** computes Kendall's \\\tau\\ between the standardized
deviates from the fixed-effect estimate and their variances.

Both tests have low power with few studies. The Cochrane Handbook
advises against them with fewer than 10 studies, and a warning is issued
in that case; a non-significant result then is not evidence of absence
of bias. In transcriptomic meta-analyses of public data, genome-wide
summary statistics are reported for all genes, so classical publication
bias (selective reporting of genes) does not apply; asymmetry instead
points to small-study differences in design, tissue or platform.

## References

Egger, M., Davey Smith, G., Schneider, M. & Minder, C. (1997) Bias in
meta-analysis detected by a simple, graphical test. *BMJ*,
**315**(7109), 629–634.
[doi:10.1136/bmj.315.7109.629](https://doi.org/10.1136/bmj.315.7109.629)

Begg, C.B. & Mazumdar, M. (1994) Operating characteristics of a rank
correlation test for publication bias. *Biometrics*, **50**(4),
1088–1101. [doi:10.2307/2533446](https://doi.org/10.2307/2533446)

Sterne, J.A.C. et al. (2011) Recommendations for examining and
interpreting funnel plot asymmetry in meta-analyses of randomised
controlled trials. *BMJ*, **343**, d4002.
[doi:10.1136/bmj.d4002](https://doi.org/10.1136/bmj.d4002)

## See also

[`mx_funnel`](https://mdjubayerhossain.com/metaXpress/reference/mx_funnel.md),
[`mx_meta`](https://mdjubayerhossain.com/metaXpress/reference/mx_meta.md)

## Examples

``` r
ad <- readRDS(system.file("extdata", "ad_brain_de.rds",
                          package = "metaXpress"))
ss <- suppressWarnings(mx_small_study(ad$de_results,
                                      genes = c("GFAP", "VGF", "SST")))
ss
#>   gene_id n_studies egger_intercept  egger_se egger_pvalue egger_padj
#> 1     SST         4       2.3710865 3.5614876    0.5740746  0.8611119
#> 2     VGF         4      -0.2196817 4.1596622    0.9626820  0.9626820
#> 3    GFAP         4       1.0189531 0.5046195    0.1809088  0.5427264
```
