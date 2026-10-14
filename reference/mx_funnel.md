# Funnel plot for one gene

Plots each study's log2 fold change against its standard error, with the
pooled estimate and pseudo 95% confidence limits. Asymmetry suggests
small-study effects; see
[`mx_small_study`](https://mdjubayerhossain.com/metaXpress/reference/mx_small_study.md)
for formal tests.

## Usage

``` r
mx_funnel(gene, de_results, meta_result = NULL)
```

## Arguments

- gene:

  Gene identifier.

- de_results:

  Named list of per-study DE `data.frame`s or
  [`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
  objects.

- meta_result:

  Optional
  [`metaXpressResult`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressResult-class.md);
  its pooled estimate is used as the funnel centre. Otherwise the
  inverse-variance fixed-effect estimate is used.

## Value

A `ggplot2` object.

## References

Sterne, J.A.C. & Egger, M. (2001) Funnel plots for detecting bias in
meta-analysis: guidelines on choice of axis. *Journal of Clinical
Epidemiology*, **54**(10), 1046–1055.
[doi:10.1016/S0895-4356(01)00377-8](https://doi.org/10.1016/S0895-4356%2801%2900377-8)

## See also

[`mx_small_study`](https://mdjubayerhossain.com/metaXpress/reference/mx_small_study.md),
[`mx_forest`](https://mdjubayerhossain.com/metaXpress/reference/mx_forest.md)

## Examples

``` r
ad <- readRDS(system.file("extdata", "ad_brain_de.rds",
                          package = "metaXpress"))
mx_funnel("GFAP", ad$de_results)

```
