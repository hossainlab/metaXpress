# Forest plot for a single gene

Displays per-study log2 fold-change estimates with 95% confidence
intervals and a pooled estimate diamond. I-squared heterogeneity is
annotated. Forest plots are the standard visualization for meta-analysis
results, showing both individual study contributions and the pooled
summary (Lewis & Clarke 2001).

## Usage

``` r
mx_forest(gene, de_results, meta_result = NULL)
```

## Arguments

- gene:

  Character scalar. Gene identifier to plot.

- de_results:

  A named list of per-study DE `data.frame`s (from
  [`mx_de_all`](https://mdjubayerhossain.com/metaXpress/reference/mx_de_all.md)
  result slots), or a list of
  [`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
  objects.

- meta_result:

  Optional. A
  [`metaXpressResult`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressResult-class.md)
  object. When provided, the pooled estimate diamond and I-squared are
  taken from the meta-analysis result.

## Value

A `ggplot2` object.

## References

Lewis, S. & Clarke, M. (2001) Forest plots: trying to see the wood and
the trees. *BMJ*, **322**(7300), 1479–1480.
[doi:10.1136/bmj.322.7300.1479](https://doi.org/10.1136/bmj.322.7300.1479)

## See also

[`mx_volcano`](https://mdjubayerhossain.com/metaXpress/reference/mx_volcano.md),
[`mx_heterogeneity_plot`](https://mdjubayerhossain.com/metaXpress/reference/mx_heterogeneity_plot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  p <- mx_forest("TP53", de_results, meta_result)
  print(p)
} # }
```
