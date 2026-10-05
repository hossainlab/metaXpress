# I-squared distribution plot

Displays the distribution of I-squared heterogeneity values across all
genes in the meta-analysis, with reference lines at 25% (low
heterogeneity) and 75% (high heterogeneity). These thresholds follow the
classification proposed by Higgins et al. (2003): I-squared values of
25%, 50%, and 75% correspond to low, moderate, and high heterogeneity,
respectively.

## Usage

``` r
mx_heterogeneity_plot(meta_result)
```

## Arguments

- meta_result:

  A
  [`metaXpressResult`](https://hossainlab.github.io/metaXpress/reference/metaXpressResult-class.md)
  object.

## Value

A `ggplot2` object.

## References

Higgins, J.P.T. et al. (2003) Measuring inconsistency in meta-analyses.
*BMJ*, **327**(7414), 557–560.
[doi:10.1136/bmj.327.7414.557](https://doi.org/10.1136/bmj.327.7414.557)

## See also

[`mx_heterogeneity`](https://hossainlab.github.io/metaXpress/reference/mx_heterogeneity.md),
[`mx_forest`](https://hossainlab.github.io/metaXpress/reference/mx_forest.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  mx_heterogeneity_plot(meta_result)
} # }
```
