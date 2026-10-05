# Compute per-gene heterogeneity statistics

Returns Cochran's Q, degrees of freedom, I-squared, tau-squared, and the
p-value for the Q test for each gene across studies.

## Usage

``` r
mx_heterogeneity(de_results)
```

## Arguments

- de_results:

  A named list of `data.frame` objects from
  [`mx_de_all`](https://hossainlab.github.io/metaXpress/reference/mx_de_all.md),
  or a list of
  [`metaXpressStudy`](https://hossainlab.github.io/metaXpress/reference/metaXpressStudy-class.md)
  objects.

## Value

A `data.frame` with columns `gene_id`, `Q`, `df`, `I_sq`, `tau_sq`,
`p_heterogeneity`.

## References

Cochran, W.G. (1954) The combination of estimates from different
experiments. *Biometrics*, **10**(1), 101–129.
[doi:10.2307/3001666](https://doi.org/10.2307/3001666)

Higgins, J.P.T. & Thompson, S.G. (2002) Quantifying heterogeneity in a
meta-analysis. *Statistics in Medicine*, **21**(11), 1539–1558.
[doi:10.1002/sim.1186](https://doi.org/10.1002/sim.1186)

## See also

[`mx_meta`](https://hossainlab.github.io/metaXpress/reference/mx_meta.md),
[`mx_heterogeneity_plot`](https://hossainlab.github.io/metaXpress/reference/mx_heterogeneity_plot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  het <- mx_heterogeneity(de_results)
  hist(het$I_sq, main = "I-squared distribution")
} # }
```
