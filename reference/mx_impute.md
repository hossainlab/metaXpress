# Impute missing gene statistics across studies

Handles genes not measured in all studies by applying one of four
imputation strategies before meta-analysis.

## Usage

``` r
mx_impute(de_results, method = c("exclude", "mean", "knn", "weighted"))
```

## Arguments

- de_results:

  A named list of `data.frame` objects from
  [`mx_de_all`](https://mdjubayerhossain.com/metaXpress/reference/mx_de_all.md).

- method:

  Character scalar. Imputation strategy. One of:

  `"exclude"`

  :   Exclude genes not present in all studies (most conservative;
      default).

  `"mean"`

  :   Impute missing log2FC with the mean across available studies;
      impute p-values as 1.

  `"knn"`

  :   K-nearest-neighbour imputation of log2FC across genes with similar
      cross-study fold-change profiles.

  `"weighted"`

  :   Mean log2FC weighted by gene coverage.

## Value

The input list with missing values filled according to the chosen
strategy.

## Caution

The `"mean"`, `"knn"` and `"weighted"` strategies fill a study's missing
fold change with values derived from the *other* studies. This is
circular: it inflates the apparent number of studies per gene, shrinks
heterogeneity (I-squared) and has no standard error, so imputed rows are
ignored by effect-size methods. For inference, prefer
[`mx_meta`](https://mdjubayerhossain.com/metaXpress/reference/mx_meta.md)
with `min_studies`, which pools each gene over the studies that measured
it, or
[`mx_filter_coverage`](https://mdjubayerhossain.com/metaXpress/reference/mx_filter_coverage.md).
Imputation is mainly useful for visualising complete matrices.

## References

Villatoro-García, J.A. et al. (2022) DExMA: an R package for performing
gene expression meta-analysis with missing genes. *Mathematics*,
**10**(18), 3376.
[doi:10.3390/math10183376](https://doi.org/10.3390/math10183376)

## See also

[`mx_missing_summary`](https://mdjubayerhossain.com/metaXpress/reference/mx_missing_summary.md),
[`mx_filter_coverage`](https://mdjubayerhossain.com/metaXpress/reference/mx_filter_coverage.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  de_results <- mx_impute(de_results, method = "mean")
} # }
```
