# Impute missing gene statistics across studies

Handles genes not measured in all studies by applying one of four
imputation strategies before meta-analysis.

## Usage

``` r
mx_impute(
  de_results,
  method = c("exclude", "mean", "knn", "weighted"),
  weights = NULL
)
```

## Arguments

- de_results:

  A named list of `data.frame` objects from
  [`mx_de_all`](https://hossainlab.github.io/metaXpress/reference/mx_de_all.md).

- method:

  Character scalar. Imputation strategy. One of:

  `"exclude"`

  :   Exclude genes not present in all studies (most conservative;
      default).

  `"mean"`

  :   Impute missing log2FC with the mean across available studies;
      impute p-values as 1.

  `"knn"`

  :   K-nearest neighbours imputation based on gene expression profiles
      (Hastie et al. 1999).

  `"weighted"`

  :   Weighted imputation proportional to study sample size.

- weights:

  Numeric vector. Study weights for `"weighted"` imputation. Default:
  `NULL` (uses equal weights, equivalent to `"mean"`).

## Value

The input list with missing values filled according to the chosen
strategy.

## References

Hastie, T. et al. (1999) Imputing missing data for gene expression
arrays. Stanford University Statistics Department Technical report.

Villatoro-García, J.A. et al. (2022) Missing gene expression data
imputation for gene-study meta-analysis. *Mathematics*, **10**(18),
3376. [doi:10.3390/math10183376](https://doi.org/10.3390/math10183376)

## See also

[`mx_missing_summary`](https://hossainlab.github.io/metaXpress/reference/mx_missing_summary.md),
[`mx_filter_coverage`](https://hossainlab.github.io/metaXpress/reference/mx_filter_coverage.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  de_results <- mx_impute(de_results, method = "mean")
} # }
```
