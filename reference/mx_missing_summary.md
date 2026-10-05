# Summarise gene coverage across studies

Returns a matrix showing which genes are present (1) or absent (0) in
each study's DE result, enabling informed decisions about how to handle
partial overlap before meta-analysis.

## Usage

``` r
mx_missing_summary(de_results)
```

## Arguments

- de_results:

  A named list of `data.frame` objects from
  [`mx_de_all`](https://hossainlab.github.io/metaXpress/reference/mx_de_all.md),
  or a list of
  [`metaXpressStudy`](https://hossainlab.github.io/metaXpress/reference/metaXpressStudy-class.md)
  objects.

## Value

A binary matrix (genes x studies) where 1 indicates the gene is present
and 0 indicates it is absent. Attribute `"coverage_pct"` contains the
percentage of studies covering each gene.

## See also

[`mx_impute`](https://hossainlab.github.io/metaXpress/reference/mx_impute.md),
[`mx_filter_coverage`](https://hossainlab.github.io/metaXpress/reference/mx_filter_coverage.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  cov_mat <- mx_missing_summary(de_results)
  hist(attr(cov_mat, "coverage_pct"), main = "Gene coverage across studies")
} # }
```
