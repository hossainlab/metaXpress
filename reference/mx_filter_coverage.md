# Filter genes by minimum study coverage

Retains only genes present in at least `min_studies` studies. This is
the simplest approach to handling missing genes and is applied
automatically inside
[`mx_meta`](https://hossainlab.github.io/metaXpress/reference/mx_meta.md)
via the `min_studies` argument.

## Usage

``` r
mx_filter_coverage(de_results, min_studies = 2)
```

## Arguments

- de_results:

  A named list of `data.frame` objects from
  [`mx_de_all`](https://hossainlab.github.io/metaXpress/reference/mx_de_all.md).

- min_studies:

  Integer. Minimum number of studies a gene must appear in. Default:
  `2`.

## Value

The input list with each `data.frame` restricted to genes meeting the
coverage threshold.

## See also

[`mx_missing_summary`](https://hossainlab.github.io/metaXpress/reference/mx_missing_summary.md),
[`mx_impute`](https://hossainlab.github.io/metaXpress/reference/mx_impute.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  de_results <- mx_filter_coverage(de_results, min_studies = 3)
} # }
```
