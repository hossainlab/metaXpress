# Filter studies below a QC threshold

Removes studies whose `qc_score` is below the specified threshold.
Studies with `NA` QC scores (i.e.,
[`mx_qc_study`](https://mdjubayerhossain.com/metaXpress/reference/mx_qc_study.md)
has not been run) are always removed with a warning.

## Usage

``` r
mx_filter_studies(studies, qc_threshold = 7)
```

## Arguments

- studies:

  A named list of
  [`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
  objects.

- qc_threshold:

  Numeric scalar. Minimum acceptable QC score (inclusive). Default: `7`.

## Value

A filtered list of
[`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
objects.

## See also

[`mx_qc_study`](https://mdjubayerhossain.com/metaXpress/reference/mx_qc_study.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  studies <- mx_fetch_geo(c("GSE12345", "GSE67890"))
  studies <- mx_filter_studies(studies, qc_threshold = 7)
} # }
```
