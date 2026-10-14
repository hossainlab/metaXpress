# Remove studies that fail quality control

Keeps studies according to the per-criterion QC results from
[`mx_qc_study`](https://mdjubayerhossain.com/metaXpress/reference/mx_qc_study.md).
Failures, not the summed score, drive the decision, because criteria
that cannot be evaluated (common for GEO data) should not count against
a study.

## Usage

``` r
mx_filter_studies(
  studies,
  required = c(1, 7, 10),
  max_failed = Inf,
  qc_threshold = NULL
)
```

## Arguments

- studies:

  A named list of
  [`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
  objects after
  [`mx_qc_study`](https://mdjubayerhossain.com/metaXpress/reference/mx_qc_study.md).

- required:

  Integer vector of criteria that must not fail. Default: `c(1, 7, 10)`
  (\>= 3 samples per group, complete condition labels, raw counts) – the
  minimum for valid count-based DE.

- max_failed:

  Integer. Maximum number of other failed criteria allowed. Default:
  `Inf` (only `required` criteria exclude).

- qc_threshold:

  Optional minimum `qc_score` (criteria passed). Default `NULL`: not
  used.

## Value

The filtered list. A message lists removed studies and reasons.

## Examples

``` r
if (requireNamespace("airway", quietly = TRUE)) {
  data("airway", package = "airway")
  s <- mx_qc_study(mx_from_se(airway, condition_col = "dex",
                              accession = "GSE52778"))
  length(mx_filter_studies(list(airway = s)))
}
#> [1] 1
```
