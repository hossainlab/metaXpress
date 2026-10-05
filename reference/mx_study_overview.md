# Study characteristics overview plot

Creates a multi-panel summary of study-level QC metrics including QC
scores, sample sizes, sequencing depths, and gene detection rates.

## Usage

``` r
mx_study_overview(studies)
```

## Arguments

- studies:

  A named list of
  [`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
  objects.

## Value

A `ggplot2` object (or list of plots).

## Examples

``` r
if (FALSE) { # \dontrun{
  mx_study_overview(studies)
} # }
```
