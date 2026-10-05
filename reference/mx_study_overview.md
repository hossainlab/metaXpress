# Study characteristics overview plot

Creates a multi-panel summary of study-level QC metrics including QC
scores, sample sizes, and sequencing depths.

## Usage

``` r
mx_study_overview(studies)
```

## Arguments

- studies:

  A named list of
  [`metaXpressStudy`](https://hossainlab.github.io/metaXpress/reference/metaXpressStudy-class.md)
  objects.

## Value

A `patchwork` object containing multiple `ggplot2` plots.

## Examples

``` r
if (FALSE) { # \dontrun{
  mx_study_overview(studies)
} # }
```
