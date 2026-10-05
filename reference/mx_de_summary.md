# Summarise per-study DE results

Returns a summary table showing the number of significant DEGs in each
study at the default threshold (`padj <= 0.05`, `|log2FC| >= 1`).

## Usage

``` r
mx_de_summary(studies, padj_threshold = 0.05, lfc_threshold = 1)
```

## Arguments

- studies:

  A named list of
  [`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
  objects after running
  [`mx_de_all`](https://mdjubayerhossain.com/metaXpress/reference/mx_de_all.md).

- padj_threshold:

  Numeric. Adjusted p-value threshold. Default: 0.05.

- lfc_threshold:

  Numeric. Absolute log2 fold-change threshold. Default: 1.

## Value

A `data.frame` with columns `study`, `n_total`, `n_up`, `n_down`,
`method`.

## See also

[`mx_de_all`](https://mdjubayerhossain.com/metaXpress/reference/mx_de_all.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  summary_tbl <- mx_de_summary(studies)
} # }
```
