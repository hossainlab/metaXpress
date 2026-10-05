# UpSet plot of DEG overlap across studies

Shows the intersection structure of significant DEGs across studies
using an UpSet plot.

## Usage

``` r
mx_upset(de_results, padj_threshold = 0.05, lfc_threshold = 1)
```

## Arguments

- de_results:

  A named list of per-study DE `data.frame`s.

- padj_threshold:

  Numeric. Adjusted p-value threshold. Default: `0.05`.

- lfc_threshold:

  Numeric. Log2 fold-change threshold. Default: `1`.

## Value

A
[`ComplexHeatmap::UpSet`](https://rdrr.io/pkg/ComplexHeatmap/man/UpSet.html)
object.

## Examples

``` r
if (FALSE) { # \dontrun{
  mx_upset(de_results)
} # }
```
