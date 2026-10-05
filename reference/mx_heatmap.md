# Heatmap of top DEGs across studies

Displays log2 fold-change values for the top differentially expressed
genes across all studies as a heatmap.

## Usage

``` r
mx_heatmap(meta_result, studies, top_n = 50)
```

## Arguments

- meta_result:

  A
  [`metaXpressResult`](https://hossainlab.github.io/metaXpress/reference/metaXpressResult-class.md)
  object.

- studies:

  A named list of
  [`metaXpressStudy`](https://hossainlab.github.io/metaXpress/reference/metaXpressStudy-class.md)
  objects (used to extract per-study log2FC if needed).

- top_n:

  Integer. Number of top genes (by `meta_padj`) to display. Default:
  `50`.

## Value

A
[`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
object.

## See also

[`mx_volcano`](https://hossainlab.github.io/metaXpress/reference/mx_volcano.md),
[`mx_upset`](https://hossainlab.github.io/metaXpress/reference/mx_upset.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  mx_heatmap(meta_result, studies, top_n = 30)
} # }
```
