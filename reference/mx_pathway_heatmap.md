# Cross-study pathway heatmap

Creates a heatmap showing enrichment scores or -log10(padj) for top
pathways across all studies.

## Usage

``` r
mx_pathway_heatmap(pathway_results, top_n = 30, value = c("padj", "NES"))
```

## Arguments

- pathway_results:

  A list of per-study pathway enrichment `data.frame` objects, or a
  single `data.frame` from
  [`mx_pathway_meta`](https://mdjubayerhossain.com/metaXpress/reference/mx_pathway_meta.md).

- top_n:

  Integer. Number of top pathways to display. Default: `30`.

- value:

  Character scalar. Value to display: `"padj"` (default) or `"NES"` (for
  GSEA results).

## Value

A `ggplot2` object.

## See also

[`mx_pathway_meta`](https://mdjubayerhossain.com/metaXpress/reference/mx_pathway_meta.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  mx_pathway_heatmap(result@pathway_result, top_n = 20)
} # }
```
