# Cross-study pathway heatmap

Creates a heatmap showing enrichment scores or -log10(padj) for top
pathways across all studies.

## Usage

``` r
mx_pathway_heatmap(pathway_results, top_n = 30, value = c("padj", "NES"))
```

## Arguments

- pathway_results:

  A list of per-study pathway enrichment `data.frame` objects.

- top_n:

  Integer. Number of top pathways to display. Default: `30`.

- value:

  Character scalar. Value to display: `"padj"` (default) or `"NES"` (for
  GSEA results).

## Value

A
[`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
object.

## See also

[`mx_pathway_meta`](https://hossainlab.github.io/metaXpress/reference/mx_pathway_meta.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  mx_pathway_heatmap(per_study_pathways, top_n = 20)
} # }
```
