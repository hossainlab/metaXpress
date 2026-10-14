# Remove redundant pathways

Greedy redundancy filter: pathways are visited from most to least
significant and a pathway is dropped when its gene-set Jaccard
similarity with an already retained pathway is at least
`jaccard_threshold`.

## Usage

``` r
mx_pathway_dedup(pathway_results, gene_sets = NULL, jaccard_threshold = 0.5)
```

## Arguments

- pathway_results:

  A `data.frame` with `pathway_id` and a p-value column (`padj` or
  `meta_padj`).

- gene_sets:

  The named list of gene sets used for the enrichment. When `NULL`,
  results are returned unchanged with a warning.

- jaccard_threshold:

  Numeric. Default: `0.5`.

## Value

The filtered `data.frame`, ordered by significance.

## Examples

``` r
sets <- list(A = letters[1:10], B = letters[1:9], C = letters[15:25])
res  <- data.frame(pathway_id = c("A", "B", "C"), padj = c(0.01, 0.02, 0.03))
mx_pathway_dedup(res, gene_sets = sets)
#>   pathway_id padj
#> 1          A 0.01
#> 3          C 0.03
```
