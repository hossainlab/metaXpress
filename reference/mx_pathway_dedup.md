# Remove redundant pathways

Reduces pathway result sets by removing pathways with high gene set
overlap, retaining the most significant representative from each cluster
of similar pathways.

## Usage

``` r
mx_pathway_dedup(pathway_results, jaccard_threshold = 0.5)
```

## Arguments

- pathway_results:

  A `data.frame` from
  [`mx_pathway_meta`](https://mdjubayerhossain.com/metaXpress/reference/mx_pathway_meta.md)
  or a list thereof.

- jaccard_threshold:

  Numeric. Jaccard similarity threshold above which two pathways are
  considered redundant. Default: `0.5`.

## Value

A `data.frame` with redundant pathways removed.

## References

Gu, Z. & Huebschmann, D. (2023) simplifyEnrichment: a Bioconductor
package for clustering and visualizing functional enrichment results.
*Genomics, Proteomics & Bioinformatics*, **21**(1), 190–202.
[doi:10.1016/j.gpb.2022.04.008](https://doi.org/10.1016/j.gpb.2022.04.008)

## Examples

``` r
if (FALSE) { # \dontrun{
  result <- mx_pathway_dedup(result@pathway_result)
} # }
```
