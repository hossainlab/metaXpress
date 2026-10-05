# Find consensus pathways across studies

Identifies pathways that are significantly enriched in the majority of
individual studies, indicating robust cross-study signal.

## Usage

``` r
mx_pathway_consensus(
  pathway_results,
  min_fraction = 0.5,
  padj_threshold = 0.05
)
```

## Arguments

- pathway_results:

  A list of `data.frame` objects from per-study enrichment analyses.

- min_fraction:

  Numeric. Minimum fraction of studies in which a pathway must be
  significant. Default: `0.5`.

- padj_threshold:

  Numeric. Significance threshold per study. Default: `0.05`.

## Value

A `data.frame` of consensus pathways with columns `pathway_id`,
`pathway_name`, `n_significant_studies`, `fraction_studies`.

## See also

[`mx_pathway_meta`](https://hossainlab.github.io/metaXpress/reference/mx_pathway_meta.md),
[`mx_pathway_dedup`](https://hossainlab.github.io/metaXpress/reference/mx_pathway_dedup.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  consensus <- mx_pathway_consensus(pathway_results)
} # }
```
