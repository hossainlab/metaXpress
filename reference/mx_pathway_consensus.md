# Find consensus pathways across studies (vote counting)

Counts in how many per-study enrichment analyses each pathway is
significant. Vote counting ignores effect size and is less powerful than
[`mx_pathway_studies`](https://mdjubayerhossain.com/metaXpress/reference/mx_pathway_studies.md);
use it as a descriptive summary.

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

  A list of per-study enrichment `data.frame`s with columns `pathway_id`
  and `padj`, e.g. `attr(mx_pathway_studies(...), "per_study")`.

- min_fraction:

  Numeric. Minimum fraction of studies in which a pathway must be
  significant. Default: `0.5`.

- padj_threshold:

  Numeric. Significance threshold per study. Default: `0.05`.

## Value

A `data.frame` with `pathway_id`, `n_significant_studies` and
`fraction_studies`.

## See also

[`mx_pathway_studies`](https://mdjubayerhossain.com/metaXpress/reference/mx_pathway_studies.md)

## Examples

``` r
pr1 <- data.frame(pathway_id = c("A", "B"), padj = c(0.01, 0.2))
pr2 <- data.frame(pathway_id = c("A", "B"), padj = c(0.03, 0.01))
mx_pathway_consensus(list(pr1, pr2), min_fraction = 1)
#>   pathway_id n_significant_studies fraction_studies
#> 1          A                     2                1
```
