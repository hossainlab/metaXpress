# Meta-analyse cell-type composition shifts across studies

Scores each sample for each cell type as the mean within-study z-score
of log-CPM expression across that cell type's marker genes (an
MCP-counter style abundance proxy), computes Hedges' g for case vs
control in every study, and combines the per-study effects with a
DerSimonian-Laird random-effects model. A significant shift indicates
that part of the bulk differential expression may reflect changes in
tissue composition.

## Usage

``` r
mx_cell_composition(studies, markers, min_markers = 3)
```

## Arguments

- studies:

  A named list of
  [`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
  objects.

- markers:

  A named list mapping cell types to marker gene identifiers, or a
  `data.frame` with columns `gene` and `cell_type`. Identifiers must
  match `rownames` of the studies.

- min_markers:

  Minimum number of expressed markers a cell type needs in a study to be
  scored there. Default 3.

## Value

A `data.frame` sorted by p-value with columns `cell_type`, `meta_g`
(pooled Hedges' g; positive = higher in cases), `pvalue`, `padj`,
`i_squared`, `n_studies`, and one `g_<study>` column per study.
Per-sample scores are stored in `attr(x, "scores")`.

## References

Becht, E. et al. (2016) Estimating the population abundance of
tissue-infiltrating immune and stromal cell populations using gene
expression. *Genome Biology*, **17**, 218.
[doi:10.1186/s13059-016-1070-5](https://doi.org/10.1186/s13059-016-1070-5)

Hedges, L. V. (1981) Distribution theory for Glass's estimator of effect
size and related estimators. *Journal of Educational Statistics*,
**6**(2), 107–128.
[doi:10.3102/10769986006002107](https://doi.org/10.3102/10769986006002107)

## See also

[`mx_coexpression_modules`](https://mdjubayerhossain.com/metaXpress/reference/mx_coexpression_modules.md),
[`mx_downstream_plot`](https://mdjubayerhossain.com/metaXpress/reference/mx_downstream_plot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  markers <- list(T_cells = c("CD3D", "CD3E", "CD2"),
                  Fibroblasts = c("COL1A1", "COL3A1", "DCN"))
  comp <- mx_cell_composition(studies, markers)
} # }
```
