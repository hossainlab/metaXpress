# Find consensus co-expression modules across studies

Builds a gene-gene Pearson correlation matrix on log-CPM expression in
each study, combines them into a consensus (by default the element-wise
minimum, as in WGCNA consensus networks, so a link must hold in every
study), clusters genes hierarchically, and cuts the tree into modules.
Each module is summarised by its eigengene (first principal component)
in every study; eigengene differences between cases and controls are
pooled across studies with a random-effects model.

## Usage

``` r
mx_coexpression_modules(
  studies,
  genes = NULL,
  n_genes = 2000,
  consensus = c("min", "mean"),
  cut_height = 0.25,
  min_module_size = 20,
  meta_result = NULL
)
```

## Arguments

- studies:

  A named list of
  [`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
  objects.

- genes:

  Optional character vector of genes to cluster (e.g. meta-DEGs).
  Default `NULL` uses the `n_genes` most variable genes shared by all
  studies.

- n_genes:

  Number of most variable genes to use when `genes` is `NULL`. Default
  2000.

- consensus:

  How to combine study correlation matrices: `"min"` (conservative,
  default) or `"mean"`.

- cut_height:

  Tree cut height on the distance \\(1 - r)/2\\. The default 0.25 groups
  genes whose average consensus correlation is at least about 0.5.

- min_module_size:

  Modules smaller than this are merged into module `"M0"` (unassigned).
  Default 20.

- meta_result:

  Optional
  [`metaXpressResult`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressResult-class.md).
  If given, each module reports the fraction of its genes that are
  significant meta-DEGs (`meta_padj <= 0.05`) up or down.

## Value

A list with

- modules:

  `data.frame` of `gene_id`, `module`, `kME` (mean correlation with the
  module eigengene across studies).

- module_stats:

  `data.frame` with one row per module: `module`, `n_genes`, `meta_g`,
  `pvalue`, `padj`, `i_squared`, per-study `g_<study>`, `hub_genes`, and
  optionally `pct_meta_up`, `pct_meta_down`.

- eigengenes:

  List of per-study sample x module eigengene matrices.

## References

Langfelder, P. and Horvath, S. (2007) Eigengene networks for studying
the relationships between co-expression modules. *BMC Systems Biology*,
**1**, 54.
[doi:10.1186/1752-0509-1-54](https://doi.org/10.1186/1752-0509-1-54)

Langfelder, P. and Horvath, S. (2008) WGCNA: an R package for weighted
correlation network analysis. *BMC Bioinformatics*, **9**, 559.
[doi:10.1186/1471-2105-9-559](https://doi.org/10.1186/1471-2105-9-559)

## See also

[`mx_cell_composition`](https://mdjubayerhossain.com/metaXpress/reference/mx_cell_composition.md),
[`mx_downstream_plot`](https://mdjubayerhossain.com/metaXpress/reference/mx_downstream_plot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  net <- mx_coexpression_modules(studies, n_genes = 1500,
                                 meta_result = meta_result)
  net$module_stats
} # }
```
