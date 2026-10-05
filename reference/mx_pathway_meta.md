# Perform pathway meta-analysis on meta-analysis results

Runs over-representation analysis (ORA) or gene set enrichment analysis
(GSEA) using the meta-analysis gene rankings, then combines pathway
results across studies.

## Usage

``` r
mx_pathway_meta(
  meta_result,
  db = c("Hallmarks", "KEGG", "Reactome", "GO_BP"),
  method = c("ORA", "GSEA"),
  padj_threshold = 0.05,
  lfc_threshold = 1,
  ...
)
```

## Arguments

- meta_result:

  A
  [`metaXpressResult`](https://hossainlab.github.io/metaXpress/reference/metaXpressResult-class.md)
  object from
  [`mx_meta`](https://hossainlab.github.io/metaXpress/reference/mx_meta.md).

- db:

  Character scalar. Gene set database. One of `"Hallmarks"` (default),
  `"KEGG"`, `"Reactome"`, or `"GO_BP"`.

- method:

  Character scalar. Enrichment method. One of `"ORA"` (default) or
  `"GSEA"`.

- padj_threshold:

  Numeric. Adjusted p-value threshold for defining the significant gene
  list (ORA only). Default: `0.05`.

- lfc_threshold:

  Numeric. Log2 fold-change threshold (ORA only). Default: `1`.

- ...:

  Additional arguments passed to clusterProfiler.

## Value

A
[`metaXpressResult`](https://hossainlab.github.io/metaXpress/reference/metaXpressResult-class.md)
object with the `pathway_result` slot populated. The `pathway_result`
`data.frame` contains: `pathway_id`, `pathway_name`, `n_genes`,
`pvalue`, `padj`, `gene_ratio`.

## References

Wu, T. et al. (2021) clusterProfiler 4.0: a universal enrichment tool
for interpreting omics data. *The Innovation*, **2**(3), 100141.
[doi:10.1016/j.xinn.2021.100141](https://doi.org/10.1016/j.xinn.2021.100141)

Liberzon, A. et al. (2015) The Molecular Signatures Database (MSigDB)
hallmark gene set collection. *Cell Systems*, **1**(6), 417–425.
[doi:10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004)

Subramanian, A. et al. (2005) Gene set enrichment analysis: a
knowledge-based approach for interpreting genome-wide expression
profiles. *Proceedings of the National Academy of Sciences*,
**102**(43), 15545–15550.
[doi:10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102)

## See also

[`mx_pathway_consensus`](https://hossainlab.github.io/metaXpress/reference/mx_pathway_consensus.md),
[`mx_pathway_dedup`](https://hossainlab.github.io/metaXpress/reference/mx_pathway_dedup.md),
[`mx_pathway_heatmap`](https://hossainlab.github.io/metaXpress/reference/mx_pathway_heatmap.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  result <- mx_pathway_meta(meta_result, db = "Hallmarks")
  head(result@pathway_result[order(result@pathway_result$padj), ])
} # }
```
