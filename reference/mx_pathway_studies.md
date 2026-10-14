# Cross-study pathway meta-analysis

Runs pre-ranked GSEA separately in every study and combines the
pathway-level evidence across studies. This answers whether a pathway is
consistently perturbed across cohorts, which is more robust than
enrichment of a single pooled gene list when studies differ in platform,
depth or gene coverage (Shen & Tseng 2010).

## Usage

``` r
mx_pathway_studies(
  de_results,
  db = c("Hallmarks", "Reactome", "KEGG", "WikiPathways", "GO_BP"),
  gene_sets = NULL,
  min_size = 15,
  max_size = 500,
  weights = "equal",
  seed = 1L,
  ...
)
```

## Arguments

- de_results:

  A named list of per-study DE `data.frame`s or
  [`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
  objects after
  [`mx_de_all`](https://mdjubayerhossain.com/metaXpress/reference/mx_de_all.md).

- db:

  Gene set collection passed to
  [`mx_gene_sets`](https://mdjubayerhossain.com/metaXpress/reference/mx_gene_sets.md).
  Ignored when `gene_sets` is supplied.

- gene_sets:

  Optional named list of gene sets (e.g. from
  [`mx_gene_sets`](https://mdjubayerhossain.com/metaXpress/reference/mx_gene_sets.md)
  or a GMT file) whose IDs match `meta_table$gene_id`.

- min_size, max_size:

  Gene set size limits after restricting to genes in the analysis.
  Defaults: `15` and `500`.

- weights:

  Per-study weights for the Stouffer combination: `"equal"` (default),
  `"sample_size"` (square root of the per-study sample size from
  `n_case` / `n_control`), or a numeric vector with one weight per
  study. Equal weights are the default because GSEA p-values of a large
  cohort are already far more extreme, and sample-size weights would let
  one cohort dominate.

- seed:

  Integer seed for the fgsea permutation sampling, for reproducibility.
  The global RNG state is restored afterwards.

- ...:

  Further arguments passed to
  [`fgsea::fgsea()`](https://rdrr.io/pkg/fgsea/man/fgsea.html).

## Value

A `data.frame` with one row per gene set: `pathway_id`, `n_studies`,
`mean_NES`, `n_up` / `n_down` (studies with positive / negative NES),
`n_sig` (studies with per-study `padj <= 0.05`), `meta_z`,
`meta_pvalue`, `meta_padj`, followed by one `NES_<study>` column per
study. Per-study fgsea tables are in `attr(x, "per_study")`.

## Details

For each study, genes are ranked by the signed z-score derived from the
per-study test (the `stat` column from
[`mx_de`](https://mdjubayerhossain.com/metaXpress/reference/mx_de.md)
when present). Each pathway's two-sided GSEA p-value is converted to a
signed z-score using the sign of its NES, and the z-scores are combined
by weighted Stouffer's method, so pathways enriched in opposite
directions cancel out. This is the p-value-combination variant of MAPE
(Shen & Tseng 2010) with directional combination.

## References

Shen, K. & Tseng, G.C. (2010) Meta-analysis for pathway enrichment
analysis when combining multiple genomic studies. *Bioinformatics*,
**26**(10), 1316–1323.
[doi:10.1093/bioinformatics/btq148](https://doi.org/10.1093/bioinformatics/btq148)

Korotkevich, G. et al. (2016) Fast gene set enrichment analysis.
*bioRxiv*. [doi:10.1101/060012](https://doi.org/10.1101/060012)

## See also

[`mx_pathway_meta`](https://mdjubayerhossain.com/metaXpress/reference/mx_pathway_meta.md),
[`mx_pathway_heatmap`](https://mdjubayerhossain.com/metaXpress/reference/mx_pathway_heatmap.md)

## Examples

``` r
if (requireNamespace("fgsea", quietly = TRUE)) {
  set.seed(2)
  genes <- paste0("G", 1:1000)
  sets  <- list(UP_SET = genes[1:40], NULL_SET = genes[501:560])
  de <- lapply(1:3, function(s) {
    lfc <- rnorm(1000, c(rep(0.8, 40), rep(0, 960)), 0.3)
    p <- 2 * pnorm(-abs(lfc / 0.3))
    data.frame(gene_id = genes, log2FC = lfc, pvalue = p,
               padj = p.adjust(p, "BH"))
  })
  names(de) <- c("S1", "S2", "S3")
  mx_pathway_studies(de, gene_sets = sets, min_size = 10)
}
#> Running per-study GSEA on 3 studies...
#>   pathway_id n_studies  mean_NES n_up n_down n_sig    meta_z  meta_pvalue
#> 1     UP_SET         3 3.0669911    3      0     3 17.732418 2.356787e-70
#> 2   NULL_SET         3 0.4940344    2      1     0  1.166798 2.432918e-01
#>      meta_padj   NES_S1     NES_S2   NES_S3
#> 1 4.713575e-70 3.037499  3.1581902 3.005284
#> 2 2.432918e-01 1.102117 -0.7690123 1.148999
```
