# Pathway enrichment of meta-analysis results

Runs pre-ranked gene set enrichment analysis (GSEA, via fgsea) or
over-representation analysis (ORA, hypergeometric test) on the
gene-level meta-analysis.

## Usage

``` r
mx_pathway_meta(
  meta_result,
  db = c("Hallmarks", "Reactome", "KEGG", "WikiPathways", "GO_BP"),
  method = c("GSEA", "ORA"),
  gene_sets = NULL,
  padj_threshold = 0.05,
  lfc_threshold = 0,
  min_size = 15,
  max_size = 500,
  seed = 1L,
  ...
)
```

## Arguments

- meta_result:

  A
  [`metaXpressResult`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressResult-class.md)
  from
  [`mx_meta`](https://mdjubayerhossain.com/metaXpress/reference/mx_meta.md).

- db:

  Gene set collection passed to
  [`mx_gene_sets`](https://mdjubayerhossain.com/metaXpress/reference/mx_gene_sets.md).
  Ignored when `gene_sets` is supplied.

- method:

  `"GSEA"` (default) or `"ORA"`.

- gene_sets:

  Optional named list of gene sets (e.g. from
  [`mx_gene_sets`](https://mdjubayerhossain.com/metaXpress/reference/mx_gene_sets.md)
  or a GMT file) whose IDs match `meta_table$gene_id`.

- padj_threshold, lfc_threshold:

  ORA only: thresholds on `meta_padj` and `|meta_log2FC|` defining the
  gene lists. Defaults: `0.05` and `0`.

- min_size, max_size:

  Gene set size limits after restricting to genes in the analysis.
  Defaults: `15` and `500`.

- seed:

  Integer seed for the fgsea permutation sampling, for reproducibility.
  The global RNG state is restored afterwards.

- ...:

  Further arguments passed to
  [`fgsea::fgsea()`](https://rdrr.io/pkg/fgsea/man/fgsea.html).

## Value

The `meta_result` with `pathway_result` filled:

- GSEA:

  `pathway_id`, `set_size`, `NES`, `pvalue`, `padj`, `leading_edge`.

- ORA:

  `pathway_id`, `direction` (up/down), `set_size`, `n_overlap`,
  `n_list`, `fold_enrichment`, `pvalue`, `padj`, `overlap_genes`.

## Details

GSEA ranks genes by the signed meta-analysis z-score
\\\mathrm{sign}(\beta)\\\Phi^{-1}(1 - p/2)\\, which uses both effect
size and its precision; this avoids arbitrary significance cut-offs and
is the recommended default. ORA tests up- and down-regulated genes
separately against the universe of all genes in the meta-analysis (not
the whole genome), as recommended to avoid inflated enrichment.

## References

Subramanian, A. et al. (2005) Gene set enrichment analysis: a
knowledge-based approach for interpreting genome-wide expression
profiles. *PNAS*, **102**(43), 15545–15550.
[doi:10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102)

Korotkevich, G. et al. (2016) Fast gene set enrichment analysis.
*bioRxiv*. [doi:10.1101/060012](https://doi.org/10.1101/060012)

Timmons, J.A., Szkop, K.J. & Gallagher, I.J. (2015) Multiple sources of
bias confound functional enrichment analysis of global -omics data.
*Genome Biology*, **16**, 186.
[doi:10.1186/s13059-015-0761-7](https://doi.org/10.1186/s13059-015-0761-7)

## See also

[`mx_pathway_studies`](https://mdjubayerhossain.com/metaXpress/reference/mx_pathway_studies.md),
[`mx_gene_sets`](https://mdjubayerhossain.com/metaXpress/reference/mx_gene_sets.md),
[`mx_pathway_dedup`](https://mdjubayerhossain.com/metaXpress/reference/mx_pathway_dedup.md),
[`mx_pathway_heatmap`](https://mdjubayerhossain.com/metaXpress/reference/mx_pathway_heatmap.md)

## Examples

``` r
if (requireNamespace("fgsea", quietly = TRUE)) {
  set.seed(1)
  genes <- paste0("G", 1:1000)
  sets  <- list(UP_SET = genes[1:40], NULL_SET = genes[501:560])
  de <- lapply(1:3, function(s) {
    lfc <- rnorm(1000, c(rep(1, 40), rep(0, 960)), 0.3)
    p <- 2 * pnorm(-abs(lfc / 0.3))
    data.frame(gene_id = genes, log2FC = lfc, lfcSE = 0.3, pvalue = p,
               padj = p.adjust(p, "BH"))
  })
  res <- mx_meta(de)
  res <- mx_pathway_meta(res, gene_sets = sets, min_size = 10)
  res@pathway_result
}
#> Running meta-analysis (random_effects) on 3 studies...
#>   pathway_id set_size       NES       pvalue         padj
#> 2     UP_SET       40  2.985097 1.000000e-50 2.000000e-50
#> 1   NULL_SET       60 -1.050589 3.446328e-01 3.446328e-01
#>                                                                                                                                             leading_edge
#> 2 G20;G4;G15;G25;G40;G30;G7;G36;G39;G16;G22;G10;G5;G23;G19;G34;G11;G33;G12;G32;G17;G1;G31;G37;G27;G26;G13;G18;G29;G35;G3;G6;G21;G38;G24;G2;G14;G28;G8;G9
#> 1                                                                        G511;G526;G513;G521;G559;G517;G502;G549;G510;G519;G533;G536;G507;G512;G546;G531
```
