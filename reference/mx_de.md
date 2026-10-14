# Run differential expression analysis on a single study

Performs differential expression analysis on one `metaXpressStudy` using
DESeq2, edgeR (quasi-likelihood), or limma-voom. Stores the full result
(all tested genes) in the `de_result` slot of the returned study object,
in a common format that
[`mx_meta`](https://mdjubayerhossain.com/metaXpress/reference/mx_meta.md)
consumes directly.

## Usage

``` r
mx_de(
  study,
  method = c("DESeq2", "edgeR", "limma-voom"),
  formula = ~condition,
  contrast = NULL,
  shrink = c("none", "apeglm", "ashr", "normal"),
  ...
)
```

## Arguments

- study:

  A
  [`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
  object with raw counts.

- method:

  Character scalar. DE method. One of `"DESeq2"` (default), `"edgeR"`,
  or `"limma-voom"`.

- formula:

  A design formula whose variables are columns of `colData(study)`. Must
  contain `condition`; any other terms (e.g. `~ sex + age + condition`)
  are adjusted for. Default: `~ condition`.

- contrast:

  Optional length-2 character vector `c(case_level, control_level)`
  naming the two `condition` levels to compare. Fold changes are
  `case / control`. When `NULL`, the control level is auto-detected
  (e.g. `"control"`, `"normal"`, `"healthy"`, `"untreated"`) and the
  first other level is the case. Set this explicitly for every study in
  a meta-analysis so fold-change directions agree.

- shrink:

  Character scalar. Optional log fold-change shrinkage for DESeq2:
  `"none"` (default), `"apeglm"`, `"ashr"` or `"normal"`. Shrunken
  values are added as `log2FC_shrunk` for ranking and plotting; the
  unshrunken maximum-likelihood `log2FC` and `lfcSE` are always kept
  because inverse-variance meta-analysis requires unbiased estimates
  with their sampling variance.

- ...:

  Additional arguments passed to
  [`DESeq2::DESeq()`](https://rdrr.io/pkg/DESeq2/man/DESeq.html),
  [`edgeR::glmQLFit()`](https://rdrr.io/pkg/edgeR/man/glmQLFit.html) or
  [`limma::lmFit()`](https://rdrr.io/pkg/limma/man/lmFit.html).

## Value

The input `study` with the `de_result` slot populated. `de_result` is a
`data.frame` with columns:

- gene_id:

  Gene identifier (rownames of the count matrix).

- log2FC:

  Log2 fold change, case vs control.

- lfcSE:

  Standard error of `log2FC` (log2 scale).

- stat:

  Signed standard-normal test statistic
  (`sign(log2FC) * qnorm(pvalue / 2, lower.tail = FALSE)` for
  edgeR/limma, the Wald statistic for DESeq2).

- pvalue, padj:

  Raw and Benjamini-Hochberg adjusted p-values.

- baseMean:

  Mean normalized expression.

- n_case, n_control:

  Number of samples in each group, used by
  [`mx_meta`](https://mdjubayerhossain.com/metaXpress/reference/mx_meta.md)
  for sample-size weighting.

- method:

  DE engine used.

The case and control labels are stored in `attr(de_result, "contrast")`.

## Details

For edgeR and limma-voom, low-count genes are removed with
[`edgeR::filterByExpr()`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html)
before testing; DESeq2 applies its own independent filtering to `padj`
but reports all genes. Genes removed by filtering are treated as *not
measured* in that study and are handled explicitly by Module 5
([`mx_missing_summary`](https://mdjubayerhossain.com/metaXpress/reference/mx_missing_summary.md)).

For edgeR, standard errors are back-calculated from the quasi-likelihood
test as `|log2FC| / |z|`, where `z` is the normal quantile of the
two-sided p-value. This places all three engines on the same log2-scale
effect/variance footing required for inverse-variance meta-analysis (Rau
et al. 2014).

Use `padj` (not raw `pvalue`) for per-study significance calls.

## References

Love, M.I., Huber, W. & Anders, S. (2014) Moderated estimation of fold
change and dispersion for RNA-seq data with DESeq2. *Genome Biology*,
**15**(12), 550.
[doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)

Chen, Y., Lun, A.T.L. & Smyth, G.K. (2016) From reads to genes to
pathways: differential expression analysis of RNA-Seq experiments using
Rsubread and the edgeR quasi-likelihood pipeline. *F1000Research*,
**5**, 1438.
[doi:10.12688/f1000research.8987.2](https://doi.org/10.12688/f1000research.8987.2)

Law, C.W. et al. (2014) voom: precision weights unlock linear model
analysis tools for RNA-seq read counts. *Genome Biology*, **15**(2),
R29.
[doi:10.1186/gb-2014-15-2-r29](https://doi.org/10.1186/gb-2014-15-2-r29)

Zhu, A., Ibrahim, J.G. & Love, M.I. (2019) Heavy-tailed prior
distributions for sequence count data: removing the noise and preserving
large differences. *Bioinformatics*, **35**(12), 2084–2092.
[doi:10.1093/bioinformatics/bty895](https://doi.org/10.1093/bioinformatics/bty895)

## See also

[`mx_de_all`](https://mdjubayerhossain.com/metaXpress/reference/mx_de_all.md),
[`mx_de_summary`](https://mdjubayerhossain.com/metaXpress/reference/mx_de_summary.md),
[`mx_meta`](https://mdjubayerhossain.com/metaXpress/reference/mx_meta.md)

## Examples

``` r
if (requireNamespace("airway", quietly = TRUE) &&
    requireNamespace("DESeq2", quietly = TRUE)) {
  # Himes et al. (2014) airway smooth muscle cells, GEO GSE52778:
  # dexamethasone vs untreated, paired by donor cell line.
  data("airway", package = "airway")
  study <- mx_from_se(airway, condition_col = "dex",
                      accession = "GSE52778")
  study <- study[rowSums(SummarizedExperiment::assay(study)) >= 10, ]
  study <- mx_de(study, formula = ~ cell + condition,
                 contrast = c("trt", "untrt"))
  de <- study@de_result
  head(de[order(de$padj), ])
}
#>   Running DESeq2 on study: GSE52778 (trt vs untrt)
#> converting counts to integer mode
#>               gene_id   log2FC     lfcSE     stat        pvalue          padj
#> 8109  ENSG00000152583 4.574967 0.1842414 24.83137 4.110667e-136 7.412355e-132
#> 9865  ENSG00000165995 3.291099 0.1330527 24.73530 4.463384e-135 4.024187e-131
#> 4359  ENSG00000120129 2.947850 0.1218764 24.18720 3.033839e-129 1.823540e-125
#> 2214  ENSG00000101347 3.767022 0.1559919 24.14882 7.682657e-129 3.463342e-125
#> 13405 ENSG00000189221 3.353655 0.1421810 23.58722 5.212705e-123 1.879910e-119
#> 14808 ENSG00000211445 3.730439 0.1663984 22.41872 2.585059e-111 7.768964e-108
#>         baseMean n_case n_control method
#> 8109    997.4447      4         4 DESeq2
#> 9865    495.0957      4         4 DESeq2
#> 4359   3409.0384      4         4 DESeq2
#> 2214  12703.4128      4         4 DESeq2
#> 13405  2341.7807      4         4 DESeq2
#> 14808 12285.7001      4         4 DESeq2
```
