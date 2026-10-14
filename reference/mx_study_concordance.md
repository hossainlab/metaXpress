# Check fold-change concordance between studies

Computes pairwise correlations of per-study log2 fold changes over genes
measured in both studies. A negative correlation usually means the
case/control orientation is flipped in one study, which would invalidate
any directional meta-analysis; low correlation flags biologically or
technically divergent studies.

## Usage

``` r
mx_study_concordance(
  de_results,
  genes = c("significant", "all"),
  padj_threshold = 0.05,
  cor_method = c("spearman", "pearson")
)
```

## Arguments

- de_results:

  A named list of `data.frame` objects (the `de_result` slots produced
  by
  [`mx_de`](https://mdjubayerhossain.com/metaXpress/reference/mx_de.md)),
  or a named list of
  [`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
  objects after
  [`mx_de_all`](https://mdjubayerhossain.com/metaXpress/reference/mx_de_all.md).
  Required columns: `gene_id`, `log2FC`, `pvalue`, `padj`; effect-size
  methods also need `lfcSE`.

- genes:

  Which genes to correlate: `"significant"` (default; genes with
  `padj <= padj_threshold` in at least one of the two studies) or
  `"all"`.

- padj_threshold:

  Numeric. Default: `0.05`.

- cor_method:

  Correlation method. Default: `"spearman"`.

## Value

A `data.frame` with one row per study pair: `study_1`, `study_2`,
`n_genes`, `correlation`, `sign_agreement` (fraction of genes with the
same log2FC sign) and a `flag`.

## See also

[`mx_meta`](https://mdjubayerhossain.com/metaXpress/reference/mx_meta.md),
[`mx_de`](https://mdjubayerhossain.com/metaXpress/reference/mx_de.md)

## Examples

``` r
set.seed(1)
base <- rnorm(300)
de <- lapply(c(1, 1, -1), function(dir) {
  lfc <- dir * base + rnorm(300, 0, 0.3); p <- 2 * pnorm(-abs(lfc / 0.3))
  data.frame(gene_id = paste0("G", 1:300), log2FC = lfc, pvalue = p,
             padj = p.adjust(p, "BH"))
})
names(de) <- c("A", "B", "C_flipped")
mx_study_concordance(de)
#> Warning: Negative fold-change correlation between studies: A vs C_flipped; B vs C_flipped. Check that every study uses the same case/control contrast.
#>   study_1   study_2 n_genes correlation sign_agreement
#> 1       A         B     177   0.8957839     0.98305085
#> 2       A C_flipped     177  -0.8905901     0.01694915
#> 3       B C_flipped     176  -0.9020361     0.01704545
#>                                   flag
#> 1                                   ok
#> 2 NEGATIVE: check contrast orientation
#> 3 NEGATIVE: check contrast orientation
```
