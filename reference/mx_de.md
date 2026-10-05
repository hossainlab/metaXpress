# Run differential expression analysis on a single study

Performs differential expression analysis on one `metaXpressStudy` using
DESeq2, edgeR, or limma-voom. Stores the full result (all genes) in the
`de_result` slot of the returned study object.

## Usage

``` r
mx_de(
  study,
  method = c("DESeq2", "edgeR", "limma-voom"),
  formula = ~condition,
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

  A formula specifying the design. Default: `~ condition`. The
  `condition` column in `study@metadata` is used as the contrast
  variable; the last level is the reference.

- ...:

  Additional arguments passed to the underlying DE function.

## Value

The input `study` with the `de_result` slot populated. `de_result` is a
`data.frame` with columns: `gene_id`, `log2FC`, `pvalue`, `padj`,
`baseMean`, `method`.

## Details

**CRITICAL:** Always use `padj` (not raw `pvalue`) for significance
filtering downstream. The default significance threshold is
`padj <= 0.05` AND `|log2FC| >= 1`.

## References

Love, M.I., Huber, W. & Anders, S. (2014) Moderated estimation of fold
change and dispersion for RNA-seq data with DESeq2. *Genome Biology*,
**15**(12), 550.
[doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)

Robinson, M.D., McCarthy, D.J. & Smyth, G.K. (2010) edgeR: a
Bioconductor package for differential expression analysis of digital
gene expression data. *Bioinformatics*, **26**(1), 139–140.
[doi:10.1093/bioinformatics/btp616](https://doi.org/10.1093/bioinformatics/btp616)

Ritchie, M.E. et al. (2015) limma powers differential expression
analyses for RNA-sequencing and microarray studies. *Nucleic Acids
Research*, **43**(7), e47.
[doi:10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007)

Law, C.W. et al. (2014) voom: precision weights unlock linear model
analysis tools for RNA-seq read counts. *Genome Biology*, **15**(2),
R29.
[doi:10.1186/gb-2014-15-2-r29](https://doi.org/10.1186/gb-2014-15-2-r29)

## See also

[`mx_de_all`](https://mdjubayerhossain.com/metaXpress/reference/mx_de_all.md),
[`mx_de_summary`](https://mdjubayerhossain.com/metaXpress/reference/mx_de_summary.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  study <- mx_de(study, method = "DESeq2", formula = ~ condition)
  head(study@de_result)
} # }
```
