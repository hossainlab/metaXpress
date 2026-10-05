# Reannotate gene IDs to a common namespace

Converts gene identifiers in one or more studies to a unified namespace
(Ensembl gene ID, gene symbol, or Entrez ID) using AnnotationDbi. Genes
that cannot be mapped are dropped with a warning. Harmonising gene
identifiers across studies is a prerequisite for any cross-study
comparison (Stokholm et al. 2024).

## Usage

``` r
mx_reannotate(
  studies,
  org = "Homo sapiens",
  target_id = c("SYMBOL", "ENSEMBL", "ENTREZID"),
  BPPARAM = BiocParallel::bpparam()
)
```

## Arguments

- studies:

  A named list of
  [`metaXpressStudy`](https://hossainlab.github.io/metaXpress/reference/metaXpressStudy-class.md)
  objects.

- org:

  Character scalar. Target organism. Currently only `"Homo sapiens"` is
  supported. Default: `"Homo sapiens"`.

- target_id:

  Character scalar. Target gene identifier type. One of `"SYMBOL"`
  (default), `"ENSEMBL"`, or `"ENTREZID"`.

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object controlling parallelization. Default:
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html).
  Set to
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
  for sequential execution.

## Value

The input list of studies with `counts` rownames translated to the
target namespace.

## References

Stokholm, A., Rabaglino, M.B. & Kadarmideen, H.N. (2024) A protocol for
cross-platform meta-analysis of microarray and RNA-seq data. *Current
Protocols*, **4**(12), e70046.
[doi:10.1002/cpz1.70046](https://doi.org/10.1002/cpz1.70046)

Pagès, H. et al. (2024) AnnotationDbi: manipulation of SQLite-based
annotations in Bioconductor. R package version 1.66.0.
<https://bioconductor.org/packages/AnnotationDbi>

## See also

[`mx_normalize`](https://hossainlab.github.io/metaXpress/reference/mx_normalize.md),
[`mx_align_genes`](https://hossainlab.github.io/metaXpress/reference/mx_align_genes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  studies <- mx_reannotate(studies, org = "Homo sapiens",
                            target_id = "SYMBOL")
} # }
```
