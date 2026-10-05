# Run differential expression on all studies

Applies
[`mx_de`](https://hossainlab.github.io/metaXpress/reference/mx_de.md) to
each study in a list, optionally in parallel using BiocParallel.

## Usage

``` r
mx_de_all(
  studies,
  method = "DESeq2",
  formula = ~condition,
  BPPARAM = BiocParallel::SerialParam(),
  ...
)
```

## Arguments

- studies:

  A named list of
  [`metaXpressStudy`](https://hossainlab.github.io/metaXpress/reference/metaXpressStudy-class.md)
  objects.

- method:

  Character scalar. DE method passed to
  [`mx_de`](https://hossainlab.github.io/metaXpress/reference/mx_de.md).
  Default: `"DESeq2"`.

- formula:

  A formula specifying the design. Default: `~ condition`.

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object controlling parallelization. Defaults to
  [`SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
  (single-core). Use
  [`MulticoreParam()`](https://rdrr.io/pkg/BiocParallel/man/MulticoreParam-class.html)
  for multi-core execution.

- ...:

  Additional arguments passed to
  [`mx_de`](https://hossainlab.github.io/metaXpress/reference/mx_de.md).

## Value

The input list of studies, each with `de_result` slot filled.

## See also

[`mx_de`](https://hossainlab.github.io/metaXpress/reference/mx_de.md),
[`mx_de_summary`](https://hossainlab.github.io/metaXpress/reference/mx_de_summary.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  studies <- mx_de_all(studies, method = "DESeq2")
  # Parallel:
  studies <- mx_de_all(studies, method = "DESeq2",
                        BPPARAM = BiocParallel::MulticoreParam(4))
} # }
```
