# Restrict all studies to a common set of genes

Finds the intersection of gene IDs across all studies and returns
studies whose count matrices are restricted to that common set. Gene IDs
are taken from the rownames of each study's count matrix.

## Usage

``` r
mx_align_genes(studies)
```

## Arguments

- studies:

  A named list of
  [`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
  objects.

## Value

The input list of studies, each with `counts` restricted to the
intersection of all gene sets.

## See also

[`mx_reannotate`](https://mdjubayerhossain.com/metaXpress/reference/mx_reannotate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  studies <- mx_align_genes(studies)
} # }
```
