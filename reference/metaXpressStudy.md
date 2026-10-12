# Construct a metaXpressStudy

Builds a
[`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
from a count matrix and a sample metadata table.

## Usage

``` r
metaXpressStudy(
  counts = matrix(0L, nrow = 0, ncol = 0),
  metadata = data.frame(),
  accession = character(0),
  organism = "Homo sapiens",
  qc_score = NA_real_,
  de_result = data.frame(),
  gene_info = NULL
)
```

## Arguments

- counts:

  A numeric matrix of raw counts (genes x samples). Rownames must be
  gene identifiers; colnames must be sample identifiers.

- metadata:

  A `data.frame` of sample metadata, one row per column of `counts`.
  Must contain a `condition` column when `accession` is set.

- accession:

  A length-1 character accession ID, or `"local"`.

- organism:

  A length-1 character species name.

- qc_score:

  A length-1 numeric QC score, `NA` if not yet run.

- de_result:

  A `data.frame` of DE results (usually empty).

- gene_info:

  Optional `data.frame` or `DataFrame` of gene annotation, one row per
  row of `counts`; stored as `rowData`.

## Value

A
[`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
object.

## Examples

``` r
counts <- matrix(rpois(24, 100), nrow = 4,
                 dimnames = list(paste0("GENE", 1:4), paste0("S", 1:6)))
meta <- data.frame(sample_id = paste0("S", 1:6),
                   condition = rep(c("control", "case"), each = 3))
study <- metaXpressStudy(counts, meta, accession = "local")
study
#> metaXpressStudy
#>   Accession : local 
#>   Organism  : Homo sapiens 
#>   Genes     : 4 
#>   Samples   : 6 
#>   QC score  : not run 
#>   DE run    : FALSE 
```
