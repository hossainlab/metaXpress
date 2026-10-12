# Convert a SummarizedExperiment to a metaXpressStudy

Wraps any `SummarizedExperiment` (including `DESeqDataSet` and
`RangedSummarizedExperiment` objects, e.g. from recount3 or tximeta) as
a
[`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
so it can enter the metaXpress pipeline. `as(se, "metaXpressStudy")`
calls this with defaults.

## Usage

``` r
mx_from_se(
  se,
  assay_name = "counts",
  condition_col = "condition",
  accession = NULL,
  organism = NULL
)
```

## Arguments

- se:

  A `SummarizedExperiment` or subclass.

- assay_name:

  Name of the assay holding raw counts. Defaults to `"counts"`, falling
  back to the first assay if absent.

- condition_col:

  Name of the `colData` column holding case/control labels; copied to
  `condition` if different.

- accession:

  Accession ID. Defaults to `metadata(se)$accession` if present,
  otherwise `"local"`.

- organism:

  Species name. Defaults to `metadata(se)$organism` if present,
  otherwise `"Homo sapiens"`.

## Value

A
[`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
object.

## See also

[`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy.md)

## Examples

``` r
se <- SummarizedExperiment::SummarizedExperiment(
  assays  = list(counts = matrix(rpois(24, 100), nrow = 4,
    dimnames = list(paste0("GENE", 1:4), paste0("S", 1:6)))),
  colData = data.frame(dex = rep(c("untrt", "trt"), each = 3))
)
study <- mx_from_se(se, condition_col = "dex")
```
