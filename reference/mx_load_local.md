# Load user-supplied count matrices as metaXpressStudy objects

Constructs `metaXpressStudy` objects from local count matrix files and
matching metadata files. Supports CSV, TSV, and RDS formats.

## Usage

``` r
mx_load_local(
  count_paths,
  metadata_paths,
  accessions = NULL,
  organism = "Homo sapiens"
)
```

## Arguments

- count_paths:

  Character vector of file paths to count matrices. Each file must have
  genes as rows (with rownames) and samples as columns.

- metadata_paths:

  Character vector of file paths to sample metadata tables. Each must
  contain at minimum `condition` and `sample_id` columns. Length must
  match `count_paths`.

- accessions:

  Character vector of accession IDs (labels), one per study. Defaults to
  the basename of `count_paths` without extension.

- organism:

  Character scalar. Species for all studies. Default: `"Homo sapiens"`.

## Value

A named list of
[`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
objects.

## See also

[`mx_fetch_geo`](https://mdjubayerhossain.com/metaXpress/reference/mx_fetch_geo.md),
[`mx_qc_study`](https://mdjubayerhossain.com/metaXpress/reference/mx_qc_study.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  studies <- mx_load_local(
    count_paths    = c("study1_counts.csv", "study2_counts.csv"),
    metadata_paths = c("study1_meta.csv",   "study2_meta.csv")
  )
} # }
```
