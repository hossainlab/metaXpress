# Fetch bulk RNA-seq studies from GEO

Downloads count matrices and sample metadata from the Gene Expression
Omnibus (GEO) for one or more accession IDs and returns a list of
`metaXpressStudy` objects ready for downstream analysis.

## Usage

``` r
mx_fetch_geo(
  accessions,
  count_type = "raw",
  cache_dir = tempdir(),
  BPPARAM = BiocParallel::bpparam()
)
```

## Arguments

- accessions:

  Character vector of GEO accession IDs (e.g.,
  `c("GSE12345", "GSE67890")`).

- count_type:

  Character scalar. One of `"raw"` (default) or `"normalized"`. Always
  prefer `"raw"` for meta-analysis.

- cache_dir:

  Character scalar. Directory for caching downloaded files. Defaults to
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object controlling parallelization. Default:
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html).
  Set to
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
  for sequential execution.

## Value

A named list of
[`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
objects, one per accession. QC scores are filled by an internal call to
[`mx_qc_study`](https://mdjubayerhossain.com/metaXpress/reference/mx_qc_study.md).
Studies failing QC receive a warning but are returned (use
[`mx_filter_studies`](https://mdjubayerhossain.com/metaXpress/reference/mx_filter_studies.md)
to remove them).

## References

Davis, S. & Meltzer, P.S. (2007) GEOquery: a bridge between the Gene
Expression Omnibus (GEO) and BioConductor. *Bioinformatics*, **23**(14),
1846–1847.
[doi:10.1093/bioinformatics/btm254](https://doi.org/10.1093/bioinformatics/btm254)

Heberle, H. et al. (2025) A comprehensive framework for quality control
and meta-analysis of bulk RNA-seq data. *Alzheimer's & Dementia*,
**21**(1), e70025.
[doi:10.1002/alz.70025](https://doi.org/10.1002/alz.70025)

## See also

[`mx_load_local`](https://mdjubayerhossain.com/metaXpress/reference/mx_load_local.md),
[`mx_qc_study`](https://mdjubayerhossain.com/metaXpress/reference/mx_qc_study.md),
[`mx_filter_studies`](https://mdjubayerhossain.com/metaXpress/reference/mx_filter_studies.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  studies <- mx_fetch_geo(c("GSE12345", "GSE67890"))
  studies <- mx_filter_studies(studies, qc_threshold = 7)
} # }
```
