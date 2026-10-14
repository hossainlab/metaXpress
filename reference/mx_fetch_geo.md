# Fetch bulk RNA-seq studies from GEO

Downloads raw count matrices and sample metadata from the Gene
Expression Omnibus (GEO) for one or more series and returns a list of
`metaXpressStudy` objects ready for downstream analysis.

## Usage

``` r
mx_fetch_geo(
  accessions,
  source = c("ncbi", "supplementary"),
  count_type = "raw",
  condition_col = NULL,
  annotate = TRUE,
  cache_dir = tempdir(),
  BPPARAM = NULL
)
```

## Arguments

- accessions:

  Character vector of GEO series accessions (e.g.,
  `c("GSE52778", "GSE95587")`).

- source:

  Character scalar. `"ncbi"` (default) uses NCBI GEO-generated raw
  counts with fallback to supplementary files; `"supplementary"` only
  uses author-supplied files.

- count_type:

  Character scalar. One of `"raw"` (default) or `"normalized"`. Always
  use `"raw"` for count-based DE.

- condition_col:

  Optional. Name of the GEO phenotype column holding the case/control
  labels (e.g. `"disease state:ch1"`), either one name for all series or
  a character vector named by accession. When `NULL`, a column is
  auto-detected; always check the result.

- annotate:

  Logical. Add Entrez ID, gene symbol and gene type from NCBI
  `gene_info` to `rowData` for NCBI counts. Default: `TRUE`.

- cache_dir:

  Character scalar. Directory for caching downloaded files. Defaults to
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html); use a persistent
  directory to avoid repeated downloads.

- BPPARAM:

  Unused; downloads run serially to respect NCBI rate limits. Kept for
  backward compatibility.

## Value

A named list of
[`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
objects, one per accession. `metadata(study)$count_source` records where
the counts came from. QC scores are filled by
[`mx_qc_study`](https://mdjubayerhossain.com/metaXpress/reference/mx_qc_study.md).

## Details

By default (`source = "ncbi"`) counts come from the NCBI GEO RNA-seq
pipeline, which re-processes the SRA reads of every eligible human and
mouse series with a single pipeline (HISAT2 alignment to GRCh38.p13 /
GRCm39 and featureCounts against NCBI RefSeq annotation). Because every
study is quantified identically, these counts remove the largest source
of technical heterogeneity in a meta-analysis (different aligners,
annotations and count definitions). Genes are NCBI Entrez Gene IDs.
Series without NCBI counts fall back to author-supplied supplementary
files.

## References

Davis, S. & Meltzer, P.S. (2007) GEOquery: a bridge between the Gene
Expression Omnibus (GEO) and BioConductor. *Bioinformatics*, **23**(14),
1846–1847.
[doi:10.1093/bioinformatics/btm254](https://doi.org/10.1093/bioinformatics/btm254)

Clough, E. et al. (2024) NCBI GEO: archive for gene expression and
epigenomics data sets: 23-year update. *Nucleic Acids Research*,
**52**(D1), D138–D144.
[doi:10.1093/nar/gkad965](https://doi.org/10.1093/nar/gkad965)

## See also

[`mx_load_local`](https://mdjubayerhossain.com/metaXpress/reference/mx_load_local.md),
[`mx_qc_study`](https://mdjubayerhossain.com/metaXpress/reference/mx_qc_study.md),
[`mx_filter_studies`](https://mdjubayerhossain.com/metaXpress/reference/mx_filter_studies.md)

## Examples

``` r
# \donttest{
# Requires internet access. Himes et al. (2014) airway smooth muscle,
# NCBI-processed counts for 16 samples.
if (requireNamespace("GEOquery", quietly = TRUE)) {
  studies <- mx_fetch_geo("GSE52778",
                          condition_col = "treatment:ch1")
  studies$GSE52778
  table(studies$GSE52778$condition)
}
#> Setting options('download.file.method.GEOquery'='auto')
#> Setting options('GEOquery.inmemory.gpl'=FALSE)
#> Fetching 1 studies from GEO...
#>   Downloading: GSE52778
#> Found 1 file(s)
#> GSE52778_series_matrix.txt.gz
#>   Using NCBI-generated raw counts (39376 genes x 16 samples)
#>   Condition from 'treatment:ch1': Albuterol=4, Albuterol_Dexamethasone=4, Dexamethasone=4, Untreated=4
#> Warning: Study GSE52778 failed QC: clear contrast. See S4Vectors::metadata(study)$qc.
#> 
#>               Albuterol Albuterol_Dexamethasone           Dexamethasone 
#>                       4                       4                       4 
#>               Untreated 
#>                       4 
# }
```
