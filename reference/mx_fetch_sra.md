# Fetch raw RNA-seq data from SRA

Downloads pre-quantified count matrices for SRA project IDs using the
\pkgrecount3 package. This avoids computationally expensive local
alignment while ensuring uniform processing across studies
(Collado-Torres et al. 2017; Wilks et al. 2021).

## Usage

``` r
mx_fetch_sra(srp_ids, cache_dir = tempdir(), organism = "human")
```

## Arguments

- srp_ids:

  Character vector of SRA project IDs (e.g., `"SRP123456"`).

- cache_dir:

  Character scalar. Directory for caching downloaded files.

- organism:

  Character scalar. Species for the studies. Default: `"human"`. (Note:
  recount3 uses `"human"` or `"mouse"`).

## Value

A named list of
[`metaXpressStudy`](https://hossainlab.github.io/metaXpress/reference/metaXpressStudy-class.md)
objects.

## References

Wilks, C. et al. (2021) recount3: summaries and queries for large-scale
RNA-seq expression and splicing. *Genome Biology*, **22**(1), 323.
[doi:10.1186/s13059-021-02533-6](https://doi.org/10.1186/s13059-021-02533-6)

## See also

[`mx_fetch_geo`](https://hossainlab.github.io/metaXpress/reference/mx_fetch_geo.md),
[`mx_load_local`](https://hossainlab.github.io/metaXpress/reference/mx_load_local.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  studies <- mx_fetch_sra("SRP123456")
} # }
```
