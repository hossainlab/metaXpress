# Remove batch effects across studies

Applies one of four batch effect removal methods to harmonize expression
levels across studies before meta-analysis.

## Usage

``` r
mx_remove_batch(
  studies,
  method = c("ComBat-seq", "ComBat", "limma", "harmony")
)
```

## Arguments

- studies:

  A named list of
  [`metaXpressStudy`](https://hossainlab.github.io/metaXpress/reference/metaXpressStudy-class.md)
  objects.

- method:

  Character scalar. Batch removal method. One of `"ComBat-seq"`
  (default), `"ComBat"`, `"limma"`, or `"harmony"`.

## Value

The input list of studies with batch effects removed from `counts`.

## Details

Batch is defined as study of origin. `"ComBat-seq"` operates on raw
counts and is preferred for count data. `"ComBat"` and `"limma"` operate
on log-transformed values. `"harmony"` operates on PCA projections of
the data, which is computationally efficient.

## References

Johnson, W.E., Li, C. & Rabinovic, A. (2007) Adjusting batch effects in
microarray expression data using empirical Bayes methods.
*Biostatistics*, **8**(1), 118–127.
[doi:10.1093/biostatistics/kxj037](https://doi.org/10.1093/biostatistics/kxj037)

Zhang, Y. et al. (2020) ComBat-seq: batch effect adjustment for RNA-seq
count data. *NAR Genomics and Bioinformatics*, **2**(3), lqaa078.
[doi:10.1093/nargab/lqaa078](https://doi.org/10.1093/nargab/lqaa078)

Ritchie, M.E. et al. (2015) limma powers differential expression
analyses for RNA-sequencing and microarray studies. *Nucleic Acids
Research*, **43**(7), e47.
[doi:10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007)

Korsunsky, I. et al. (2019) Fast, sensitive and accurate integration of
single-cell data with Harmony. *Nature Methods*, **16**(12), 1289-1296.
[doi:10.1038/s41592-019-0619-0](https://doi.org/10.1038/s41592-019-0619-0)

## See also

[`mx_normalize`](https://hossainlab.github.io/metaXpress/reference/mx_normalize.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  studies <- mx_remove_batch(studies, method = "ComBat-seq")
} # }
```
