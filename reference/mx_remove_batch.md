# Remove between-study batch effects for pooled expression analyses

Combines all studies into one expression matrix and removes study-level
batch effects while preserving the biological `condition` effect.

## Usage

``` r
mx_remove_batch(
  studies,
  method = c("ComBat-seq", "ComBat", "limma"),
  covariate = "condition"
)
```

## Arguments

- studies:

  A named list of
  [`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
  objects.

- method:

  Character scalar. `"ComBat-seq"` (default; negative binomial model on
  counts), `"ComBat"` (empirical Bayes on log2 CPM) or `"limma"`
  (`removeBatchEffect` on log2 CPM).

- covariate:

  Name of a `colData` column whose effect must be preserved. Default:
  `"condition"`. Use `NULL` only if you understand that condition
  effects confounded with study will be removed.

## Value

The input list of studies, restricted to their common genes, with the
corrected values in a new assay: `"combat_seq_counts"` for ComBat-seq or
`"log2cpm_batch_corrected"` otherwise. The raw `"counts"` assay is
untouched.

## Details

**Batch correction is not required, and should not be used, before
per-study DE**:
[`mx_de`](https://mdjubayerhossain.com/metaXpress/reference/mx_de.md)
analyses each study separately, so between-study batch effects never
enter the per-study fold changes, and feeding batch-corrected values
into count models invalidates their variance assumptions. Use this
function only for pooled visualizations (PCA, heatmaps), consensus
co-expression or other analyses on a merged matrix. Batch is defined as
study of origin.

## References

Johnson, W.E., Li, C. & Rabinovic, A. (2007) Adjusting batch effects in
microarray expression data using empirical Bayes methods.
*Biostatistics*, **8**(1), 118–127.
[doi:10.1093/biostatistics/kxj037](https://doi.org/10.1093/biostatistics/kxj037)

Zhang, Y., Parmigiani, G. & Johnson, W.E. (2020) ComBat-seq: batch
effect adjustment for RNA-seq count data. *NAR Genomics and
Bioinformatics*, **2**(3), lqaa078.
[doi:10.1093/nargab/lqaa078](https://doi.org/10.1093/nargab/lqaa078)

Nygaard, V., Rødland, E.A. & Hovig, E. (2016) Methods that remove batch
effects while retaining group differences may lead to exaggerated
confidence in downstream analyses. *Biostatistics*, **17**(1), 29–39.
[doi:10.1093/biostatistics/kxv027](https://doi.org/10.1093/biostatistics/kxv027)

## See also

[`mx_normalize`](https://mdjubayerhossain.com/metaXpress/reference/mx_normalize.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  pooled <- mx_remove_batch(studies, method = "limma")
  SummarizedExperiment::assay(pooled[[1]], "log2cpm_batch_corrected")
} # }
```
