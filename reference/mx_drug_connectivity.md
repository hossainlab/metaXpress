# Rank drugs by how strongly they reverse the meta-signature

Scores each drug perturbation signature (sets of genes the drug pushes
up and down) against the meta-analysis ranking. The reversal score is a
standardised difference of the mean signed meta z-score of the drug's up
and down genes; strongly negative values mean the drug moves genes in
the opposite direction to the disease and are candidate therapeutics.
The Connectivity Map Kolmogorov-Smirnov score (Lamb et al. 2006) is
reported alongside for comparison.

## Usage

``` r
mx_drug_connectivity(
  meta_result,
  signatures,
  id_col = "gene_id",
  min_size = 10
)
```

## Arguments

- meta_result:

  A
  [`metaXpressResult`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressResult-class.md)
  object.

- signatures:

  A named list of drug signatures, each a list with character vectors
  `up` and `down`, as returned by
  [`mx_read_signatures`](https://mdjubayerhossain.com/metaXpress/reference/mx_read_signatures.md).

- id_col:

  Column of `meta_table` whose values match the signature genes. Default
  `"gene_id"`.

- min_size:

  Minimum number of measured genes in each of the up and down sets.
  Default 10.

## Value

A `data.frame` sorted by reversal score with columns `drug`,
`reversal_score` (negative = reverses the disease signature), `pvalue`,
`padj`, `cmap_score`, `n_up`, `n_down`. P-values assume independent
genes and are best used for ranking rather than as calibrated error
rates.

## References

Lamb, J. et al. (2006) The Connectivity Map: using gene-expression
signatures to connect small molecules, genes, and disease. *Science*,
**313**(5795), 1929–1935.
[doi:10.1126/science.1132939](https://doi.org/10.1126/science.1132939)

Subramanian, A. et al. (2017) A next generation Connectivity Map: L1000
platform and the first 1,000,000 profiles. *Cell*, **171**(6),
1437–1452.
[doi:10.1016/j.cell.2017.10.049](https://doi.org/10.1016/j.cell.2017.10.049)

## See also

[`mx_read_signatures`](https://mdjubayerhossain.com/metaXpress/reference/mx_read_signatures.md),
[`mx_downstream_plot`](https://mdjubayerhossain.com/metaXpress/reference/mx_downstream_plot.md)

## Examples

``` r
mt <- data.frame(gene_id = paste0("G", 1:300),
                 meta_log2FC = c(rep(2, 30), rep(-2, 30), rnorm(240, sd = 0.2)),
                 meta_pvalue = c(rep(1e-8, 60), runif(240)),
                 meta_padj = 1, i_squared = 0, n_studies = 3L)
res <- new("metaXpressResult", meta_table = mt, method = "fisher",
           n_studies = 3L)
sigs <- list(reverser = list(up = paste0("G", 31:60), down = paste0("G", 1:30)),
             mimic    = list(up = paste0("G", 1:30), down = paste0("G", 31:60)))
mx_drug_connectivity(res, sigs)
#>       drug reversal_score       pvalue         padj cmap_score n_up n_down
#> 1 reverser      -16.43202 1.128432e-60 1.128432e-60  -1.803333   30     30
#> 2    mimic       16.43202 1.128432e-60 1.128432e-60   1.803333   30     30
```
