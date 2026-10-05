# Run meta-analysis across per-study DE results

Integrates per-study differential expression results from
[`mx_de_all`](https://mdjubayerhossain.com/metaXpress/reference/mx_de_all.md)
into a single meta-analysis result using one of six statistical methods.

## Usage

``` r
mx_meta(
  de_results,
  method = c("random_effects", "fisher", "stouffer", "inverse_normal", "fixed_effects",
    "awmeta"),
  min_studies = 2,
  alpha = 0.05,
  n_samples = NULL,
  ...
)
```

## Arguments

- de_results:

  A named list of `data.frame` objects from
  [`mx_de_all`](https://mdjubayerhossain.com/metaXpress/reference/mx_de_all.md)
  (i.e., the `de_result` slots extracted from each study), or a list of
  [`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
  objects.

- method:

  Character scalar. Meta-analysis method. One of:

  `"fisher"`

  :   Fisher's combined p-value (Rau et al. 2013)

  `"stouffer"`

  :   Stouffer's Z-score weighted by sample size

  `"inverse_normal"`

  :   Fused inverse-normal (Prasad & Li 2021)

  `"fixed_effects"`

  :   Fixed effects, inverse-variance weighting

  `"random_effects"`

  :   DerSimonian-Laird random effects (default)

  `"awmeta"`

  :   Adaptive weighting (Hu et al. 2025)

- min_studies:

  Integer. Minimum number of studies a gene must appear in to be
  included. Default: `2`.

- alpha:

  Numeric. Significance level for the BH FDR correction. Default:
  `0.05`.

- n_samples:

  Integer vector. Sample sizes (total n per study) used for Stouffer
  weighting. If `NULL`, uniform weights are used.

- ...:

  Additional arguments (currently unused).

## Value

A
[`metaXpressResult`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressResult-class.md)
object.

## Details

**CRITICAL:** Use `meta_padj` (BH-adjusted meta p-value) for
significance filtering, not `meta_pvalue`. Default significance
threshold: `meta_padj <= 0.05` AND `|meta_log2FC| >= 1`.

## References

Fisher, R.A. (1932) *Statistical Methods for Research Workers*. 4th edn.
Edinburgh: Oliver & Boyd.

Stouffer, S.A. et al. (1949) *The American Soldier: Adjustment During
Army Life*. Princeton University Press.

DerSimonian, R. & Laird, N. (1986) Meta-analysis in clinical trials.
*Controlled Clinical Trials*, **7**(3), 177–188.
[doi:10.1016/0197-2456(86)90046-2](https://doi.org/10.1016/0197-2456%2886%2990046-2)

Cochran, W.G. (1954) The combination of estimates from different
experiments. *Biometrics*, **10**(1), 101–129.
[doi:10.2307/3001666](https://doi.org/10.2307/3001666)

Higgins, J.P.T. & Thompson, S.G. (2002) Quantifying heterogeneity in a
meta-analysis. *Statistics in Medicine*, **21**(11), 1539–1558.
[doi:10.1002/sim.1186](https://doi.org/10.1002/sim.1186)

Rau, A., Marot, G. & Jaffrézic, F. (2014) Differential meta-analysis of
RNA-seq data from multiple studies. *BMC Bioinformatics*, **15**, 91.
[doi:10.1186/1471-2105-15-91](https://doi.org/10.1186/1471-2105-15-91)

Prasad, A. & Li, Q. (2022) Fused inverse-normal method for integrated
differential expression analysis of RNA-seq data. *BMC Bioinformatics*,
**23**, 371.
[doi:10.1186/s12859-022-04859-9](https://doi.org/10.1186/s12859-022-04859-9)

Keel, B.N. & Lindholm-Perry, A.K. (2022) Recent developments and future
directions in meta-analysis of differential gene expression. *Frontiers
in Genetics*, **13**, 983043.
[doi:10.3389/fgene.2022.983043](https://doi.org/10.3389/fgene.2022.983043)

Hu, J. et al. (2025) AWmeta: adaptive weighting meta-analysis for
combining heterogeneous genomic studies. *bioRxiv*.
[doi:10.1101/2025.05.06.650408](https://doi.org/10.1101/2025.05.06.650408)

## See also

[`mx_heterogeneity`](https://mdjubayerhossain.com/metaXpress/reference/mx_heterogeneity.md),
[`mx_sensitivity`](https://mdjubayerhossain.com/metaXpress/reference/mx_sensitivity.md),
[`mx_forest`](https://mdjubayerhossain.com/metaXpress/reference/mx_forest.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  de_list <- lapply(studies, function(s) s@de_result)
  result  <- mx_meta(de_list, method = "random_effects")
  head(result@meta_table[order(result@meta_table$meta_padj), ])
} # }
```
