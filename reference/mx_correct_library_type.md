# Correct for library type bias (polyA vs rRNA-depleted)

Adjusts for systematic expression differences between polyA-selected and
rRNA-depleted RNA-seq libraries with a per-gene ratio correction in the
spirit of Bush et al. (2017). The ratio is estimated from *control*
samples only, so disease effects do not leak into it.

## Usage

``` r
mx_correct_library_type(studies)
```

## Arguments

- studies:

  A named list of
  [`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
  objects. Studies must have a `library_type` column in `metadata` with
  values `"polyA"` or `"rRNA_depleted"`.

## Value

The input list of studies, restricted to their common genes, with an
added `"cpm_libtype_corrected"` assay in rRNA-depleted studies (and
plain CPM in polyA studies).

## Details

Per-gene multiplicative scaling cancels out of within-study log fold
changes, so this correction is **not needed** for the standard per-study
DE + meta-analysis workflow. It matters only when expression levels are
compared across studies (pooled heatmaps, co-expression,
cell-composition scores). The corrected values are stored in a new
`"cpm_libtype_corrected"` assay; raw counts are untouched.

## References

Bush, S.J. et al. (2017) Integration of quantitated expression estimates
from polyA-selected and rRNA-depleted RNA-seq libraries. *BMC
Bioinformatics*, **18**, 301.
[doi:10.1186/s12859-017-1714-9](https://doi.org/10.1186/s12859-017-1714-9)

## Examples

``` r
if (FALSE) { # \dontrun{
  studies <- mx_correct_library_type(studies)
} # }
```
