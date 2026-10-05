# Correct for library type bias (polyA vs rRNA-depleted)

Adjusts for systematic expression differences between polyA-selected and
rRNA-depleted RNA-seq libraries using the ratio-based correction
described by Bush et al. (2017).

## Usage

``` r
mx_correct_library_type(studies)
```

## Arguments

- studies:

  A named list of
  [`metaXpressStudy`](https://hossainlab.github.io/metaXpress/reference/metaXpressStudy-class.md)
  objects. Studies must have a `library_type` column in `metadata` with
  values `"polyA"` or `"rRNA_depleted"`.

## Value

The input list of studies with library type bias corrected.

## References

Bush, S.J. et al. (2017) Cross-species inference of long non-coding RNAs
greatly expands the ruminant transcriptome. *Genetics Selection
Evolution*, **50**(1), 20.
[doi:10.1186/s12859-017-1714-9](https://doi.org/10.1186/s12859-017-1714-9)

## Examples

``` r
if (FALSE) { # \dontrun{
  studies <- mx_correct_library_type(studies)
} # }
```
