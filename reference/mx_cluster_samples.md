# Automatically cluster samples by metadata

Uses text mining on metadata fields to assign samples to case/control
groups, mimicking the \pkgsampleclusteR approach for automated group
assignment from GEO metadata. It computes string distances between
sample descriptions and uses hierarchical clustering to partition them
into two groups.

## Usage

``` r
mx_cluster_samples(study)
```

## Arguments

- study:

  A
  [`metaXpressStudy`](https://hossainlab.github.io/metaXpress/reference/metaXpressStudy-class.md)
  object.

## Value

The input `study` with the `condition` column updated in `metadata`.

## References

Coke, T., Niranjan, M. & Ewing, R.M. (2025) sampleclusteR: automated
case/control group assignment from GEO metadata. *bioRxiv*.
[doi:10.1101/2025.04.10.648129](https://doi.org/10.1101/2025.04.10.648129)

## Examples

``` r
if (FALSE) { # \dontrun{
  study <- mx_cluster_samples(study)
} # }
```
