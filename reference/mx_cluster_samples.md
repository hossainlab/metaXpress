# Suggest sample groups from metadata text

Exploratory helper for GEO series whose phenotype columns do not state
the case/control group cleanly. Concatenates each sample's text
metadata, computes Jaro-Winkler string distances and cuts an
average-linkage hierarchical clustering into `k` groups. The result is
written to a new `cluster` column; `condition` is **not** changed,
because automatically derived labels must be checked by a person before
they define a contrast.

## Usage

``` r
mx_cluster_samples(study, k = 2, columns = NULL)
```

## Arguments

- study:

  A
  [`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
  object.

- k:

  Integer. Number of groups. Default: `2`.

- columns:

  Optional character vector of `colData` columns to use. Default: all
  character/factor columns except identifiers.

## Value

The input `study` with `colData(study)$cluster`. A cross-tabulation
against `condition` (if present) is printed.

## Examples

``` r
counts <- matrix(rpois(60, 50), nrow = 10,
                 dimnames = list(paste0("G", 1:10), paste0("S", 1:6)))
meta <- data.frame(title = c(paste("healthy donor", 1:3),
                             paste("tumour patient", 1:3)))
if (requireNamespace("stringdist", quietly = TRUE)) {
  s <- mx_cluster_samples(metaXpressStudy(counts, meta))
  s$cluster
}
#> Assigned 6 samples to 2 clusters in 'cluster'. Inspect before using as condition labels.
#> [1] "cluster1" "cluster1" "cluster1" "cluster2" "cluster2" "cluster2"
```
