# Leave-one-out sensitivity analysis

Reruns the meta-analysis iteratively, each time excluding one study, to
assess whether any single study drives the overall result. This is the
standard sensitivity analysis recommended by the Cochrane Handbook
(Higgins et al. 2019, Section 10.14).

## Usage

``` r
mx_sensitivity(
  de_results,
  method = c("random_effects", "fisher", "stouffer", "inverse_normal", "fixed_effects",
    "awmeta")
)
```

## Arguments

- de_results:

  A named list of `data.frame` objects from
  [`mx_de_all`](https://hossainlab.github.io/metaXpress/reference/mx_de_all.md).

- method:

  Character scalar. Meta-analysis method. Default: `"random_effects"`.

## Value

A list of
[`metaXpressResult`](https://hossainlab.github.io/metaXpress/reference/metaXpressResult-class.md)
objects, one per leave-one-out iteration. Names correspond to the
excluded study.

## References

Higgins, J.P.T. et al. (eds.) (2019) *Cochrane Handbook for Systematic
Reviews of Interventions*. 2nd edn. Chichester: John Wiley & Sons.
[doi:10.1002/9781119536604](https://doi.org/10.1002/9781119536604)

## See also

[`mx_meta`](https://hossainlab.github.io/metaXpress/reference/mx_meta.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  loo <- mx_sensitivity(de_results, method = "random_effects")
  # Compare n significant genes across LOO runs
  vapply(loo, function(r) sum(r@meta_table$meta_padj <= 0.05, na.rm = TRUE),
         integer(1))
} # }
```
