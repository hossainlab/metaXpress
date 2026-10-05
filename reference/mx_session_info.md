# Capture session information for reproducibility

Returns a structured list containing the R session info, package
versions, and metaXpress analysis parameters for inclusion in reports
and supplementary materials.

## Usage

``` r
mx_session_info()
```

## Value

A list with elements `session_info` (from
[`utils::sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html)),
`timestamp`, and `metaXpress_version`.

## Examples

``` r
info <- mx_session_info()
info$timestamp
#> [1] "2026-10-05 15:21:56 UTC"
```
