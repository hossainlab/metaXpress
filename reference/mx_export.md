# Export meta-analysis results to file

Saves the gene-level meta-analysis table from a
[`metaXpressResult`](https://hossainlab.github.io/metaXpress/reference/metaXpressResult-class.md)
object to CSV, Excel, or RDS format.

## Usage

``` r
mx_export(
  meta_result,
  format = c("csv", "excel", "rds"),
  output_dir = ".",
  prefix = "metaXpress_results"
)
```

## Arguments

- meta_result:

  A
  [`metaXpressResult`](https://hossainlab.github.io/metaXpress/reference/metaXpressResult-class.md)
  object.

- format:

  Character scalar. Export format. One of `"csv"` (default), `"excel"`,
  or `"rds"`.

- output_dir:

  Character scalar. Directory for output files. Default: current working
  directory (`"."`).

- prefix:

  Character scalar. Filename prefix. Default: `"metaXpress_results"`.

## Value

The path to the exported file (invisibly).

## See also

[`mx_report`](https://hossainlab.github.io/metaXpress/reference/mx_report.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  mx_export(meta_result, format = "csv", output_dir = "results/")
  mx_export(meta_result, format = "excel")
} # }
```
