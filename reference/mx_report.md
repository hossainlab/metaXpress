# Generate a reproducible analysis report

Renders a full HTML or PDF report summarising the meta-analysis
pipeline, including QC metrics, per-study DE summaries, meta-analysis
results, and pathway enrichment results.

## Usage

``` r
mx_report(
  meta_result,
  studies,
  de_results,
  format = c("html", "pdf"),
  output_dir = ".",
  output_file = "metaXpress_report"
)
```

## Arguments

- meta_result:

  A
  [`metaXpressResult`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressResult-class.md)
  object from
  [`mx_meta`](https://mdjubayerhossain.com/metaXpress/reference/mx_meta.md).

- studies:

  A named list of
  [`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
  objects.

- de_results:

  A named list of per-study DE `data.frame`s from
  [`mx_de_all`](https://mdjubayerhossain.com/metaXpress/reference/mx_de_all.md).

- format:

  Character scalar. Output format. One of `"html"` (default) or `"pdf"`.

- output_dir:

  Character scalar. Directory for the output report. Default: current
  working directory (`"."`).

- output_file:

  Character scalar. Filename (without extension). Default:
  `"metaXpress_report"`.

## Value

The path to the generated report file (invisibly).

## See also

[`mx_export`](https://mdjubayerhossain.com/metaXpress/reference/mx_export.md),
[`mx_session_info`](https://mdjubayerhossain.com/metaXpress/reference/mx_session_info.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  mx_report(meta_result, studies, de_results,
            format = "html", output_dir = "results/")
} # }
```
