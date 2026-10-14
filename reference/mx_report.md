# Generate a reproducible analysis report

Renders an HTML (or PDF) report with methods text, study table,
fold-change concordance, volcano plot, top-gene table with confidence
and prediction intervals, heterogeneity, forest plots, leave-one-out
robustness, pathway results (if present) and session information.

## Usage

``` r
mx_report(
  meta_result,
  de_results,
  studies = NULL,
  format = c("html", "pdf"),
  output_dir = ".",
  output_file = "metaXpress_report",
  padj_threshold = 0.05,
  lfc_threshold = 0
)
```

## Arguments

- meta_result:

  A
  [`metaXpressResult`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressResult-class.md)
  object from
  [`mx_meta`](https://mdjubayerhossain.com/metaXpress/reference/mx_meta.md).

- de_results:

  The named list of per-study DE `data.frame`s (or
  [`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
  objects) used for the meta-analysis.

- studies:

  Optional list of
  [`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
  objects (currently unused; kept for API compatibility).

- format:

  `"html"` (default) or `"pdf"` (requires LaTeX).

- output_dir:

  Output directory. Default: `"."`.

- output_file:

  File name without extension. Default: `"metaXpress_report"`.

- padj_threshold, lfc_threshold:

  Significance thresholds used in the report. Defaults: `0.05` and `0`.

## Value

The path to the generated report (invisibly).

## See also

[`mx_export`](https://mdjubayerhossain.com/metaXpress/reference/mx_export.md),
[`mx_session_info`](https://mdjubayerhossain.com/metaXpress/reference/mx_session_info.md)

## Examples

``` r
# \donttest{
if (rmarkdown::pandoc_available()) {
  ad <- readRDS(system.file("extdata", "ad_brain_de.rds",
                            package = "metaXpress"))
  res <- mx_meta(ad$de_results)
  mx_report(res, ad$de_results, output_dir = tempdir())
}
#> Running meta-analysis (random_effects) on 4 studies...
#> Report written to: /tmp/RtmpCvQhZS/metaXpress_report.html
# }
```
