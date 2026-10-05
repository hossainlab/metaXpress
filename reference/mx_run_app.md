# Launch the metaXpress Interactive Explorer Shiny Application

Opens the interactive web dashboard for exploring multi-study bulk
RNA-seq meta-analysis pipelines, visualizing volcano plots, generating
gene forest plots, and evaluating meta-analysis statistics.

## Usage

``` r
mx_run_app(launch.browser = TRUE, port = NULL, ...)
```

## Arguments

- launch.browser:

  Logical. Whether to open the application automatically in the default
  web browser. Default: `TRUE`.

- port:

  Integer. Optional network port to listen on.

- ...:

  Additional arguments passed to
  [`runApp`](https://rdrr.io/pkg/shiny/man/runApp.html).

## Value

No return value, called for side effects (starts local web server).

## Examples

``` r
if (FALSE) { # \dontrun{
  mx_run_app()
} # }
```
