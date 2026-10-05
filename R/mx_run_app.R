#' Launch the metaXpress Interactive Explorer Shiny Application
#'
#' Opens the interactive web dashboard for exploring multi-study bulk RNA-seq
#' meta-analysis pipelines, visualizing volcano plots, generating gene forest
#' plots, and evaluating meta-analysis statistics.
#'
#' @param launch.browser Logical. Whether to open the application automatically
#'   in the default web browser. Default: \code{TRUE}.
#' @param port Integer. Optional network port to listen on.
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}.
#'
#' @return No return value, called for side effects (starts local web server).
#'
#' @examples
#' \dontrun{
#'   mx_run_app()
#' }
#'
#' @export
mx_run_app <- function(launch.browser = TRUE, port = NULL, ...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("The 'shiny' package is required to run the interactive explorer. ",
         "Please install it with install.packages('shiny').")
  }

  app_dir <- system.file("shiny", package = "metaXpress")
  if (app_dir == "" || !dir.exists(app_dir)) {
    # Fallback to local development path if running from source
    if (dir.exists("inst/shiny")) {
      app_dir <- "inst/shiny"
    } else {
      stop("Could not locate the Shiny app directory in package 'metaXpress'.")
    }
  }

  shiny::runApp(appDir = app_dir, launch.browser = launch.browser,
                port = port, ...)
}
