#' @include AllClasses.R
NULL

# ============================================================================
# Module 8: Report Generation & Export
# ============================================================================

#' Generate a reproducible analysis report
#'
#' Renders a full HTML or PDF report summarising the meta-analysis pipeline,
#' including QC metrics, per-study DE summaries, meta-analysis results, and
#' pathway enrichment results.
#'
#' @param meta_result A \code{\linkS4class{metaXpressResult}} object from
#'   \code{\link{mx_meta}}.
#' @param studies A named list of \code{\linkS4class{metaXpressStudy}} objects.
#' @param de_results A named list of per-study DE \code{data.frame}s from
#'   \code{\link{mx_de_all}}.
#' @param format Character scalar. Output format. One of \code{"html"}
#'   (default) or \code{"pdf"}.
#' @param output_dir Character scalar. Directory for the output report.
#'   Default: current working directory (\code{"."}).
#' @param output_file Character scalar. Filename (without extension).
#'   Default: \code{"metaXpress_report"}.
#'
#' @return The path to the generated report file (invisibly).
#'
#' @examples
#' \dontrun{
#'   mx_report(meta_result, studies, de_results,
#'             format = "html", output_dir = "results/")
#' }
#'
#' @seealso \code{\link{mx_export}}, \code{\link{mx_session_info}}
#'
#' @importFrom methods is
#' @export
mx_report <- function(meta_result, studies, de_results,
                       format      = c("html", "pdf"),
                       output_dir  = ".",
                       output_file = "metaXpress_report") {
  if (!is(meta_result, "metaXpressResult"))
    stop("'meta_result' must be a metaXpressResult object")
  format <- match.arg(format)

  if (!requireNamespace("rmarkdown", quietly = TRUE))
    stop("Package 'rmarkdown' is required for report generation. ",
         "Install with: install.packages('rmarkdown')")

  template <- system.file("rmarkdown", "templates", "report.Rmd",
                           package = "metaXpress")
  if (!file.exists(template))
    stop("Report template not found. Package may not be fully installed.")

  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)

  output_format <- if (format == "html")
    rmarkdown::html_document(toc = TRUE, toc_float = TRUE,
                              code_folding = "hide")
  else
    rmarkdown::pdf_document(toc = TRUE)

  out_path <- file.path(output_dir,
                         paste0(output_file, ".", format))

  rmarkdown::render(
    input         = template,
    output_format = output_format,
    output_file   = out_path,
    params        = list(
      meta_result = meta_result,
      studies     = studies,
      de_results  = de_results
    ),
    envir = new.env(parent = globalenv())
  )

  message("Report written to: ", out_path)
  invisible(out_path)
}

#' Export meta-analysis results to file
#'
#' Saves the gene-level meta-analysis table from a
#' \code{\linkS4class{metaXpressResult}} object to CSV, Excel, or RDS format.
#'
#' @param meta_result A \code{\linkS4class{metaXpressResult}} object.
#' @param format Character scalar. Export format. One of \code{"csv"}
#'   (default), \code{"excel"}, or \code{"rds"}.
#' @param output_dir Character scalar. Directory for output files.
#'   Default: current working directory (\code{"."}).
#' @param prefix Character scalar. Filename prefix. Default:
#'   \code{"metaXpress_results"}.
#'
#' @return The path to the exported file (invisibly).
#'
#' @examples
#' \dontrun{
#'   mx_export(meta_result, format = "csv", output_dir = "results/")
#'   mx_export(meta_result, format = "excel")
#' }
#'
#' @seealso \code{\link{mx_report}}
#'
#' @importFrom methods is
#' @export
mx_export <- function(meta_result,
                       format     = c("csv", "excel", "rds"),
                       output_dir = ".",
                       prefix     = "metaXpress_results") {
  if (!is(meta_result, "metaXpressResult"))
    stop("'meta_result' must be a metaXpressResult object")
  format <- match.arg(format)

  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)

  out_path <- switch(format,
    csv   = {
      path <- file.path(output_dir, paste0(prefix, ".csv"))
      utils::write.csv(meta_result@meta_table, path, row.names = FALSE)
      path
    },
    excel = {
      if (!requireNamespace("openxlsx", quietly = TRUE))
        stop("Package 'openxlsx' is required for Excel export. ",
             "Install with: install.packages('openxlsx')")
      path <- file.path(output_dir, paste0(prefix, ".xlsx"))
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "meta_results")
      openxlsx::writeData(wb, "meta_results", meta_result@meta_table)
      if (nrow(meta_result@heterogeneity) > 0) {
        openxlsx::addWorksheet(wb, "heterogeneity")
        openxlsx::writeData(wb, "heterogeneity", meta_result@heterogeneity)
      }
      if (nrow(meta_result@pathway_result) > 0) {
        openxlsx::addWorksheet(wb, "pathway_results")
        openxlsx::writeData(wb, "pathway_results",
                             meta_result@pathway_result)
      }
      openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
      path
    },
    rds   = {
      path <- file.path(output_dir, paste0(prefix, ".rds"))
      saveRDS(meta_result, path)
      path
    }
  )

  message("Results exported to: ", out_path)
  invisible(out_path)
}

#' Capture session information for reproducibility
#'
#' Returns a structured list containing the R session info, package versions,
#' and metaXpress analysis parameters for inclusion in reports and
#' supplementary materials.
#'
#' @return A list with elements \code{session_info} (from
#'   \code{utils::sessionInfo()}), \code{timestamp}, and
#'   \code{metaXpress_version}.
#'
#' @examples
#' info <- mx_session_info()
#' info$timestamp
#'
#' @export
mx_session_info <- function() {
  list(
    session_info      = utils::sessionInfo(),
    timestamp         = Sys.time(),
    metaXpress_version = utils::packageVersion("metaXpress")
  )
}
