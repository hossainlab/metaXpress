#' @include AllClasses.R
NULL

# ============================================================================
# Module 8: Report Generation & Export
# ============================================================================

#' Generate a reproducible analysis report
#'
#' Renders an HTML (or PDF) report with methods text, study table,
#' fold-change concordance, volcano plot, top-gene table with confidence and
#' prediction intervals, heterogeneity, forest plots, leave-one-out
#' robustness, pathway results (if present) and session information.
#'
#' @param meta_result A \code{\linkS4class{metaXpressResult}} object from
#'   \code{\link{mx_meta}}.
#' @param de_results The named list of per-study DE \code{data.frame}s (or
#'   \code{\linkS4class{metaXpressStudy}} objects) used for the meta-analysis.
#' @param studies Optional list of \code{\linkS4class{metaXpressStudy}}
#'   objects (currently unused; kept for API compatibility).
#' @param format \code{"html"} (default) or \code{"pdf"} (requires LaTeX).
#' @param output_dir Output directory. Default: \code{"."}.
#' @param output_file File name without extension. Default:
#'   \code{"metaXpress_report"}.
#' @param padj_threshold,lfc_threshold Significance thresholds used in the
#'   report. Defaults: \code{0.05} and \code{0}.
#'
#' @return The path to the generated report (invisibly).
#'
#' @examples
#' \donttest{
#' if (rmarkdown::pandoc_available()) {
#'   ad <- readRDS(system.file("extdata", "ad_brain_de.rds",
#'                             package = "metaXpress"))
#'   res <- mx_meta(ad$de_results)
#'   mx_report(res, ad$de_results, output_dir = tempdir())
#' }
#' }
#'
#' @seealso \code{\link{mx_export}}, \code{\link{mx_session_info}}
#'
#' @importFrom methods is
#' @export
mx_report <- function(meta_result, de_results, studies = NULL,
                      format         = c("html", "pdf"),
                      output_dir     = ".",
                      output_file    = "metaXpress_report",
                      padj_threshold = 0.05,
                      lfc_threshold  = 0) {
  if (!is(meta_result, "metaXpressResult"))
    stop("'meta_result' must be a metaXpressResult object")
  if (is.list(de_results) && length(de_results) > 0 &&
      is(de_results[[1]], "metaXpressStudy")) {
    if (is.null(studies)) studies <- de_results
    de_results <- .extract_de_results(de_results)
  }
  .validate_de_list(de_results)
  if (is.null(names(de_results)))
    names(de_results) <- paste0("study_", seq_along(de_results))
  format <- match.arg(format)
  .require_package("rmarkdown", "report generation")

  template <- system.file("report", "metaXpress_report.Rmd",
                          package = "metaXpress")
  if (!nzchar(template))
    stop("Report template not found. Reinstall metaXpress.")

  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)
  output_dir <- normalizePath(output_dir)

  work <- tempfile("mx_report_")
  dir.create(work)
  file.copy(template, work)

  output_format <- if (format == "html")
    rmarkdown::html_document(toc = TRUE, toc_float = TRUE)
  else
    rmarkdown::pdf_document(toc = TRUE)

  out_path <- rmarkdown::render(
    input         = file.path(work, basename(template)),
    output_format = output_format,
    output_file   = paste0(output_file, ".", format),
    output_dir    = output_dir,
    params        = list(meta_result = meta_result,
                         de_results = de_results,
                         studies = studies,
                         padj_threshold = padj_threshold,
                         lfc_threshold = lfc_threshold),
    envir = new.env(parent = globalenv()),
    quiet = TRUE
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
