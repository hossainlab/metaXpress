#' Example metaXpress dataset
#'
#' A minimal simulated dataset for use in examples, vignettes, and tests.
#' Contains three \code{metaXpressStudy} objects representing hypothetical
#' case-vs-control bulk RNA-seq studies, each with 500 genes and 6 samples
#' (3 case, 3 control).
#'
#' @format A named list with three elements:
#' \describe{
#'   \item{\code{study1}}{A \code{\linkS4class{metaXpressStudy}} object
#'     (accession: \code{"SIM001"}, 500 genes, 6 samples)}
#'   \item{\code{study2}}{A \code{\linkS4class{metaXpressStudy}} object
#'     (accession: \code{"SIM002"}, 500 genes, 6 samples)}
#'   \item{\code{study3}}{A \code{\linkS4class{metaXpressStudy}} object
#'     (accession: \code{"SIM003"}, 480 genes, 8 samples — partial overlap)}
#' }
#'
#' @usage data(metaXpress_example)
#'
#' @examples
#' \dontrun{
#'   data(metaXpress_example)
#'   metaXpress_example$study1
#' }
#'
#' @name metaXpress_example
#' @docType data
NULL
