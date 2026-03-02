#' @include AllClasses.R
NULL

# ============================================================================
# show methods
# ============================================================================

#' Display a metaXpressStudy object
#'
#' @param object A \code{metaXpressStudy} object.
#' @return Invisibly returns \code{object}. Called for its side-effect of
#'   printing a summary to the console.
#' @exportMethod show
setMethod("show", "metaXpressStudy", function(object) {
  cat("metaXpressStudy\n")
  cat("  Accession :", object@accession, "\n")
  cat("  Organism  :", object@organism, "\n")
  cat("  Genes     :", nrow(object@counts), "\n")
  cat("  Samples   :", ncol(object@counts), "\n")
  qc <- if (is.na(object@qc_score)) "not run" else
    paste0(object@qc_score, " / 10")
  cat("  QC score  :", qc, "\n")
  cat("  DE run    :", nrow(object@de_result) > 0, "\n")
  invisible(object)
})

#' Display a metaXpressResult object
#'
#' @param object A \code{metaXpressResult} object.
#' @return Invisibly returns \code{object}. Called for its side-effect of
#'   printing a summary to the console.
#' @exportMethod show
setMethod("show", "metaXpressResult", function(object) {
  cat("metaXpressResult\n")
  cat("  Method    :", object@method, "\n")
  cat("  Studies   :", object@n_studies, "\n")
  cat("  Genes     :", nrow(object@meta_table), "\n")
  if (nrow(object@meta_table) > 0 &&
      "meta_padj" %in% colnames(object@meta_table)) {
    n_sig <- sum(object@meta_table$meta_padj <= 0.05, na.rm = TRUE)
    cat("  Sig genes :", n_sig, "(meta_padj <= 0.05)\n")
  }
  invisible(object)
})
