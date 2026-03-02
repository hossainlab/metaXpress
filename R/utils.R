# ============================================================================
# Internal utility functions shared across modules
# ============================================================================
# These functions are NOT exported. Use the package:: prefix if calling
# from vignettes or tests.

#' Validate that an object is a non-empty list of metaXpressStudy objects
#' @noRd
.validate_study_list <- function(studies, arg_name = "studies") {
  if (!is.list(studies) || length(studies) == 0)
    stop("'", arg_name, "' must be a non-empty list of metaXpressStudy objects")
  not_study <- !vapply(studies, is, logical(1), "metaXpressStudy")
  if (any(not_study))
    stop("'", arg_name, "' element(s) ",
         paste(which(not_study), collapse = ", "),
         " are not metaXpressStudy objects")
  invisible(TRUE)
}

#' Validate that an object is a list of DE result data.frames
#' @noRd
.validate_de_list <- function(de_results, arg_name = "de_results") {
  if (!is.list(de_results) || length(de_results) == 0)
    stop("'", arg_name, "' must be a non-empty list of data.frames")

  required_cols <- c("gene_id", "log2FC", "pvalue", "padj")
  for (i in seq_along(de_results)) {
    if (!is.data.frame(de_results[[i]]))
      stop("'", arg_name, "[[", i, "]]' must be a data.frame")
    missing_cols <- setdiff(required_cols, colnames(de_results[[i]]))
    if (length(missing_cols) > 0)
      stop(sprintf("'%s[[%d]]' is missing required columns: %s",
                   arg_name, i, paste(missing_cols, collapse = ", ")))
  }
  invisible(TRUE)
}

#' Extract DE results from a list that may contain either studies or data.frames
#' @noRd
.extract_de_results <- function(x) {
  if (is.list(x) && length(x) > 0 && is(x[[1]], "metaXpressStudy"))
    return(lapply(x, function(s) s@de_result))
  x
}

#' Safe log transform: log1p with a floor at 0
#' @noRd
.safe_log1p <- function(x) {
  x[x < 0] <- 0
  log1p(x)
}

#' Suppress messages and warnings during package loading
#' @noRd
.quietly <- function(expr) {
  suppressMessages(suppressWarnings(expr))
}

#' Check that a package is installed, with an informative error
#' @noRd
.require_package <- function(pkg, purpose = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- sprintf("Package '%s' is required", pkg)
    if (!is.null(purpose)) msg <- paste0(msg, " for ", purpose)
    msg <- paste0(msg, ". Install with: BiocManager::install('", pkg, "')")
    stop(msg)
  }
  invisible(TRUE)
}

#' Build a gene x studies matrix from a list of DE data.frames
#' @noRd
.build_gene_matrix <- function(de_results, value_col = "log2FC") {
  all_genes <- Reduce(union, lapply(de_results, function(d) d$gene_id))
  k         <- length(de_results)

  mat <- matrix(NA_real_, nrow = length(all_genes), ncol = k,
                 dimnames = list(all_genes, names(de_results)))

  for (i in seq_len(k)) {
    d   <- de_results[[i]]
    idx <- match(d$gene_id, all_genes)
    mat[idx, i] <- d[[value_col]]
  }
  mat
}
