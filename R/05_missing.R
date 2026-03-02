#' @include AllClasses.R
NULL

# ============================================================================
# Module 5: Missing Gene Handling
# ============================================================================

#' Summarise gene coverage across studies
#'
#' Returns a matrix showing which genes are present (1) or absent (0) in each
#' study's DE result, enabling informed decisions about how to handle partial
#' overlap before meta-analysis.
#'
#' @param de_results A named list of \code{data.frame} objects from
#'   \code{\link{mx_de_all}}, or a list of
#'   \code{\linkS4class{metaXpressStudy}} objects.
#'
#' @return A binary matrix (genes x studies) where 1 indicates the gene is
#'   present and 0 indicates it is absent. Attribute \code{"coverage_pct"}
#'   contains the percentage of studies covering each gene.
#'
#' @examples
#' \dontrun{
#'   cov_mat <- mx_missing_summary(de_results)
#'   hist(attr(cov_mat, "coverage_pct"), main = "Gene coverage across studies")
#' }
#'
#' @seealso \code{\link{mx_impute}}, \code{\link{mx_filter_coverage}}
#'
#' @export
mx_missing_summary <- function(de_results) {
  if (is.list(de_results) && all(vapply(de_results, is, logical(1),
                                         "metaXpressStudy")))
    de_results <- lapply(de_results, function(s) s@de_result)

  if (!is.list(de_results))
    stop("'de_results' must be a list of data.frames from mx_de_all()")

  all_genes <- Reduce(union, lapply(de_results, function(d) d$gene_id))
  k         <- length(de_results)

  cov_mat <- matrix(0L, nrow = length(all_genes), ncol = k,
                     dimnames = list(all_genes, names(de_results)))

  for (i in seq_len(k)) {
    present <- de_results[[i]]$gene_id
    cov_mat[present, i] <- 1L
  }

  attr(cov_mat, "coverage_pct") <- rowMeans(cov_mat) * 100
  cov_mat
}

#' Impute missing gene statistics across studies
#'
#' Handles genes not measured in all studies by applying one of four
#' imputation strategies before meta-analysis.
#'
#' @param de_results A named list of \code{data.frame} objects from
#'   \code{\link{mx_de_all}}.
#' @param method Character scalar. Imputation strategy. One of:
#'   \describe{
#'     \item{\code{"exclude"}}{Exclude genes not present in all studies
#'       (most conservative; default).}
#'     \item{\code{"mean"}}{Impute missing log2FC with the mean across
#'       available studies; impute p-values as 1.}
#'     \item{\code{"knn"}}{K-nearest neighbours imputation based on
#'       gene expression profiles.}
#'     \item{\code{"weighted"}}{Weighted imputation proportional to study
#'       sample size.}
#'   }
#'
#' @return The input list with missing values filled according to the chosen
#'   strategy.
#'
#' @references
#' Villatoro-García et al. (2022) Mathematics.
#' \doi{10.3390/math10183376}
#'
#' @examples
#' \dontrun{
#'   de_results <- mx_impute(de_results, method = "mean")
#' }
#'
#' @seealso \code{\link{mx_missing_summary}}, \code{\link{mx_filter_coverage}}
#'
#' @export
mx_impute <- function(de_results,
                       method = c("exclude", "mean", "knn", "weighted")) {
  if (is.list(de_results) && all(vapply(de_results, is, logical(1),
                                         "metaXpressStudy")))
    de_results <- lapply(de_results, function(s) s@de_result)

  if (!is.list(de_results))
    stop("'de_results' must be a list of data.frames from mx_de_all()")

  method <- match.arg(method)

  switch(method,
    exclude  = .impute_exclude(de_results),
    mean     = .impute_mean(de_results),
    knn      = stop("KNN imputation is not yet implemented."),
    weighted = stop("Weighted imputation is not yet implemented.")
  )
}

#' Filter genes by minimum study coverage
#'
#' Retains only genes present in at least \code{min_studies} studies.
#' This is the simplest approach to handling missing genes and is applied
#' automatically inside \code{\link{mx_meta}} via the \code{min_studies}
#' argument.
#'
#' @param de_results A named list of \code{data.frame} objects from
#'   \code{\link{mx_de_all}}.
#' @param min_studies Integer. Minimum number of studies a gene must appear in.
#'   Default: \code{2}.
#'
#' @return The input list with each \code{data.frame} restricted to genes
#'   meeting the coverage threshold.
#'
#' @examples
#' \dontrun{
#'   de_results <- mx_filter_coverage(de_results, min_studies = 3)
#' }
#'
#' @seealso \code{\link{mx_missing_summary}}, \code{\link{mx_impute}}
#'
#' @export
mx_filter_coverage <- function(de_results, min_studies = 2) {
  if (!is.list(de_results))
    stop("'de_results' must be a list of data.frames from mx_de_all()")
  if (!is.numeric(min_studies) || min_studies < 1)
    stop("'min_studies' must be a positive integer")

  cov_mat   <- mx_missing_summary(de_results)
  n_covered <- rowSums(cov_mat)
  keep_genes <- names(n_covered)[n_covered >= min_studies]

  n_kept <- length(keep_genes)
  n_total <- nrow(cov_mat)
  message(sprintf("Keeping %d / %d genes present in >= %d studies (%.1f%%).",
                   n_kept, n_total, min_studies,
                   n_kept / n_total * 100))

  lapply(de_results, function(d) {
    d[d$gene_id %in% keep_genes, , drop = FALSE]
  })
}

# ============================================================================
# Internal helpers
# ============================================================================

.impute_exclude <- function(de_results) {
  all_gene_sets <- lapply(de_results, function(d) d$gene_id)
  common_genes  <- Reduce(intersect, all_gene_sets)

  n_all    <- length(Reduce(union, all_gene_sets))
  n_common <- length(common_genes)
  message(sprintf("Excluding strategy: keeping %d common genes (%.1f%% of all %d).",
                   n_common, n_common / n_all * 100, n_all))

  lapply(de_results, function(d) {
    d[d$gene_id %in% common_genes, , drop = FALSE]
  })
}

.impute_mean <- function(de_results) {
  all_genes <- Reduce(union, lapply(de_results, function(d) d$gene_id))
  k         <- length(de_results)

  lfc_mat  <- matrix(NA_real_, nrow = length(all_genes), ncol = k,
                      dimnames = list(all_genes, NULL))
  pval_mat <- lfc_mat

  for (i in seq_len(k)) {
    d <- de_results[[i]]
    idx <- match(d$gene_id, all_genes)
    lfc_mat[idx, i]  <- d$log2FC
    pval_mat[idx, i] <- d$pvalue
  }

  mean_lfc <- rowMeans(lfc_mat, na.rm = TRUE)
  for (i in seq_len(k)) {
    na_rows <- is.na(lfc_mat[, i])
    lfc_mat[na_rows, i]  <- mean_lfc[na_rows]
    pval_mat[na_rows, i] <- 1.0
  }

  lapply(seq_len(k), function(i) {
    data.frame(
      gene_id  = all_genes,
      log2FC   = lfc_mat[, i],
      pvalue   = pval_mat[, i],
      padj     = p.adjust(pval_mat[, i], method = "BH"),
      baseMean = NA_real_,
      stringsAsFactors = FALSE
    )
  })
}
