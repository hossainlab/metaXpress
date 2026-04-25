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
#'       gene expression profiles (Hastie et al. 1999).}
#'     \item{\code{"weighted"}}{Weighted imputation proportional to study
#'       sample size.}
#'   }
#' @param weights Numeric vector. Study weights for \code{"weighted"} imputation.
#'   Default: \code{NULL} (uses equal weights, equivalent to \code{"mean"}).
#'
#' @return The input list with missing values filled according to the chosen
#'   strategy.
#'
#' @references
#' Hastie, T. et al. (1999) Imputing missing data for gene expression arrays.
#' Stanford University Statistics Department Technical report.
#'
#' Villatoro-García, J.A. et al. (2022) Missing gene expression data
#' imputation for gene-study meta-analysis. \emph{Mathematics},
#' \strong{10}(18), 3376. \doi{10.3390/math10183376}
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
                       method = c("exclude", "mean", "knn", "weighted"),
                       weights = NULL) {
  if (is.list(de_results) && all(vapply(de_results, is, logical(1),
                                         "metaXpressStudy")))
    de_results <- lapply(de_results, function(s) s@de_result)

  if (!is.list(de_results))
    stop("'de_results' must be a list of data.frames from mx_de_all()")

  method <- match.arg(method)

  switch(method,
    exclude  = .impute_exclude(de_results),
    mean     = .impute_mean(de_results),
    knn      = .impute_knn(de_results),
    weighted = .impute_weighted(de_results, weights)
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
  mats      <- .build_gene_matrices(de_results)
  all_genes <- mats$all_genes
  lfc_mat   <- mats$lfc_mat
  pval_mat  <- mats$pval_mat
  k         <- ncol(lfc_mat)

  mean_lfc <- rowMeans(lfc_mat, na.rm = TRUE)
  for (i in seq_len(k)) {
    na_rows <- is.na(lfc_mat[, i])
    lfc_mat[na_rows, i]  <- mean_lfc[na_rows]
    pval_mat[na_rows, i] <- 1.0
  }

  result <- lapply(seq_len(k), function(i) {
    data.frame(
      gene_id  = all_genes,
      log2FC   = lfc_mat[, i],
      pvalue   = pval_mat[, i],
      padj     = p.adjust(pval_mat[, i], method = "BH"),
      baseMean = NA_real_,
      stringsAsFactors = FALSE
    )
  })
  names(result) <- names(de_results)
  result
}

.impute_knn <- function(de_results) {
  .require_package("impute", "KNN imputation")
  
  mats      <- .build_gene_matrices(de_results)
  all_genes <- mats$all_genes
  lfc_mat   <- mats$lfc_mat
  pval_mat  <- mats$pval_mat
  k         <- ncol(lfc_mat)

  # impute.knn requires genes as rows and samples (studies) as columns
  # It throws error if any row has > 50% missing or any column has > 80% missing
  # Default is k=10, but we might have fewer than 10 studies
  knn_k <- min(10, k - 1)
  if (knn_k < 1) stop("Not enough studies for KNN imputation")
  
  # Filter out genes with > 50% missing or we'll get an error
  na_prop <- rowMeans(is.na(lfc_mat))
  valid_rows <- na_prop <= 0.5
  
  if (sum(!valid_rows) > 0) {
    message(sum(!valid_rows), " genes have > 50% missing values and will be excluded before KNN.")
  }
  
  lfc_valid <- lfc_mat[valid_rows, , drop = FALSE]
  pval_valid <- pval_mat[valid_rows, , drop = FALSE]
  valid_genes <- all_genes[valid_rows]
  
  imputed <- impute::impute.knn(as.matrix(lfc_valid), k = knn_k)
  lfc_imputed <- imputed$data
  
  for (i in seq_len(k)) {
    na_rows <- is.na(pval_valid[, i])
    pval_valid[na_rows, i] <- 1.0
  }

  result <- lapply(seq_len(k), function(i) {
    data.frame(
      gene_id  = valid_genes,
      log2FC   = lfc_imputed[, i],
      pvalue   = pval_valid[, i],
      padj     = p.adjust(pval_valid[, i], method = "BH"),
      baseMean = NA_real_,
      stringsAsFactors = FALSE
    )
  })
  names(result) <- names(de_results)
  result
}

.impute_weighted <- function(de_results, weights) {
  mats      <- .build_gene_matrices(de_results)
  all_genes <- mats$all_genes
  lfc_mat   <- mats$lfc_mat
  pval_mat  <- mats$pval_mat
  k         <- ncol(lfc_mat)

  if (is.null(weights)) {
    warning("'weights' not provided for weighted imputation; using uniform weights (mean).")
    weights <- rep(1, k)
  }
  if (length(weights) != k) {
    stop("Length of 'weights' must equal the number of studies")
  }
  
  # Normalize weights
  w_mat <- matrix(weights, nrow = nrow(lfc_mat), ncol = k, byrow = TRUE)
  w_mat[is.na(lfc_mat)] <- NA_real_
  
  row_w_sums <- rowSums(w_mat, na.rm = TRUE)
  weighted_lfc <- rowSums(lfc_mat * w_mat, na.rm = TRUE) / row_w_sums
  
  for (i in seq_len(k)) {
    na_rows <- is.na(lfc_mat[, i])
    lfc_mat[na_rows, i]  <- weighted_lfc[na_rows]
    pval_mat[na_rows, i] <- 1.0
  }

  result <- lapply(seq_len(k), function(i) {
    data.frame(
      gene_id  = all_genes,
      log2FC   = lfc_mat[, i],
      pvalue   = pval_mat[, i],
      padj     = p.adjust(pval_mat[, i], method = "BH"),
      baseMean = NA_real_,
      stringsAsFactors = FALSE
    )
  })
  names(result) <- names(de_results)
  result
}
