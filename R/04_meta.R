#' @include AllClasses.R
NULL

# ============================================================================
# Module 4: Meta-Analysis Statistics
# ============================================================================

#' Run meta-analysis across per-study DE results
#'
#' Integrates per-study differential expression results from
#' \code{\link{mx_de_all}} into a single meta-analysis result using one of
#' six statistical methods.
#'
#' @param de_results A named list of \code{data.frame} objects from
#'   \code{\link{mx_de_all}} (i.e., the \code{de_result} slots extracted from
#'   each study), or a list of \code{\linkS4class{metaXpressStudy}} objects.
#' @param method Character scalar. Meta-analysis method. One of:
#'   \describe{
#'     \item{\code{"fisher"}}{Fisher's combined p-value (Rau et al. 2013)}
#'     \item{\code{"stouffer"}}{Stouffer's Z-score weighted by sample size}
#'     \item{\code{"inverse_normal"}}{Fused inverse-normal (Prasad & Li 2021)}
#'     \item{\code{"fixed_effects"}}{Fixed effects, inverse-variance weighting}
#'     \item{\code{"random_effects"}}{DerSimonian-Laird random effects (default)}
#'     \item{\code{"awmeta"}}{Adaptive weighting (Hu et al. 2025)}
#'   }
#' @param min_studies Integer. Minimum number of studies a gene must appear in
#'   to be included. Default: \code{2}.
#' @param alpha Numeric. Significance level for the BH FDR correction.
#'   Default: \code{0.05}.
#' @param n_samples Integer vector. Sample sizes (total n per study) used for
#'   Stouffer weighting. If \code{NULL}, uniform weights are used.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{\linkS4class{metaXpressResult}} object.
#'
#' @details
#' \strong{CRITICAL:} Use \code{meta_padj} (BH-adjusted meta p-value) for
#' significance filtering, not \code{meta_pvalue}. Default significance
#' threshold: \code{meta_padj <= 0.05} AND \code{|meta_log2FC| >= 1}.
#'
#' @references
#' Fisher, R.A. (1932) \emph{Statistical Methods for Research Workers}.
#' 4th edn. Edinburgh: Oliver & Boyd.
#'
#' Stouffer, S.A. et al. (1949) \emph{The American Soldier: Adjustment During
#' Army Life}. Princeton University Press.
#'
#' DerSimonian, R. & Laird, N. (1986) Meta-analysis in clinical trials.
#' \emph{Controlled Clinical Trials}, \strong{7}(3), 177--188.
#' \doi{10.1016/0197-2456(86)90046-2}
#'
#' Cochran, W.G. (1954) The combination of estimates from different
#' experiments. \emph{Biometrics}, \strong{10}(1), 101--129.
#' \doi{10.2307/3001666}
#'
#' Higgins, J.P.T. & Thompson, S.G. (2002) Quantifying heterogeneity in a
#' meta-analysis. \emph{Statistics in Medicine}, \strong{21}(11), 1539--1558.
#' \doi{10.1002/sim.1186}
#'
#' Rau, A., Marot, G. & Jaffrézic, F. (2014) Differential meta-analysis of
#' RNA-seq data from multiple studies. \emph{BMC Bioinformatics},
#' \strong{15}, 91. \doi{10.1186/1471-2105-15-91}
#'
#' Prasad, A. & Li, Q. (2022) Fused inverse-normal method for integrated
#' differential expression analysis of RNA-seq data. \emph{BMC
#' Bioinformatics}, \strong{23}, 371. \doi{10.1186/s12859-022-04859-9}
#'
#' Keel, B.N. & Lindholm-Perry, A.K. (2022) Recent developments and future
#' directions in meta-analysis of differential gene expression.
#' \emph{Frontiers in Genetics}, \strong{13}, 983043.
#' \doi{10.3389/fgene.2022.983043}
#'
#' Hu, J. et al. (2025) AWmeta: adaptive weighting meta-analysis for
#' combining heterogeneous genomic studies. \emph{bioRxiv}.
#' \doi{10.1101/2025.05.06.650408}
#'
#' @examples
#' \dontrun{
#'   de_list <- lapply(studies, function(s) s@de_result)
#'   result  <- mx_meta(de_list, method = "random_effects")
#'   head(result@meta_table[order(result@meta_table$meta_padj), ])
#' }
#'
#' @seealso \code{\link{mx_heterogeneity}}, \code{\link{mx_sensitivity}},
#'   \code{\link{mx_forest}}
#'
#' @importFrom methods new
#' @importFrom stats p.adjust pchisq pnorm qnorm
#' @export
mx_meta <- function(de_results,
                     method      = c("random_effects", "fisher", "stouffer",
                                     "inverse_normal", "fixed_effects",
                                     "awmeta"),
                     min_studies = 2,
                     alpha       = 0.05,
                     n_samples   = NULL,
                     ...) {
  # --- Input validation ---
  if (is.list(de_results) && all(vapply(de_results, is, logical(1),
                                         "metaXpressStudy")))
    de_results <- lapply(de_results, function(s) s@de_result)

  if (!is.list(de_results))
    stop("'de_results' must be a list of data.frames from mx_de_all()")

  method <- match.arg(method)

  if (!is.numeric(min_studies) || length(min_studies) != 1 || min_studies < 2)
    stop("'min_studies' must be a single integer >= 2")

  if (method %in% c("random_effects", "awmeta") && length(de_results) < 3)
    warning("'", method, "' is unreliable with fewer than 3 studies. ",
            "Consider method = 'fisher' or method = 'stouffer'.")

  required_cols <- c("gene_id", "log2FC", "pvalue", "padj")
  for (i in seq_along(de_results)) {
    missing_cols <- setdiff(required_cols, colnames(de_results[[i]]))
    if (length(missing_cols) > 0)
      stop(sprintf("de_results[[%d]] is missing columns: %s",
                   i, paste(missing_cols, collapse = ", ")))
  }

  message("Running meta-analysis (", method, ") on ",
          length(de_results), " studies...")

  # --- Build aligned gene matrix ---
  all_genes <- Reduce(union, lapply(de_results, function(d) d$gene_id))
  k         <- length(de_results)

  lfc_mat   <- matrix(NA_real_, nrow = length(all_genes), ncol = k,
                       dimnames = list(all_genes, names(de_results)))
  pval_mat  <- lfc_mat
  se_mat    <- lfc_mat

  for (i in seq_len(k)) {
    d <- de_results[[i]]
    idx <- match(d$gene_id, all_genes)
    lfc_mat[idx, i]  <- d$log2FC
    pval_mat[idx, i] <- d$pvalue
    if ("lfcSE" %in% colnames(d))
      se_mat[idx, i] <- d$lfcSE
  }

  # Filter by coverage
  n_studies_per_gene <- rowSums(!is.na(pval_mat))
  keep <- n_studies_per_gene >= min_studies
  if (sum(keep) == 0)
    stop("No genes found in >= ", min_studies, " studies. ",
         "Lower 'min_studies' or check gene ID alignment.")

  lfc_mat  <- lfc_mat[keep, , drop = FALSE]
  pval_mat <- pval_mat[keep, , drop = FALSE]
  se_mat   <- se_mat[keep, , drop = FALSE]

  # --- Compute heterogeneity ---
  het <- t(apply(cbind(lfc_mat, se_mat), 1, function(row) {
    lfc <- row[seq_len(k)]
    se  <- row[(k + 1):(2 * k)]
    res <- .compute_heterogeneity(lfc, se)
    c(Q = res$Q, df = res$df, I_sq = res$I_sq, tau_sq = res$tau_sq,
      p_het = res$p_heterogeneity)
  }))

  # --- Compute per-sample n for weighting ---
  if (is.null(n_samples))
    n_samples <- rep(10L, k)

  # --- Run chosen method ---
  meta_stats <- switch(method,
    fisher         = .meta_fisher(pval_mat, lfc_mat),
    stouffer       = .meta_stouffer(pval_mat, lfc_mat, n_samples),
    inverse_normal = .meta_inverse_normal(pval_mat, lfc_mat, n_samples),
    fixed_effects  = .meta_fixed_effects(lfc_mat, se_mat),
    random_effects = .meta_random_effects(lfc_mat, se_mat, het[, "tau_sq"]),
    awmeta         = .meta_awmeta(pval_mat, lfc_mat, se_mat, het[, "I_sq"])
  )

  # --- Assemble result table ---
  direction_consistency <- rowMeans(
    sign(lfc_mat) == sign(meta_stats$meta_log2FC), na.rm = TRUE)

  meta_table <- data.frame(
    gene_id               = rownames(lfc_mat),
    meta_log2FC           = meta_stats$meta_log2FC,
    meta_pvalue           = meta_stats$meta_pvalue,
    meta_padj             = p.adjust(meta_stats$meta_pvalue, method = "BH"),
    i_squared             = het[, "I_sq"],
    q_stat                = het[, "Q"],
    n_studies             = as.integer(rowSums(!is.na(pval_mat))),
    direction_consistency = direction_consistency,
    stringsAsFactors      = FALSE,
    row.names             = NULL
  )

  heterogeneity_df <- data.frame(
    gene_id         = rownames(lfc_mat),
    Q               = het[, "Q"],
    df              = het[, "df"],
    I_sq            = het[, "I_sq"],
    tau_sq          = het[, "tau_sq"],
    p_heterogeneity = het[, "p_het"],
    stringsAsFactors = FALSE,
    row.names        = NULL
  )

  new("metaXpressResult",
      meta_table     = meta_table,
      method         = method,
      n_studies      = as.integer(k),
      heterogeneity  = heterogeneity_df,
      pathway_result = data.frame())
}

#' Compute per-gene heterogeneity statistics
#'
#' Returns Cochran's Q, degrees of freedom, I-squared, tau-squared, and the
#' p-value for the Q test for each gene across studies.
#'
#' @param de_results A named list of \code{data.frame} objects from
#'   \code{\link{mx_de_all}}, or a list of
#'   \code{\linkS4class{metaXpressStudy}} objects.
#'
#' @return A \code{data.frame} with columns \code{gene_id}, \code{Q},
#'   \code{df}, \code{I_sq}, \code{tau_sq}, \code{p_heterogeneity}.
#'
#' @references
#' Cochran, W.G. (1954) The combination of estimates from different
#' experiments. \emph{Biometrics}, \strong{10}(1), 101--129.
#' \doi{10.2307/3001666}
#'
#' Higgins, J.P.T. & Thompson, S.G. (2002) Quantifying heterogeneity in a
#' meta-analysis. \emph{Statistics in Medicine}, \strong{21}(11), 1539--1558.
#' \doi{10.1002/sim.1186}
#'
#' @examples
#' \dontrun{
#'   het <- mx_heterogeneity(de_results)
#'   hist(het$I_sq, main = "I-squared distribution")
#' }
#'
#' @seealso \code{\link{mx_meta}}, \code{\link{mx_heterogeneity_plot}}
#'
#' @export
mx_heterogeneity <- function(de_results) {
  if (is.list(de_results) && all(vapply(de_results, is, logical(1),
                                         "metaXpressStudy")))
    de_results <- lapply(de_results, function(s) s@de_result)

  if (!is.list(de_results))
    stop("'de_results' must be a list of data.frames from mx_de_all()")

  all_genes <- Reduce(union, lapply(de_results, function(d) d$gene_id))
  k         <- length(de_results)
  lfc_mat   <- matrix(NA_real_, nrow = length(all_genes), ncol = k,
                       dimnames = list(all_genes, names(de_results)))
  se_mat    <- lfc_mat

  for (i in seq_len(k)) {
    d <- de_results[[i]]
    idx <- match(d$gene_id, all_genes)
    lfc_mat[idx, i] <- d$log2FC
    if ("lfcSE" %in% colnames(d))
      se_mat[idx, i] <- d$lfcSE
  }

  het_list <- apply(cbind(lfc_mat, se_mat), 1, function(row) {
    lfc <- row[seq_len(k)]
    se  <- row[(k + 1):(2 * k)]
    .compute_heterogeneity(lfc, se)
  })

  data.frame(
    gene_id         = all_genes,
    Q               = vapply(het_list, `[[`, numeric(1), "Q"),
    df              = vapply(het_list, `[[`, numeric(1), "df"),
    I_sq            = vapply(het_list, `[[`, numeric(1), "I_sq"),
    tau_sq          = vapply(het_list, `[[`, numeric(1), "tau_sq"),
    p_heterogeneity = vapply(het_list, `[[`, numeric(1), "p_heterogeneity"),
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
}

#' Leave-one-out sensitivity analysis
#'
#' Reruns the meta-analysis iteratively, each time excluding one study, to
#' assess whether any single study drives the overall result. This is the
#' standard sensitivity analysis recommended by the Cochrane Handbook
#' (Higgins et al. 2019, Section 10.14).
#'
#' @param de_results A named list of \code{data.frame} objects from
#'   \code{\link{mx_de_all}}.
#' @param method Character scalar. Meta-analysis method. Default:
#'   \code{"random_effects"}.
#'
#' @return A list of \code{\linkS4class{metaXpressResult}} objects, one per
#'   leave-one-out iteration. Names correspond to the excluded study.
#'
#' @references
#' Higgins, J.P.T. et al. (eds.) (2019) \emph{Cochrane Handbook for
#' Systematic Reviews of Interventions}. 2nd edn. Chichester: John Wiley
#' & Sons. \doi{10.1002/9781119536604}
#'
#' @examples
#' \dontrun{
#'   loo <- mx_sensitivity(de_results, method = "random_effects")
#'   # Compare n significant genes across LOO runs
#'   vapply(loo, function(r) sum(r@meta_table$meta_padj <= 0.05, na.rm = TRUE),
#'          integer(1))
#' }
#'
#' @seealso \code{\link{mx_meta}}
#'
#' @export
mx_sensitivity <- function(de_results,
                             method = c("random_effects", "fisher", "stouffer",
                                        "inverse_normal", "fixed_effects",
                                        "awmeta")) {
  if (is.list(de_results) && all(vapply(de_results, is, logical(1),
                                         "metaXpressStudy")))
    de_results <- lapply(de_results, function(s) s@de_result)

  method <- match.arg(method)

  if (length(de_results) < 3)
    stop("Sensitivity analysis requires at least 3 studies")

  message("Running leave-one-out sensitivity analysis...")

  study_names <- names(de_results)
  if (is.null(study_names))
    study_names <- paste0("study_", seq_along(de_results))

  loo_results <- lapply(seq_along(de_results), function(i) {
    message("  Excluding: ", study_names[i])
    mx_meta(de_results[-i], method = method, min_studies = 2)
  })

  names(loo_results) <- study_names
  loo_results
}

# ============================================================================
# Internal statistical helpers
# ============================================================================

.fisher_combine <- function(pvalues) {
  pvalues <- pvalues[!is.na(pvalues)]
  k <- length(pvalues)
  if (k < 2) return(NA_real_)
  T_stat <- -2 * sum(log(pvalues))
  pchisq(T_stat, df = 2 * k, lower.tail = FALSE)
}

.stouffer_combine <- function(pvalues, log2fc, n_samples) {
  valid   <- !is.na(pvalues) & !is.na(log2fc) & !is.na(n_samples)
  pvalues <- pvalues[valid]; log2fc <- log2fc[valid]; n <- n_samples[valid]
  if (length(pvalues) < 2) return(list(meta_p = NA_real_, meta_z = NA_real_))
  z_scores <- qnorm(1 - pvalues / 2) * sign(log2fc)
  weights  <- sqrt(n)
  z_comb   <- sum(weights * z_scores) / sqrt(sum(weights^2))
  list(meta_p = 2 * pnorm(-abs(z_comb)), meta_z = z_comb)
}

.inverse_normal_combine <- function(pvalues, log2fc, n_samples) {
  valid   <- !is.na(pvalues) & !is.na(log2fc) & !is.na(n_samples)
  pvalues <- pvalues[valid]; log2fc <- log2fc[valid]; n <- n_samples[valid]
  if (length(pvalues) < 2) return(NA_real_)
  z_scores <- qnorm(1 - pvalues) * sign(log2fc)
  weights  <- sqrt(n) * abs(log2fc)
  weights[weights == 0] <- .Machine$double.eps
  z_fused  <- sum(weights * z_scores) / sqrt(sum(weights^2))
  2 * pnorm(-abs(z_fused))
}

.compute_heterogeneity <- function(log2fc_vec, se_vec) {
  valid <- !is.na(log2fc_vec) & !is.na(se_vec) & se_vec > 0
  lfc   <- log2fc_vec[valid]; se <- se_vec[valid]
  k     <- length(lfc)
  if (k < 2) return(list(Q = NA_real_, df = NA_real_, I_sq = NA_real_,
                          tau_sq = NA_real_, p_heterogeneity = NA_real_))

  w        <- 1 / se^2
  theta_FE <- sum(w * lfc) / sum(w)
  Q        <- sum(w * (lfc - theta_FE)^2)
  df       <- k - 1
  I_sq     <- max(0, (Q - df) / Q) * 100
  C        <- sum(w) - sum(w^2) / sum(w)
  tau_sq   <- max(0, (Q - df) / C)
  p_het    <- pchisq(Q, df = df, lower.tail = FALSE)

  list(Q = Q, df = df, I_sq = I_sq, tau_sq = tau_sq,
       p_heterogeneity = p_het)
}

.meta_fisher <- function(pval_mat, lfc_mat) {
  meta_pvalue <- apply(pval_mat, 1, .fisher_combine)
  # Fisher combines p-values only; use unweighted mean log2FC as the

  # summary effect direction (not a formal estimate)
  meta_log2FC <- rowMeans(lfc_mat, na.rm = TRUE)
  list(meta_pvalue = meta_pvalue, meta_log2FC = meta_log2FC)
}

.meta_stouffer <- function(pval_mat, lfc_mat, n_samples) {
  meta_pvalue <- numeric(nrow(pval_mat))
  meta_z      <- numeric(nrow(pval_mat))

  for (g in seq_len(nrow(pval_mat))) {
    res            <- .stouffer_combine(pval_mat[g, ], lfc_mat[g, ], n_samples)
    meta_pvalue[g] <- res$meta_p
    meta_z[g]      <- res$meta_z
  }

  meta_log2FC <- rowMeans(lfc_mat, na.rm = TRUE)
  list(meta_pvalue = meta_pvalue, meta_log2FC = meta_log2FC)
}

.meta_inverse_normal <- function(pval_mat, lfc_mat, n_samples) {
  meta_pvalue <- vapply(seq_len(nrow(pval_mat)), function(g) {
    .inverse_normal_combine(pval_mat[g, ], lfc_mat[g, ], n_samples)
  }, numeric(1))
  meta_log2FC <- rowMeans(lfc_mat, na.rm = TRUE)
  list(meta_pvalue = meta_pvalue, meta_log2FC = meta_log2FC)
}

.meta_fixed_effects <- function(lfc_mat, se_mat) {
  apply_fe <- function(g) {
    lfc <- lfc_mat[g, ]; se <- se_mat[g, ]
    valid <- !is.na(lfc) & !is.na(se) & se > 0
    lfc <- lfc[valid]; se <- se[valid]
    if (length(lfc) < 2) return(c(log2FC = NA_real_, pvalue = NA_real_))
    w        <- 1 / se^2
    theta_FE <- sum(w * lfc) / sum(w)
    SE_FE    <- 1 / sqrt(sum(w))
    Z_FE     <- theta_FE / SE_FE
    c(log2FC = theta_FE, pvalue = 2 * pnorm(-abs(Z_FE)))
  }

  res_mat     <- t(vapply(seq_len(nrow(lfc_mat)), apply_fe, numeric(2)))
  list(meta_pvalue = res_mat[, "pvalue"], meta_log2FC = res_mat[, "log2FC"])
}

.meta_random_effects <- function(lfc_mat, se_mat, tau_sq_vec) {
  apply_re <- function(g) {
    lfc <- lfc_mat[g, ]; se <- se_mat[g, ]
    tau_sq <- tau_sq_vec[g]
    valid  <- !is.na(lfc) & !is.na(se) & se > 0
    lfc <- lfc[valid]; se <- se[valid]
    if (length(lfc) < 2) return(c(log2FC = NA_real_, pvalue = NA_real_))
    if (is.na(tau_sq)) tau_sq <- 0
    w_re     <- 1 / (se^2 + tau_sq)
    theta_RE <- sum(w_re * lfc) / sum(w_re)
    SE_RE    <- 1 / sqrt(sum(w_re))
    Z_RE     <- theta_RE / SE_RE
    c(log2FC = theta_RE, pvalue = 2 * pnorm(-abs(Z_RE)))
  }

  res_mat <- t(vapply(seq_len(nrow(lfc_mat)), apply_re, numeric(2)))
  list(meta_pvalue = res_mat[, "pvalue"], meta_log2FC = res_mat[, "log2FC"])
}

.meta_awmeta <- function(pval_mat, lfc_mat, se_mat, i_sq_vec) {
  apply_aw <- function(g) {
    i_sq <- i_sq_vec[g]
    if (is.na(i_sq)) i_sq <- 50
    lambda <- 1 - i_sq / 100

    # Effect size Z
    lfc <- lfc_mat[g, ]; se <- se_mat[g, ]
    valid <- !is.na(lfc) & !is.na(se) & se > 0
    if (sum(valid) < 2) return(c(log2FC = NA_real_, pvalue = NA_real_))
    w_re     <- 1 / se[valid]^2
    theta_RE <- sum(w_re * lfc[valid]) / sum(w_re)
    SE_RE    <- 1 / sqrt(sum(w_re))
    Z_effect <- theta_RE / SE_RE

    # P-value Z (Stouffer)
    pv <- pval_mat[g, ]; lf <- lfc_mat[g, ]
    pv_valid <- !is.na(pv) & !is.na(lf)
    pv <- pv[pv_valid]; lf <- lf[pv_valid]
    if (length(pv) < 2) return(c(log2FC = NA_real_, pvalue = NA_real_))
    z_pv <- qnorm(1 - pv / 2) * sign(lf)
    Z_pvalue <- mean(z_pv)

    Z_adaptive <- lambda * Z_effect + (1 - lambda) * Z_pvalue
    c(log2FC = theta_RE, pvalue = 2 * pnorm(-abs(Z_adaptive)))
  }

  res_mat <- t(vapply(seq_len(nrow(lfc_mat)), apply_aw, numeric(2)))
  list(meta_pvalue = res_mat[, "pvalue"], meta_log2FC = res_mat[, "log2FC"])
}
