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

  # --- Build aligned gene matrices (single pass) ---
  mats     <- .build_gene_matrices(de_results)
  lfc_mat  <- mats$lfc_mat
  pval_mat <- mats$pval_mat
  se_mat   <- mats$se_mat
  k        <- ncol(lfc_mat)

  # Clamp p-values away from 0 and 1 to avoid -Inf/Inf in log/qnorm
  pval_mat[pval_mat == 0] <- .Machine$double.xmin
  pval_mat[pval_mat == 1] <- 1 - .Machine$double.eps

  # Filter by coverage
  n_studies_per_gene <- rowSums(!is.na(pval_mat))
  keep <- n_studies_per_gene >= min_studies
  if (sum(keep) == 0)
    stop("No genes found in >= ", min_studies, " studies. ",
         "Lower 'min_studies' or check gene ID alignment.")

  lfc_mat  <- lfc_mat[keep, , drop = FALSE]
  pval_mat <- pval_mat[keep, , drop = FALSE]
  se_mat   <- se_mat[keep, , drop = FALSE]

  # --- Compute heterogeneity (vectorized) ---
  het <- .compute_heterogeneity_mat(lfc_mat, se_mat)

  # --- Compute per-sample n for weighting ---
  if (is.null(n_samples)) {
    if (method %in% c("stouffer", "inverse_normal"))
      message("  Note: n_samples not provided; using uniform weights. ",
              "Pass n_samples for sample-size-weighted combination.")
    n_samples <- rep(1L, k)
  }

  # --- Run chosen method ---
  meta_stats <- switch(method,
    fisher         = .meta_fisher(pval_mat, lfc_mat),
    stouffer       = .meta_stouffer(pval_mat, lfc_mat, n_samples),
    inverse_normal = .meta_inverse_normal(pval_mat, lfc_mat, n_samples),
    fixed_effects  = .meta_fixed_effects(lfc_mat, se_mat),
    random_effects = .meta_random_effects(lfc_mat, se_mat, het$tau_sq),
    awmeta         = .meta_awmeta(pval_mat, lfc_mat, se_mat, het$I_sq)
  )

  # --- Assemble result table ---
  direction_consistency <- rowMeans(
    sign(lfc_mat) == sign(meta_stats$meta_log2FC), na.rm = TRUE)

  meta_table <- data.frame(
    gene_id               = rownames(lfc_mat),
    meta_log2FC           = meta_stats$meta_log2FC,
    meta_pvalue           = meta_stats$meta_pvalue,
    meta_padj             = p.adjust(meta_stats$meta_pvalue, method = "BH"),
    i_squared             = het$I_sq,
    q_stat                = het$Q,
    n_studies             = as.integer(rowSums(!is.na(pval_mat))),
    direction_consistency = direction_consistency,
    stringsAsFactors      = FALSE,
    row.names             = NULL
  )

  heterogeneity_df <- data.frame(
    gene_id         = rownames(lfc_mat),
    Q               = het$Q,
    df              = het$df,
    I_sq            = het$I_sq,
    tau_sq          = het$tau_sq,
    p_heterogeneity = het$p_heterogeneity,
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

  mats    <- .build_gene_matrices(de_results)
  het     <- .compute_heterogeneity_mat(mats$lfc_mat, mats$se_mat)

  data.frame(
    gene_id         = mats$all_genes,
    Q               = het$Q,
    df              = het$df,
    I_sq            = het$I_sq,
    tau_sq          = het$tau_sq,
    p_heterogeneity = het$p_heterogeneity,
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
# Internal statistical helpers — fully vectorized
# ============================================================================

#' Vectorized heterogeneity computation across all genes
#'
#' Operates on entire lfc_mat and se_mat matrices using matrix algebra.
#' @param lfc_mat Matrix of log2 fold changes (genes x studies).
#' @param se_mat Matrix of standard errors (genes x studies).
#' @return A list with vectors Q, df, I_sq, tau_sq, p_heterogeneity.
#' @noRd
.compute_heterogeneity_mat <- function(lfc_mat, se_mat) {
  n <- nrow(lfc_mat)

  # Mask invalid entries (NA or se <= 0)
  valid <- !is.na(lfc_mat) & !is.na(se_mat) & se_mat > 0
  lfc_v <- lfc_mat; lfc_v[!valid] <- NA_real_
  se_v  <- se_mat;  se_v[!valid]  <- NA_real_

  w <- 1 / se_v^2
  w[!valid] <- NA_real_

  # Per-gene count of valid studies
  k_vec <- rowSums(valid)

  # Fixed-effect estimate: theta_FE = sum(w*lfc) / sum(w)
  sum_w     <- rowSums(w, na.rm = TRUE)
  sum_w_lfc <- rowSums(w * lfc_v, na.rm = TRUE)
  theta_FE  <- sum_w_lfc / sum_w

  # Q = sum(w * (lfc - theta_FE)^2)
  Q <- rowSums(w * (lfc_v - theta_FE)^2, na.rm = TRUE)
  df <- k_vec - 1L

  # I_sq = max(0, (Q - df) / Q) * 100
  I_sq <- pmax(0, (Q - df) / Q) * 100

  # C = sum(w) - sum(w^2)/sum(w)
  sum_w2 <- rowSums(w^2, na.rm = TRUE)
  C_val  <- sum_w - sum_w2 / sum_w

  # tau_sq = max(0, (Q - df) / C)
  tau_sq <- pmax(0, (Q - df) / C_val)

  # p-value
  p_het <- pchisq(Q, df = df, lower.tail = FALSE)

  # Genes with < 2 valid studies get NA
  too_few <- k_vec < 2
  Q[too_few]     <- NA_real_
  df[too_few]    <- NA_real_
  I_sq[too_few]  <- NA_real_
  tau_sq[too_few] <- NA_real_
  p_het[too_few] <- NA_real_

  list(Q = Q, df = df, I_sq = I_sq, tau_sq = tau_sq,
       p_heterogeneity = p_het)
}

#' Scalar heterogeneity (kept for backward compatibility in edge cases)
#' @noRd
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

# --- Vectorized Fisher ---
.meta_fisher <- function(pval_mat, lfc_mat) {
  # T = -2 * sum(log(p)) per gene; df = 2k per gene
  log_p  <- log(pval_mat)
  log_p[is.na(pval_mat)] <- NA_real_
  T_stat <- -2 * rowSums(log_p, na.rm = TRUE)
  k_vec  <- rowSums(!is.na(pval_mat))

  meta_pvalue <- pchisq(T_stat, df = 2 * k_vec, lower.tail = FALSE)
  # Genes with < 2 studies -> NA
  meta_pvalue[k_vec < 2] <- NA_real_

  meta_log2FC <- rowMeans(lfc_mat, na.rm = TRUE)
  list(meta_pvalue = meta_pvalue, meta_log2FC = meta_log2FC)
}

# --- Vectorized Stouffer ---
.meta_stouffer <- function(pval_mat, lfc_mat, n_samples) {
  # z = qnorm(1 - p/2) * sign(lfc)
  valid <- !is.na(pval_mat) & !is.na(lfc_mat)
  z_mat <- matrix(NA_real_, nrow = nrow(pval_mat), ncol = ncol(pval_mat))
  z_mat[valid] <- qnorm(1 - pval_mat[valid] / 2) * sign(lfc_mat[valid])

  # Weight matrix: sqrt(n_samples), broadcast across rows
  w_mat <- matrix(sqrt(n_samples), nrow = nrow(pval_mat),
                   ncol = ncol(pval_mat), byrow = TRUE)
  w_mat[!valid] <- NA_real_

  # z_comb = sum(w * z) / sqrt(sum(w^2))
  z_comb <- rowSums(w_mat * z_mat, na.rm = TRUE) /
            sqrt(rowSums(w_mat^2, na.rm = TRUE))

  meta_pvalue <- 2 * pnorm(-abs(z_comb))

  k_vec <- rowSums(valid)
  meta_pvalue[k_vec < 2] <- NA_real_

  meta_log2FC <- rowMeans(lfc_mat, na.rm = TRUE)
  list(meta_pvalue = meta_pvalue, meta_log2FC = meta_log2FC)
}

# --- Vectorized Inverse Normal ---
.meta_inverse_normal <- function(pval_mat, lfc_mat, n_samples) {
  valid <- !is.na(pval_mat) & !is.na(lfc_mat)
  z_mat <- matrix(NA_real_, nrow = nrow(pval_mat), ncol = ncol(pval_mat))
  z_mat[valid] <- qnorm(1 - pval_mat[valid]) * sign(lfc_mat[valid])

  # Weights: sqrt(n) * |lfc|, with floor at eps
  w_mat <- matrix(sqrt(n_samples), nrow = nrow(pval_mat),
                   ncol = ncol(pval_mat), byrow = TRUE)
  w_mat <- w_mat * abs(lfc_mat)
  w_mat[w_mat == 0] <- .Machine$double.eps
  w_mat[!valid] <- NA_real_

  z_fused <- rowSums(w_mat * z_mat, na.rm = TRUE) /
             sqrt(rowSums(w_mat^2, na.rm = TRUE))

  meta_pvalue <- 2 * pnorm(-abs(z_fused))

  k_vec <- rowSums(valid)
  meta_pvalue[k_vec < 2] <- NA_real_

  meta_log2FC <- rowMeans(lfc_mat, na.rm = TRUE)
  list(meta_pvalue = meta_pvalue, meta_log2FC = meta_log2FC)
}

# --- Vectorized Fixed Effects ---
.meta_fixed_effects <- function(lfc_mat, se_mat) {
  valid <- !is.na(lfc_mat) & !is.na(se_mat) & se_mat > 0
  lfc_v <- lfc_mat; lfc_v[!valid] <- NA_real_
  se_v  <- se_mat;  se_v[!valid]  <- NA_real_

  w <- 1 / se_v^2
  w[!valid] <- NA_real_

  sum_w    <- rowSums(w, na.rm = TRUE)
  theta_FE <- rowSums(w * lfc_v, na.rm = TRUE) / sum_w
  SE_FE    <- 1 / sqrt(sum_w)
  Z_FE     <- theta_FE / SE_FE

  meta_pvalue <- 2 * pnorm(-abs(Z_FE))

  k_vec <- rowSums(valid)
  meta_pvalue[k_vec < 2] <- NA_real_
  theta_FE[k_vec < 2]    <- NA_real_

  list(meta_pvalue = meta_pvalue, meta_log2FC = theta_FE)
}

# --- Vectorized Random Effects ---
.meta_random_effects <- function(lfc_mat, se_mat, tau_sq_vec) {
  valid <- !is.na(lfc_mat) & !is.na(se_mat) & se_mat > 0
  lfc_v <- lfc_mat; lfc_v[!valid] <- NA_real_
  se_v  <- se_mat;  se_v[!valid]  <- NA_real_

  tau_sq_vec[is.na(tau_sq_vec)] <- 0

  # w_re = 1 / (se^2 + tau_sq), broadcast tau_sq across columns
  w_re <- 1 / (se_v^2 + tau_sq_vec)
  w_re[!valid] <- NA_real_

  sum_w_re <- rowSums(w_re, na.rm = TRUE)
  theta_RE <- rowSums(w_re * lfc_v, na.rm = TRUE) / sum_w_re
  SE_RE    <- 1 / sqrt(sum_w_re)
  Z_RE     <- theta_RE / SE_RE

  meta_pvalue <- 2 * pnorm(-abs(Z_RE))

  k_vec <- rowSums(valid)
  meta_pvalue[k_vec < 2] <- NA_real_
  theta_RE[k_vec < 2]    <- NA_real_

  list(meta_pvalue = meta_pvalue, meta_log2FC = theta_RE)
}

# --- Vectorized AWmeta ---
.meta_awmeta <- function(pval_mat, lfc_mat, se_mat, i_sq_vec) {
  i_sq_vec[is.na(i_sq_vec)] <- 50
  lambda <- 1 - i_sq_vec / 100

  # --- Effect-size Z (inverse-variance weighted) ---
  valid_es <- !is.na(lfc_mat) & !is.na(se_mat) & se_mat > 0
  lfc_v <- lfc_mat; lfc_v[!valid_es] <- NA_real_
  se_v  <- se_mat;  se_v[!valid_es]  <- NA_real_

  w_re <- 1 / se_v^2
  w_re[!valid_es] <- NA_real_

  sum_w_re <- rowSums(w_re, na.rm = TRUE)
  theta_RE <- rowSums(w_re * lfc_v, na.rm = TRUE) / sum_w_re
  SE_RE    <- 1 / sqrt(sum_w_re)
  Z_effect <- theta_RE / SE_RE

  # --- P-value Z (Stouffer-style mean) ---
  valid_pv <- !is.na(pval_mat) & !is.na(lfc_mat)
  z_pv_mat <- matrix(NA_real_, nrow = nrow(pval_mat), ncol = ncol(pval_mat))
  z_pv_mat[valid_pv] <- qnorm(1 - pval_mat[valid_pv] / 2) *
                         sign(lfc_mat[valid_pv])

  k_pv     <- rowSums(valid_pv)
  Z_pvalue <- rowSums(z_pv_mat, na.rm = TRUE) / k_pv

  # --- Combine ---
  Z_adaptive  <- lambda * Z_effect + (1 - lambda) * Z_pvalue
  meta_pvalue <- 2 * pnorm(-abs(Z_adaptive))

  k_es <- rowSums(valid_es)
  too_few <- k_es < 2 | k_pv < 2
  meta_pvalue[too_few] <- NA_real_
  theta_RE[too_few]    <- NA_real_

  list(meta_pvalue = meta_pvalue, meta_log2FC = theta_RE)
}
