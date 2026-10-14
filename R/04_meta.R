#' @include AllClasses.R
NULL

# ============================================================================
# Module 4: Meta-Analysis Statistics
# ============================================================================

.mx_meta_methods <- c("random_effects", "fixed_effects", "fisher", "stouffer",
                      "inverse_normal", "minp", "maxp", "rop", "awmeta")

#' Run meta-analysis across per-study DE results
#'
#' Integrates per-study differential expression results from
#' \code{\link{mx_de_all}} into a single gene-level meta-analysis. Two
#' families of methods are available: \emph{effect-size} models, which pool
#' the per-study log2 fold changes weighted by their sampling variance, and
#' \emph{p-value combination} methods, which combine per-study evidence
#' without pooling effects.
#'
#' @param de_results A named list of \code{data.frame} objects (the
#'   \code{de_result} slots produced by \code{\link{mx_de}}), or a named list
#'   of \code{\linkS4class{metaXpressStudy}} objects after
#'   \code{\link{mx_de_all}}. Required columns: \code{gene_id},
#'   \code{log2FC}, \code{pvalue}, \code{padj}; effect-size methods also need
#'   \code{lfcSE}.
#' @param method Character scalar. One of:
#'   \describe{
#'     \item{\code{"random_effects"}}{(default) Inverse-variance random-effects
#'       model; between-study variance \eqn{\tau^2} estimated by
#'       \code{tau_method}. Estimates the \emph{average} effect across
#'       populations and reports a prediction interval.}
#'     \item{\code{"fixed_effects"}}{Inverse-variance common-effect model.
#'       Appropriate only when studies are near-replicates (low \eqn{I^2}).}
#'     \item{\code{"fisher"}}{Fisher's sum of \eqn{-2\log p}. With
#'       \code{directional = TRUE} (default) the two one-sided Fisher tests
#'       are computed and combined as in Owen (2009), so evidence in
#'       opposite directions does not accumulate. Like all Fisher-type tests
#'       it can still be driven by one extreme study; check
#'       \code{direction_consistency}.}
#'     \item{\code{"stouffer"}}{Unweighted Stouffer Z on signed z-scores.}
#'     \item{\code{"inverse_normal"}}{Weighted inverse-normal (Stouffer)
#'       combination with weights \eqn{\sqrt{n_s}}, as in Marot et al. (2009)
#'       and Rau et al. (2014).}
#'     \item{\code{"minp"}}{Tippett's minimum p: DE in \emph{at least one}
#'       study. Non-directional.}
#'     \item{\code{"maxp"}}{Wilkinson's maximum p: DE in \emph{all} studies.
#'       Non-directional.}
#'     \item{\code{"rop"}}{r-th ordered p-value (Song & Tseng 2014): DE in
#'       at least \code{r} studies. Non-directional.}
#'     \item{\code{"awmeta"}}{AWmeta (Hu et al. 2026). The p-value is the
#'       adaptively weighted Fisher statistic (AW-Fisher; Li & Tseng 2011),
#'       calibrated with the \pkg{AWFisher} package (Huo et al. 2020), which
#'       also returns binary per-study weights indicating which studies
#'       contribute. The effect size is an REML random-effects estimate over
#'       the selected studies only (AW-REM). Requires \pkg{AWFisher} and
#'       \code{lfcSE}. Non-directional, as in the original method; check
#'       \code{direction_consistency}. The selected studies are reported in
#'       \code{aw_studies}.}
#'   }
#' @param min_studies Integer. Minimum number of studies in which a gene must
#'   be measured. Default: \code{2}.
#' @param tau_method Estimator of the between-study variance for
#'   \code{"random_effects"}: \code{"REML"} (default; restricted maximum
#'   likelihood), \code{"PM"} (Paule-Mandel) or \code{"DL"}
#'   (DerSimonian-Laird). REML and PM are preferred over DL, which
#'   underestimates \eqn{\tau^2} (Veroniki et al. 2016; Langan et al. 2019).
#' @param hksj Logical. Apply the Hartung-Knapp-Sidik-Jonkman adjustment to
#'   the random-effects standard error and use a \eqn{t_{k-1}} reference
#'   distribution. Recommended when few studies are pooled because the
#'   standard Wald test is anti-conservative (IntHout et al. 2014). Default:
#'   \code{FALSE}.
#' @param directional Logical. For \code{"fisher"}, combine one-sided
#'   p-values in each direction (default \code{TRUE}). \code{FALSE} gives the
#'   classical Fisher test on two-sided p-values, which ignores sign.
#' @param weights For \code{"stouffer"}-type methods and the effect estimate
#'   reported by p-value methods: \code{NULL} (default) uses per-study sample
#'   sizes from the \code{n_case}/\code{n_control} columns written by
#'   \code{\link{mx_de}} for \code{"inverse_normal"} and equal weights for
#'   \code{"stouffer"}. A numeric vector (one per study) overrides this.
#' @param r Integer. Order statistic for \code{"rop"}. Default:
#'   \code{ceiling(k / 2)} for a gene measured in \code{k} studies.
#' @param alpha Numeric. Confidence level is \code{1 - alpha} for confidence
#'   and prediction intervals. Default: \code{0.05}.
#' @param n_samples Deprecated alias for \code{weights} supplied as total
#'   sample sizes; its square root is used.
#' @param ... Unused.
#'
#' @return A \code{\linkS4class{metaXpressResult}}. Its \code{meta_table} has
#'   one row per gene with columns:
#'   \describe{
#'     \item{gene_id}{Gene identifier.}
#'     \item{meta_log2FC, meta_se}{Pooled log2 fold change and its standard
#'       error. For p-value methods this is the inverse-variance common-effect
#'       estimate, reported as a descriptive effect size.}
#'     \item{ci_lower, ci_upper}{\code{1 - alpha} confidence interval.}
#'     \item{pi_lower, pi_upper}{Random-effects prediction interval for the
#'       effect in a new study (Higgins et al. 2009); requires \eqn{k \ge 3}.}
#'     \item{meta_stat}{Test statistic (Z, t or \eqn{\chi^2} as appropriate).}
#'     \item{meta_pvalue, meta_padj}{Combined p-value and its
#'       Benjamini-Hochberg adjustment across genes.}
#'     \item{tau2, i_squared, q_stat, q_pvalue}{Between-study variance,
#'       \eqn{I^2} (\%), Cochran's Q and its p-value.}
#'     \item{n_studies, n_up, n_down}{Studies measuring the gene, and those
#'       with positive / negative log2FC.}
#'     \item{direction_consistency}{Fraction of studies whose log2FC sign
#'       agrees with \code{meta_log2FC}.}
#'   }
#'
#' @details
#' Use \code{meta_padj} (BH-adjusted) for significance. For a typical
#' disease-signature meta-analysis we recommend \code{"random_effects"} with
#' \code{tau_method = "REML"}, and reporting \eqn{I^2} and the prediction
#' interval alongside the pooled effect. P-value combination methods answer
#' different questions (any vs. all vs. most studies) and should be chosen
#' a priori, not by comparing numbers of significant genes.
#'
#' A study contributes to a gene only when it reports \code{log2FC} and
#' \code{pvalue} (plus \code{lfcSE} for effect-size methods). DESeq2 sets
#' \code{pvalue} to \code{NA} for count outliers flagged by Cook's distance,
#' so such genes are treated as not measured in that study. (DESeq2's
#' independent filtering affects only \code{padj}, which is not used.)
#'
#' All methods assume studies are independent. Fold-change direction must be
#' harmonised across studies (same case/control orientation); check with
#' \code{\link{mx_study_concordance}}.
#'
#' @references
#' DerSimonian, R. & Laird, N. (1986) Meta-analysis in clinical trials.
#' \emph{Controlled Clinical Trials}, \strong{7}(3), 177--188.
#' \doi{10.1016/0197-2456(86)90046-2}
#'
#' Viechtbauer, W. (2005) Bias and efficiency of meta-analytic variance
#' estimators in the random-effects model. \emph{Journal of Educational and
#' Behavioral Statistics}, \strong{30}(3), 261--293.
#' \doi{10.3102/10769986030003261}
#'
#' Paule, R.C. & Mandel, J. (1982) Consensus values and weighting factors.
#' \emph{Journal of Research of the National Bureau of Standards},
#' \strong{87}(5), 377--385. \doi{10.6028/jres.087.022}
#'
#' Hartung, J. & Knapp, G. (2001) A refined method for the meta-analysis of
#' controlled clinical trials with binary outcome. \emph{Statistics in
#' Medicine}, \strong{20}(24), 3875--3889. \doi{10.1002/sim.1009}
#'
#' IntHout, J., Ioannidis, J.P.A. & Borm, G.F. (2014) The
#' Hartung-Knapp-Sidik-Jonkman method for random effects meta-analysis is
#' straightforward and considerably outperforms the standard
#' DerSimonian-Laird method. \emph{BMC Medical Research Methodology},
#' \strong{14}, 25. \doi{10.1186/1471-2288-14-25}
#'
#' Higgins, J.P.T., Thompson, S.G. & Spiegelhalter, D.J. (2009) A
#' re-evaluation of random-effects meta-analysis. \emph{Journal of the Royal
#' Statistical Society: Series A}, \strong{172}(1), 137--159.
#' \doi{10.1111/j.1467-985X.2008.00552.x}
#'
#' Higgins, J.P.T. & Thompson, S.G. (2002) Quantifying heterogeneity in a
#' meta-analysis. \emph{Statistics in Medicine}, \strong{21}(11), 1539--1558.
#' \doi{10.1002/sim.1186}
#'
#' Veroniki, A.A. et al. (2016) Methods to estimate the between-study
#' variance and its uncertainty in meta-analysis. \emph{Research Synthesis
#' Methods}, \strong{7}(1), 55--79. \doi{10.1002/jrsm.1164}
#'
#' Langan, D. et al. (2019) A comparison of heterogeneity variance
#' estimators in simulated random-effects meta-analyses. \emph{Research
#' Synthesis Methods}, \strong{10}(1), 83--98. \doi{10.1002/jrsm.1316}
#'
#' Owen, A.B. (2009) Karl Pearson's meta-analysis revisited. \emph{Annals of
#' Statistics}, \strong{37}(6B), 3867--3892. \doi{10.1214/09-AOS697}
#'
#' Marot, G., Foulley, J.-L., Mayer, C.-D. & Jaffrezic, F. (2009) Moderated
#' effect size and P-value combinations for microarray meta-analyses.
#' \emph{Bioinformatics}, \strong{25}(20), 2692--2699.
#' \doi{10.1093/bioinformatics/btp444}
#'
#' Rau, A., Marot, G. & Jaffrezic, F. (2014) Differential meta-analysis of
#' RNA-seq data from multiple studies. \emph{BMC Bioinformatics},
#' \strong{15}, 91. \doi{10.1186/1471-2105-15-91}
#'
#' Song, C. & Tseng, G.C. (2014) Hypothesis setting and order statistic for
#' robust genomic meta-analysis. \emph{Annals of Applied Statistics},
#' \strong{8}(2), 777--800. \doi{10.1214/13-AOAS683}
#'
#' Li, J. & Tseng, G.C. (2011) An adaptively weighted statistic for detecting
#' differential gene expression when combining multiple transcriptomic
#' studies. \emph{Annals of Applied Statistics}, \strong{5}(2A), 994--1019.
#' \doi{10.1214/10-AOAS393}
#'
#' Huo, Z., Tang, S., Park, Y. & Tseng, G. (2020) P-value evaluation,
#' variability index and biomarker categorization for adaptively weighted
#' Fisher's meta-analysis method in omics applications.
#' \emph{Bioinformatics}, \strong{36}(2), 524--532.
#' \doi{10.1093/bioinformatics/btz589}
#'
#' Hu, Y. et al. (2026) AWmeta empowers adaptively weighted transcriptomic
#' meta-analysis. \emph{Current Issues in Molecular Biology},
#' \strong{48}(5), 530. \doi{10.3390/cimb48050530}
#'
#' @examples
#' # Three per-study results in the format written by mx_de()
#' set.seed(1)
#' genes <- paste0("G", 1:200)
#' de <- lapply(1:3, function(s) {
#'   lfc <- rnorm(200, c(rep(1, 20), rep(0, 180)), 0.3)
#'   se  <- runif(200, 0.2, 0.4)
#'   p   <- 2 * pnorm(-abs(lfc / se))
#'   data.frame(gene_id = genes, log2FC = lfc, lfcSE = se, pvalue = p,
#'              padj = p.adjust(p, "BH"), n_case = 5, n_control = 5)
#' })
#' names(de) <- c("S1", "S2", "S3")
#' res <- mx_meta(de, method = "random_effects", tau_method = "REML")
#' head(res@meta_table[order(res@meta_table$meta_pvalue), ])
#'
#' @seealso \code{\link{mx_heterogeneity}}, \code{\link{mx_sensitivity}},
#'   \code{\link{mx_study_concordance}}, \code{\link{mx_forest}}
#'
#' @importFrom methods new
#' @importFrom stats p.adjust pchisq pnorm qnorm pt qt pbeta
#' @export
mx_meta <- function(de_results,
                    method      = c("random_effects", "fixed_effects",
                                    "fisher", "stouffer", "inverse_normal",
                                    "minp", "maxp", "rop", "awmeta"),
                    min_studies = 2,
                    tau_method  = c("REML", "PM", "DL"),
                    hksj        = FALSE,
                    directional = TRUE,
                    weights     = NULL,
                    r           = NULL,
                    alpha       = 0.05,
                    n_samples   = NULL,
                    ...) {
  de_results <- .extract_de_results(de_results)
  .validate_de_list(de_results)
  method     <- match.arg(method)
  tau_method <- match.arg(tau_method)

  if (!is.numeric(min_studies) || length(min_studies) != 1 || min_studies < 2)
    stop("'min_studies' must be a single integer >= 2")
  if (is.null(names(de_results)))
    names(de_results) <- paste0("study_", seq_along(de_results))
  k_total <- length(de_results)
  if (k_total < 2) stop("At least 2 studies are required")

  if (method == "awmeta") .require_package("AWFisher", "AWmeta")
  if (method %in% c("random_effects", "fixed_effects", "awmeta")) {
    no_se <- !vapply(de_results, function(d) "lfcSE" %in% colnames(d),
                     logical(1))
    if (any(no_se))
      stop("Effect-size methods need an 'lfcSE' column; missing in: ",
           paste(names(de_results)[no_se], collapse = ", "),
           ". Re-run mx_de() or use a p-value method.")
  }
  if (method == "random_effects" && k_total < 3)
    warning("'random_effects' is unreliable with fewer than 3 studies: ",
            "tau^2 cannot be estimated precisely.")

  message("Running meta-analysis (", method, ") on ", k_total, " studies...")

  # --- Aligned gene x study matrices --------------------------------------
  mats <- .build_gene_matrices(de_results)
  lfc  <- mats$lfc_mat
  se   <- mats$se_mat
  pv   <- mats$pval_mat
  se[!is.finite(se) | se <= 0] <- NA_real_
  pv[!is.na(pv)] <- pmin(pmax(pv[!is.na(pv)], .Machine$double.xmin), 1)

  measured <- !is.na(lfc) & !is.na(pv)
  if (method %in% c("random_effects", "fixed_effects", "awmeta"))
    measured <- measured & !is.na(se)
  lfc[!measured] <- NA_real_
  se[!measured]  <- NA_real_
  pv[!measured]  <- NA_real_

  k_gene <- rowSums(measured)
  keep   <- k_gene >= min_studies
  if (!any(keep))
    stop("No genes measured in >= ", min_studies, " studies. ",
         "Lower 'min_studies' or check gene ID alignment.")
  lfc <- lfc[keep, , drop = FALSE]
  se  <- se[keep, , drop = FALSE]
  pv  <- pv[keep, , drop = FALSE]
  k_gene <- k_gene[keep]

  # --- Study weights for Stouffer-type combinations ------------------------
  study_w <- .study_weights(de_results, method, weights, n_samples)

  # --- Heterogeneity (always reported when SEs exist) ----------------------
  het <- .heterogeneity(lfc, se, tau_method)

  zcrit <- qnorm(1 - alpha / 2)
  res <- switch(method,
    random_effects = .meta_re(lfc, se, het$tau2, hksj, alpha),
    fixed_effects  = .meta_fe(lfc, se, alpha),
    fisher         = .meta_fisher(pv, lfc, directional),
    stouffer       = .meta_stouffer(pv, lfc, study_w),
    inverse_normal = .meta_stouffer(pv, lfc, study_w),
    minp           = .meta_order(pv, rep(1L, nrow(pv))),
    maxp           = .meta_order(pv, k_gene),
    rop            = .meta_order(pv, if (is.null(r)) ceiling(k_gene / 2)
                                      else pmin(as.integer(r), k_gene)),
    awmeta         = .meta_awmeta(pv, lfc, se, alpha)
  )

  # P-value methods: report a descriptive pooled effect
  if (is.null(res$est)) {
    if (all(is.na(se))) {
      w <- matrix(study_w, nrow(lfc), ncol(lfc), byrow = TRUE)
      w[is.na(lfc)] <- NA_real_
      res$est <- rowSums(w * lfc, na.rm = TRUE) / rowSums(w, na.rm = TRUE)
      res$se  <- rep(NA_real_, nrow(lfc))
    } else {
      fe <- .meta_fe(lfc, se, alpha)
      res$est <- fe$est
      res$se  <- fe$se
    }
    res$ci_lower <- res$est - zcrit * res$se
    res$ci_upper <- res$est + zcrit * res$se
  }

  pi <- if (method == "random_effects") {
    .prediction_interval(res$est, res$se, het$tau2, k_gene, alpha)
  } else {
    list(lower = rep(NA_real_, nrow(lfc)), upper = rep(NA_real_, nrow(lfc)))
  }

  s_meta <- sign(res$est)
  n_up   <- rowSums(lfc > 0, na.rm = TRUE)
  n_down <- rowSums(lfc < 0, na.rm = TRUE)

  meta_table <- data.frame(
    gene_id               = rownames(lfc),
    meta_log2FC           = res$est,
    meta_se               = res$se,
    ci_lower              = res$ci_lower,
    ci_upper              = res$ci_upper,
    pi_lower              = pi$lower,
    pi_upper              = pi$upper,
    meta_stat             = res$stat,
    meta_pvalue           = res$p,
    meta_padj             = p.adjust(res$p, method = "BH"),
    tau2                  = het$tau2,
    i_squared             = het$I_sq,
    q_stat                = het$Q,
    q_pvalue              = het$p_heterogeneity,
    n_studies             = as.integer(k_gene),
    n_up                  = as.integer(n_up),
    n_down                = as.integer(n_down),
    direction_consistency = ifelse(s_meta > 0, n_up, n_down) / k_gene,
    stringsAsFactors      = FALSE,
    row.names             = NULL
  )
  if (method == "awmeta") {
    w <- res$weights
    meta_table$aw_n_selected <- as.integer(rowSums(w == 1, na.rm = TRUE))
    meta_table$aw_studies <- apply(w, 1, function(x)
      paste(colnames(lfc)[!is.na(x) & x == 1], collapse = ";"))
  }

  heterogeneity_df <- data.frame(
    gene_id         = rownames(lfc),
    Q               = het$Q,
    df              = het$df,
    I_sq            = het$I_sq,
    tau_sq          = het$tau2,
    p_heterogeneity = het$p_heterogeneity,
    stringsAsFactors = FALSE,
    row.names        = NULL
  )

  out <- new("metaXpressResult",
             meta_table     = meta_table,
             method         = method,
             n_studies      = as.integer(k_total),
             heterogeneity  = heterogeneity_df,
             pathway_result = data.frame())
  attr(out@meta_table, "settings") <- list(
    method = method, tau_method = tau_method, hksj = hksj,
    directional = directional, min_studies = min_studies, alpha = alpha,
    studies = names(de_results), weights = study_w)
  out
}

#' Compute per-gene heterogeneity statistics
#'
#' Returns Cochran's Q, degrees of freedom, \eqn{I^2}, the between-study
#' variance \eqn{\tau^2} and the p-value of the Q test for each gene measured
#' in at least two studies.
#'
#' @inheritParams mx_meta
#'
#' @return A \code{data.frame} with columns \code{gene_id}, \code{Q},
#'   \code{df}, \code{I_sq}, \code{tau_sq}, \code{p_heterogeneity}.
#'
#' @details \eqn{I^2 = \max(0, (Q - df)/Q)} describes the proportion of
#' variability due to between-study differences rather than sampling error
#' (Higgins & Thompson 2002). With few studies the Q test has low power and
#' \eqn{I^2} is imprecise, so heterogeneity should be interpreted alongside
#' \eqn{\tau^2} and the prediction interval from \code{\link{mx_meta}}.
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
#' set.seed(1)
#' de <- lapply(1:3, function(s) {
#'   lfc <- rnorm(50, s * 0.3, 0.2); se <- rep(0.2, 50)
#'   p <- 2 * pnorm(-abs(lfc / se))
#'   data.frame(gene_id = paste0("G", 1:50), log2FC = lfc, lfcSE = se,
#'              pvalue = p, padj = p.adjust(p, "BH"))
#' })
#' het <- mx_heterogeneity(de)
#' summary(het$I_sq)
#'
#' @seealso \code{\link{mx_meta}}, \code{\link{mx_heterogeneity_plot}}
#' @export
mx_heterogeneity <- function(de_results, tau_method = c("REML", "PM", "DL")) {
  de_results <- .extract_de_results(de_results)
  .validate_de_list(de_results)
  tau_method <- match.arg(tau_method)

  mats <- .build_gene_matrices(de_results)
  se   <- mats$se_mat
  se[!is.finite(se) | se <= 0] <- NA_real_
  het  <- .heterogeneity(mats$lfc_mat, se, tau_method)

  data.frame(
    gene_id         = mats$all_genes,
    Q               = het$Q,
    df              = het$df,
    I_sq            = het$I_sq,
    tau_sq          = het$tau2,
    p_heterogeneity = het$p_heterogeneity,
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
}

#' Leave-one-out sensitivity analysis
#'
#' Reruns the meta-analysis \eqn{k} times, each time excluding one study, to
#' assess whether any single study drives the result (Cochrane Handbook,
#' Section 10.14).
#'
#' @inheritParams mx_meta
#' @param padj_threshold Numeric. FDR threshold used to call genes in the
#'   robustness summary. Default: \code{0.05}.
#' @param ... Further arguments passed to \code{\link{mx_meta}} (e.g.
#'   \code{tau_method}, \code{hksj}).
#'
#' @return A named list of \code{\linkS4class{metaXpressResult}} objects, one
#'   per excluded study. The attribute \code{"robustness"} holds a
#'   \code{data.frame} with, for every gene, the number and fraction of
#'   leave-one-out runs in which it remained significant
#'   (\code{meta_padj <= padj_threshold}) with the same direction as in the
#'   full analysis.
#'
#' @references
#' Higgins, J.P.T. et al. (eds.) (2019) \emph{Cochrane Handbook for
#' Systematic Reviews of Interventions}. 2nd edn. Chichester: John Wiley
#' & Sons. \doi{10.1002/9781119536604}
#'
#' @examples
#' set.seed(1)
#' de <- lapply(1:4, function(s) {
#'   lfc <- rnorm(100, c(rep(1, 10), rep(0, 90)), 0.2); se <- rep(0.25, 100)
#'   p <- 2 * pnorm(-abs(lfc / se))
#'   data.frame(gene_id = paste0("G", 1:100), log2FC = lfc, lfcSE = se,
#'              pvalue = p, padj = p.adjust(p, "BH"))
#' })
#' names(de) <- paste0("S", 1:4)
#' loo <- mx_sensitivity(de, method = "random_effects")
#' head(attr(loo, "robustness"))
#'
#' @seealso \code{\link{mx_meta}}
#' @export
mx_sensitivity <- function(de_results,
                           method = c("random_effects", "fixed_effects",
                                      "fisher", "stouffer", "inverse_normal",
                                      "minp", "maxp", "rop", "awmeta"),
                           padj_threshold = 0.05,
                           ...) {
  de_results <- .extract_de_results(de_results)
  method <- match.arg(method)

  if (length(de_results) < 3)
    stop("Sensitivity analysis requires at least 3 studies")
  if (is.null(names(de_results)))
    names(de_results) <- paste0("study_", seq_along(de_results))

  full <- suppressMessages(mx_meta(de_results, method = method, ...))
  message("Running leave-one-out sensitivity analysis...")

  loo <- lapply(seq_along(de_results), function(i) {
    message("  Excluding: ", names(de_results)[i])
    suppressMessages(suppressWarnings(
      mx_meta(de_results[-i], method = method, ...)))
  })
  names(loo) <- names(de_results)

  ft  <- full@meta_table
  hit <- vapply(loo, function(res) {
    mt <- res@meta_table[match(ft$gene_id, res@meta_table$gene_id), ]
    !is.na(mt$meta_padj) & mt$meta_padj <= padj_threshold &
      sign(mt$meta_log2FC) == sign(ft$meta_log2FC)
  }, logical(nrow(ft)))
  hit <- matrix(hit, nrow = nrow(ft))

  attr(loo, "robustness") <- data.frame(
    gene_id          = ft$gene_id,
    full_padj        = ft$meta_padj,
    n_loo_sig        = rowSums(hit),
    frac_loo_sig     = rowMeans(hit),
    stringsAsFactors = FALSE
  )
  loo
}

#' Check fold-change concordance between studies
#'
#' Computes pairwise correlations of per-study log2 fold changes over genes
#' measured in both studies. A negative correlation usually means the
#' case/control orientation is flipped in one study, which would invalidate
#' any directional meta-analysis; low correlation flags biologically or
#' technically divergent studies.
#'
#' @inheritParams mx_meta
#' @param genes Which genes to correlate: \code{"significant"} (default;
#'   genes with \code{padj <= padj_threshold} in at least one of the two
#'   studies) or \code{"all"}.
#' @param padj_threshold Numeric. Default: \code{0.05}.
#' @param cor_method Correlation method. Default: \code{"spearman"}.
#'
#' @return A \code{data.frame} with one row per study pair: \code{study_1},
#'   \code{study_2}, \code{n_genes}, \code{correlation}, \code{sign_agreement}
#'   (fraction of genes with the same log2FC sign) and a \code{flag}.
#'
#' @examples
#' set.seed(1)
#' base <- rnorm(300)
#' de <- lapply(c(1, 1, -1), function(dir) {
#'   lfc <- dir * base + rnorm(300, 0, 0.3); p <- 2 * pnorm(-abs(lfc / 0.3))
#'   data.frame(gene_id = paste0("G", 1:300), log2FC = lfc, pvalue = p,
#'              padj = p.adjust(p, "BH"))
#' })
#' names(de) <- c("A", "B", "C_flipped")
#' mx_study_concordance(de)
#'
#' @seealso \code{\link{mx_meta}}, \code{\link{mx_de}}
#' @importFrom stats cor
#' @export
mx_study_concordance <- function(de_results,
                                 genes = c("significant", "all"),
                                 padj_threshold = 0.05,
                                 cor_method = c("spearman", "pearson")) {
  de_results <- .extract_de_results(de_results)
  .validate_de_list(de_results)
  genes      <- match.arg(genes)
  cor_method <- match.arg(cor_method)
  if (is.null(names(de_results)))
    names(de_results) <- paste0("study_", seq_along(de_results))

  lfc  <- .build_gene_matrix(de_results, "log2FC")
  padj <- .build_gene_matrix(de_results, "padj")
  pairs <- utils::combn(colnames(lfc), 2)

  rows <- lapply(seq_len(ncol(pairs)), function(j) {
    a <- pairs[1, j]; b <- pairs[2, j]
    use <- !is.na(lfc[, a]) & !is.na(lfc[, b])
    if (genes == "significant")
      use <- use & ((!is.na(padj[, a]) & padj[, a] <= padj_threshold) |
                    (!is.na(padj[, b]) & padj[, b] <= padj_threshold))
    n  <- sum(use)
    rr <- if (n >= 10) cor(lfc[use, a], lfc[use, b], method = cor_method)
          else NA_real_
    sa <- if (n > 0) mean(sign(lfc[use, a]) == sign(lfc[use, b])) else NA_real_
    flag <- if (is.na(rr)) "too few genes"
            else if (rr < 0) "NEGATIVE: check contrast orientation"
            else if (rr < 0.1) "weak concordance"
            else "ok"
    data.frame(study_1 = a, study_2 = b, n_genes = n, correlation = rr,
               sign_agreement = sa, flag = flag, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  neg <- !is.na(out$correlation) & out$correlation < 0
  if (any(neg))
    warning("Negative fold-change correlation between studies: ",
            paste(out$study_1[neg], out$study_2[neg], sep = " vs ",
                  collapse = "; "),
            ". Check that every study uses the same case/control contrast.")
  out
}

# ============================================================================
# Internal statistical helpers -- vectorized over genes (rows)
# ============================================================================

# Per-study weights for Stouffer-type combination
.study_weights <- function(de_results, method, weights, n_samples) {
  k <- length(de_results)
  if (!is.null(weights)) {
    if (length(weights) != k || any(!is.finite(weights)) || any(weights <= 0))
      stop("'weights' must be ", k, " positive numbers (one per study)")
    return(as.numeric(weights))
  }
  if (!is.null(n_samples)) {
    if (length(n_samples) != k) stop("'n_samples' must have length ", k)
    return(sqrt(as.numeric(n_samples)))
  }
  n <- vapply(de_results, function(d) {
    if (all(c("n_case", "n_control") %in% colnames(d)))
      as.numeric(d$n_case[1] + d$n_control[1]) else NA_real_
  }, numeric(1))
  if (method == "inverse_normal") {
    if (anyNA(n))
      stop("'inverse_normal' needs per-study sample sizes: supply 'weights' ",
           "or run mx_de(), which records n_case/n_control")
    return(sqrt(n))
  }
  if (method %in% c("fisher", "minp", "maxp", "rop", "stouffer") || anyNA(n))
    return(rep(1, k))
  sqrt(n)
}

# Inverse-variance quantities shared by several estimators
.iv <- function(y, v, tau2 = 0) {
  w   <- 1 / (v + tau2)
  w[is.na(y)] <- NA_real_
  sw  <- rowSums(w, na.rm = TRUE)
  mu  <- rowSums(w * y, na.rm = TRUE) / sw
  list(w = w, sw = sw, mu = mu)
}

# Q, I^2 and tau^2 (DL / REML / PM)
.heterogeneity <- function(lfc, se, tau_method = "REML") {
  y <- lfc; v <- se^2
  ok <- !is.na(y) & !is.na(v)
  y[!ok] <- NA_real_; v[!ok] <- NA_real_
  k  <- rowSums(ok)
  df <- k - 1

  fe <- .iv(y, v)
  Q  <- rowSums(fe$w * (y - fe$mu)^2, na.rm = TRUE)
  sw2 <- rowSums(fe$w^2, na.rm = TRUE)
  C   <- fe$sw - sw2 / fe$sw
  tau2_dl <- pmax(0, (Q - df) / C)

  tau2 <- switch(tau_method,
    DL   = tau2_dl,
    REML = .tau2_reml(y, v, tau2_dl),
    PM   = .tau2_pm(y, v, df)
  )
  I_sq <- ifelse(Q > 0, pmax(0, (Q - df) / Q) * 100, 0)
  p_het <- pchisq(Q, df = df, lower.tail = FALSE)

  few <- k < 2
  Q[few] <- NA_real_; I_sq[few] <- NA_real_; tau2[few] <- NA_real_
  p_het[few] <- NA_real_; df[few] <- NA_real_
  list(Q = Q, df = df, I_sq = I_sq, tau2 = tau2, p_heterogeneity = p_het)
}

# REML by Fisher scoring (Viechtbauer 2005), started from DL
.tau2_reml <- function(y, v, start, maxit = 200, tol = 1e-10) {
  tau2 <- start
  tau2[is.na(tau2)] <- 0
  active <- rowSums(!is.na(y)) >= 2
  for (it in seq_len(maxit)) {
    iv  <- .iv(y, v, tau2)
    sw2 <- rowSums(iv$w^2, na.rm = TRUE)
    sw3 <- rowSums(iv$w^3, na.rm = TRUE)
    # Score and expected information of the REML log-likelihood, with
    # P = W - W 1 1' W / sum(W): score = -tr(P)/2 + y'PPy/2, info = tr(PP)/2
    score <- -0.5 * (iv$sw - sw2 / iv$sw) +
      0.5 * rowSums(iv$w^2 * (y - iv$mu)^2, na.rm = TRUE)
    info  <- 0.5 * (sw2 - 2 * sw3 / iv$sw + sw2^2 / iv$sw^2)
    new   <- pmax(0, tau2 + score / info)
    new[!active] <- NA_real_
    delta <- abs(new - tau2)
    tau2  <- new
    if (all(delta[active] < tol, na.rm = TRUE)) break
  }
  tau2
}

# Paule-Mandel: solve generalized Q(tau2) = k - 1 by bisection
.tau2_pm <- function(y, v, df, iters = 100) {
  qgen <- function(t2) {
    iv <- .iv(y, v, t2)
    rowSums(iv$w * (y - iv$mu)^2, na.rm = TRUE)
  }
  lo <- rep(0, nrow(y))
  hi <- pmax(1, 10 * apply(y, 1, stats::var, na.rm = TRUE), na.rm = TRUE)
  hi[!is.finite(hi)] <- 1
  # Expand the upper bound until Q(hi) < df
  for (j in seq_len(50)) {
    grow <- qgen(hi) > df
    grow[is.na(grow)] <- FALSE
    if (!any(grow)) break
    hi[grow] <- hi[grow] * 10
  }
  zero <- qgen(lo) <= df
  for (it in seq_len(iters)) {
    mid <- (lo + hi) / 2
    above <- qgen(mid) > df
    lo <- ifelse(above, mid, lo)
    hi <- ifelse(above, hi, mid)
  }
  out <- (lo + hi) / 2
  out[zero] <- 0
  out[df < 1] <- NA_real_
  out
}

.meta_fe <- function(lfc, se, alpha) {
  iv   <- .iv(lfc, se^2)
  est  <- iv$mu
  s    <- 1 / sqrt(iv$sw)
  z    <- est / s
  zc   <- qnorm(1 - alpha / 2)
  list(est = est, se = s, stat = z, p = 2 * pnorm(-abs(z)),
       ci_lower = est - zc * s, ci_upper = est + zc * s)
}

.meta_re <- function(lfc, se, tau2, hksj, alpha) {
  t2 <- tau2; t2[is.na(t2)] <- 0
  iv  <- .iv(lfc, se^2, t2)
  est <- iv$mu
  s   <- 1 / sqrt(iv$sw)
  k   <- rowSums(!is.na(lfc))
  if (hksj) {
    q   <- rowSums(iv$w * (lfc - est)^2, na.rm = TRUE) / (k - 1)
    s   <- sqrt(q / iv$sw)
    stat <- est / s
    p    <- 2 * pt(-abs(stat), df = k - 1)
    crit <- qt(1 - alpha / 2, df = k - 1)
  } else {
    stat <- est / s
    p    <- 2 * pnorm(-abs(stat))
    crit <- qnorm(1 - alpha / 2)
  }
  list(est = est, se = s, stat = stat, p = p,
       ci_lower = est - crit * s, ci_upper = est + crit * s)
}

# Prediction interval: est +/- t_{k-2} * sqrt(tau2 + se^2)
.prediction_interval <- function(est, se, tau2, k, alpha) {
  half <- ifelse(k >= 3,
                 qt(1 - alpha / 2, df = pmax(k - 2, 1)) * sqrt(tau2 + se^2),
                 NA_real_)
  list(lower = est - half, upper = est + half)
}

# One-sided p-values (H1: up) from two-sided p and the sign of the effect
.one_sided_up <- function(pv, lfc) {
  ifelse(lfc > 0, pv / 2, 1 - pv / 2)
}

.fisher_p <- function(p) {
  k <- rowSums(!is.na(p))
  pchisq(-2 * rowSums(log(p), na.rm = TRUE), df = 2 * k, lower.tail = FALSE)
}

.meta_fisher <- function(pv, lfc, directional) {
  if (!directional) {
    stat <- -2 * rowSums(log(pv), na.rm = TRUE)
    return(list(stat = stat, p = .fisher_p(pv)))
  }
  p_up <- .one_sided_up(pv, lfc)
  p_dn <- 1 - p_up
  p_dn[!is.na(p_dn)] <- pmax(p_dn[!is.na(p_dn)], .Machine$double.xmin)
  f_up <- .fisher_p(p_up)
  f_dn <- .fisher_p(p_dn)
  p    <- pmin(1, 2 * pmin(f_up, f_dn))
  stat <- ifelse(f_up <= f_dn,
                 -2 * rowSums(log(p_up), na.rm = TRUE),
                 2 * rowSums(log(p_dn), na.rm = TRUE))
  list(stat = stat, p = p)
}

.meta_stouffer <- function(pv, lfc, study_w) {
  z <- qnorm(.one_sided_up(pv, lfc), lower.tail = FALSE)
  w <- matrix(study_w, nrow(pv), ncol(pv), byrow = TRUE)
  w[is.na(z)] <- NA_real_
  zc <- rowSums(w * z, na.rm = TRUE) / sqrt(rowSums(w^2, na.rm = TRUE))
  list(stat = zc, p = 2 * pnorm(-abs(zc)))
}

# Order-statistic tests on two-sided p: p_(r) ~ Beta(r, k - r + 1) under H0
.meta_order <- function(pv, r) {
  k <- rowSums(!is.na(pv))
  sorted <- t(apply(pv, 1, function(x) sort(x, na.last = TRUE)))
  if (ncol(pv) == 1) sorted <- t(sorted)
  p_r <- sorted[cbind(seq_len(nrow(pv)), pmax(1L, r))]
  list(stat = p_r, p = pbeta(p_r, r, k - r + 1))
}

#' Small-study effects and publication bias tests
#'
#' Tests, gene by gene, whether smaller (less precise) studies report
#' systematically larger effects, the pattern produced by publication bias
#' or by other small-study effects. Two standard tests are provided:
#' Egger's regression test and Begg and Mazumdar's rank correlation test.
#'
#' @inheritParams mx_meta
#' @param genes Optional character vector of genes to test (e.g. the
#'   significant genes of \code{\link{mx_meta}}). Default: all genes measured
#'   in at least \code{min_studies} studies.
#' @param method \code{"egger"} (default), \code{"begg"} or \code{"both"}.
#' @param min_studies Minimum number of studies with a standard error for a
#'   gene to be tested. Default \code{3} (the minimum for Egger's test).
#'
#' @return A \code{data.frame} with \code{gene_id}, \code{n_studies} and, for
#'   Egger, \code{egger_intercept}, \code{egger_se}, \code{egger_pvalue};
#'   for Begg, \code{begg_tau}, \code{begg_pvalue}. BH-adjusted columns
#'   (\code{*_padj}) are added across the genes tested.
#'
#' @details
#' \strong{Egger's test} regresses the standardized effect
#' \eqn{y_i / s_i} on precision \eqn{1/s_i} by ordinary least squares; an
#' intercept different from zero indicates funnel-plot asymmetry
#' (t-test with \eqn{k - 2} df). \strong{Begg's test} computes Kendall's
#' \eqn{\tau} between the standardized deviates from the fixed-effect
#' estimate and their variances.
#'
#' Both tests have low power with few studies. The Cochrane Handbook advises
#' against them with fewer than 10 studies, and a warning is issued in that
#' case; a non-significant result then is not evidence of absence of bias.
#' In transcriptomic meta-analyses of public data, genome-wide summary
#' statistics are reported for all genes, so classical publication bias
#' (selective reporting of genes) does not apply; asymmetry instead points to
#' small-study differences in design, tissue or platform.
#'
#' @references
#' Egger, M., Davey Smith, G., Schneider, M. & Minder, C. (1997) Bias in
#' meta-analysis detected by a simple, graphical test. \emph{BMJ},
#' \strong{315}(7109), 629--634. \doi{10.1136/bmj.315.7109.629}
#'
#' Begg, C.B. & Mazumdar, M. (1994) Operating characteristics of a rank
#' correlation test for publication bias. \emph{Biometrics}, \strong{50}(4),
#' 1088--1101. \doi{10.2307/2533446}
#'
#' Sterne, J.A.C. et al. (2011) Recommendations for examining and
#' interpreting funnel plot asymmetry in meta-analyses of randomised
#' controlled trials. \emph{BMJ}, \strong{343}, d4002.
#' \doi{10.1136/bmj.d4002}
#'
#' @examples
#' ad <- readRDS(system.file("extdata", "ad_brain_de.rds",
#'                           package = "metaXpress"))
#' ss <- suppressWarnings(mx_small_study(ad$de_results,
#'                                       genes = c("GFAP", "VGF", "SST")))
#' ss
#'
#' @seealso \code{\link{mx_funnel}}, \code{\link{mx_meta}}
#' @importFrom stats pt cor.test
#' @export
mx_small_study <- function(de_results, genes = NULL,
                           method = c("egger", "begg", "both"),
                           min_studies = 3) {
  de_results <- .extract_de_results(de_results)
  .validate_de_list(de_results)
  method <- match.arg(method)
  if (min_studies < 3) stop("'min_studies' must be >= 3")

  mats <- .build_gene_matrices(de_results)
  y  <- mats$lfc_mat
  se <- mats$se_mat
  se[!is.finite(se) | se <= 0] <- NA_real_
  y[is.na(se)] <- NA_real_
  if (all(is.na(se)))
    stop("Small-study tests need an 'lfcSE' column in the DE results")

  k <- rowSums(!is.na(y))
  keep <- k >= min_studies
  if (!is.null(genes)) keep <- keep & rownames(y) %in% genes
  if (!any(keep)) stop("No genes measured with SE in >= ", min_studies,
                       " studies")
  y <- y[keep, , drop = FALSE]; se <- se[keep, , drop = FALSE]
  k <- k[keep]
  if (max(k) < 10)
    warning("Only ", max(k), " studies: small-study tests have very low ",
            "power with fewer than 10 studies (Sterne et al. 2011). ",
            "Interpret non-significant results with caution.")

  out <- data.frame(gene_id = rownames(y), n_studies = as.integer(k),
                    stringsAsFactors = FALSE)

  if (method %in% c("egger", "both")) {
    eg <- .egger(y, se)
    out$egger_intercept <- eg$intercept
    out$egger_se        <- eg$se
    out$egger_pvalue    <- eg$p
    out$egger_padj      <- stats::p.adjust(eg$p, "BH")
  }
  if (method %in% c("begg", "both")) {
    bg <- t(vapply(seq_len(nrow(y)), function(i) .begg(y[i, ], se[i, ]),
                   numeric(2)))
    out$begg_tau    <- bg[, 1]
    out$begg_pvalue <- bg[, 2]
    out$begg_padj   <- stats::p.adjust(bg[, 2], "BH")
  }
  out
}

# Vectorized Egger regression: z = y/se on x = 1/se (OLS), test intercept
.egger <- function(y, se) {
  z <- y / se
  x <- 1 / se
  ok <- !is.na(z)
  k  <- rowSums(ok)
  xbar <- rowSums(x, na.rm = TRUE) / k
  zbar <- rowSums(z, na.rm = TRUE) / k
  dx <- x - xbar; dz <- z - zbar
  sxx <- rowSums(dx^2, na.rm = TRUE)
  sxz <- rowSums(dx * dz, na.rm = TRUE)
  b <- sxz / sxx
  a <- zbar - b * xbar
  res <- z - (a + b * x)
  s2  <- rowSums(res^2, na.rm = TRUE) / (k - 2)
  se_a <- sqrt(s2 * (1 / k + xbar^2 / sxx))
  t <- a / se_a
  p <- 2 * pt(-abs(t), df = k - 2)
  p[k < 3 | !is.finite(p)] <- NA_real_
  list(intercept = a, se = se_a, p = p)
}

# Begg & Mazumdar rank correlation for one gene
.begg <- function(y, se) {
  ok <- !is.na(y) & !is.na(se)
  y <- y[ok]; v <- se[ok]^2
  if (length(y) < 3) return(c(NA_real_, NA_real_))
  w  <- 1 / v
  mu <- sum(w * y) / sum(w)
  vstar <- v - 1 / sum(w)
  zs <- (y - mu) / sqrt(vstar)
  ct <- suppressWarnings(cor.test(zs, v, method = "kendall",
                                  exact = length(y) < 50))
  c(unname(ct$estimate), ct$p.value)
}

# AWmeta (Hu et al. 2026): AW-Fisher p-value and weights (AWFisher package),
# REML random-effects effect size over the AW-selected studies.
.meta_awmeta <- function(pv, lfc, se, alpha) {
  n <- nrow(pv)
  p_aw <- rep(NA_real_, n)
  w    <- matrix(NA_integer_, n, ncol(pv), dimnames = dimnames(pv))
  pattern <- apply(!is.na(pv), 1, function(x) paste(as.integer(x),
                                                     collapse = ""))
  for (pt in unique(pattern)) {
    rows <- which(pattern == pt)
    cols <- which(strsplit(pt, "")[[1]] == "1")
    if (length(cols) < 2) next
    aw <- AWFisher::AWFisher_pvalue(pv[rows, cols, drop = FALSE])
    p_aw[rows] <- aw$pvalues
    w[rows, cols] <- matrix(as.integer(aw$weights), nrow = length(rows))
  }

  # AW-REM: random effects (REML tau^2) restricted to selected studies
  sel <- !is.na(w) & w == 1
  lfc_s <- lfc; lfc_s[!sel] <- NA_real_
  se_s  <- se;  se_s[!sel]  <- NA_real_
  het <- .heterogeneity(lfc_s, se_s, "REML")
  re  <- .meta_re(lfc_s, se_s, het$tau2, hksj = FALSE, alpha = alpha)
  list(est = re$est, se = re$se, ci_lower = re$ci_lower,
       ci_upper = re$ci_upper, stat = rep(NA_real_, n), p = p_aw,
       weights = w)
}
