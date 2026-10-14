#' @include AllClasses.R
NULL

# ============================================================================
# Module 3: Per-Study Differential Expression
# ============================================================================

#' Run differential expression analysis on a single study
#'
#' Performs differential expression analysis on one
#' \code{metaXpressStudy} using DESeq2, edgeR (quasi-likelihood), or
#' limma-voom. Stores the full result (all tested genes) in the
#' \code{de_result} slot of the returned study object, in a common format
#' that \code{\link{mx_meta}} consumes directly.
#'
#' @param study A \code{\linkS4class{metaXpressStudy}} object with raw counts.
#' @param method Character scalar. DE method. One of \code{"DESeq2"}
#'   (default), \code{"edgeR"}, or \code{"limma-voom"}.
#' @param formula A design formula whose variables are columns of
#'   \code{colData(study)}. Must contain \code{condition}; any other terms
#'   (e.g. \code{~ sex + age + condition}) are adjusted for. Default:
#'   \code{~ condition}.
#' @param contrast Optional length-2 character vector
#'   \code{c(case_level, control_level)} naming the two \code{condition}
#'   levels to compare. Fold changes are \code{case / control}. When
#'   \code{NULL}, the control level is auto-detected (e.g. \code{"control"},
#'   \code{"normal"}, \code{"healthy"}, \code{"untreated"}) and the first
#'   other level is the case. Set this explicitly for every study in a
#'   meta-analysis so fold-change directions agree.
#' @param shrink Character scalar. Optional log fold-change shrinkage for
#'   DESeq2: \code{"none"} (default), \code{"apeglm"}, \code{"ashr"} or
#'   \code{"normal"}. Shrunken values are added as \code{log2FC_shrunk} for
#'   ranking and plotting; the unshrunken maximum-likelihood \code{log2FC}
#'   and \code{lfcSE} are always kept because inverse-variance meta-analysis
#'   requires unbiased estimates with their sampling variance.
#' @param ... Additional arguments passed to \code{DESeq2::DESeq()},
#'   \code{edgeR::glmQLFit()} or \code{limma::lmFit()}.
#'
#' @return The input \code{study} with the \code{de_result} slot populated.
#'   \code{de_result} is a \code{data.frame} with columns:
#'   \describe{
#'     \item{gene_id}{Gene identifier (rownames of the count matrix).}
#'     \item{log2FC}{Log2 fold change, case vs control.}
#'     \item{lfcSE}{Standard error of \code{log2FC} (log2 scale).}
#'     \item{stat}{Signed standard-normal test statistic
#'       (\code{sign(log2FC) * qnorm(pvalue / 2, lower.tail = FALSE)}
#'       for edgeR/limma, the Wald statistic for DESeq2).}
#'     \item{pvalue, padj}{Raw and Benjamini-Hochberg adjusted p-values.}
#'     \item{baseMean}{Mean normalized expression.}
#'     \item{n_case, n_control}{Number of samples in each group, used by
#'       \code{\link{mx_meta}} for sample-size weighting.}
#'     \item{method}{DE engine used.}
#'   }
#'   The case and control labels are stored in
#'   \code{attr(de_result, "contrast")}.
#'
#' @details
#' For edgeR and limma-voom, low-count genes are removed with
#' \code{edgeR::filterByExpr()} before testing; DESeq2 applies its own
#' independent filtering to \code{padj} but reports all genes. Genes removed
#' by filtering are treated as \emph{not measured} in that study and are
#' handled explicitly by Module 5 (\code{\link{mx_missing_summary}}).
#'
#' For edgeR, standard errors are back-calculated from the quasi-likelihood
#' test as \code{|log2FC| / |z|}, where \code{z} is the normal quantile of
#' the two-sided p-value. This places all three engines on the same
#' log2-scale effect/variance footing required for inverse-variance
#' meta-analysis (Rau et al. 2014).
#'
#' Use \code{padj} (not raw \code{pvalue}) for per-study significance calls.
#'
#' @references
#' Love, M.I., Huber, W. & Anders, S. (2014) Moderated estimation of fold
#' change and dispersion for RNA-seq data with DESeq2. \emph{Genome Biology},
#' \strong{15}(12), 550. \doi{10.1186/s13059-014-0550-8}
#'
#' Chen, Y., Lun, A.T.L. & Smyth, G.K. (2016) From reads to genes to
#' pathways: differential expression analysis of RNA-Seq experiments using
#' Rsubread and the edgeR quasi-likelihood pipeline. \emph{F1000Research},
#' \strong{5}, 1438. \doi{10.12688/f1000research.8987.2}
#'
#' Law, C.W. et al. (2014) voom: precision weights unlock linear model
#' analysis tools for RNA-seq read counts. \emph{Genome Biology},
#' \strong{15}(2), R29. \doi{10.1186/gb-2014-15-2-r29}
#'
#' Zhu, A., Ibrahim, J.G. & Love, M.I. (2019) Heavy-tailed prior
#' distributions for sequence count data: removing the noise and preserving
#' large differences. \emph{Bioinformatics}, \strong{35}(12), 2084--2092.
#' \doi{10.1093/bioinformatics/bty895}
#'
#' @examples
#' if (requireNamespace("airway", quietly = TRUE) &&
#'     requireNamespace("DESeq2", quietly = TRUE)) {
#'   # Himes et al. (2014) airway smooth muscle cells, GEO GSE52778:
#'   # dexamethasone vs untreated, paired by donor cell line.
#'   data("airway", package = "airway")
#'   study <- mx_from_se(airway, condition_col = "dex",
#'                       accession = "GSE52778")
#'   study <- study[rowSums(SummarizedExperiment::assay(study)) >= 10, ]
#'   study <- mx_de(study, formula = ~ cell + condition,
#'                  contrast = c("trt", "untrt"))
#'   de <- study@de_result
#'   head(de[order(de$padj), ])
#' }
#'
#' @seealso \code{\link{mx_de_all}}, \code{\link{mx_de_summary}},
#'   \code{\link{mx_meta}}
#'
#' @importFrom methods is
#' @export
mx_de <- function(study,
                  method   = c("DESeq2", "edgeR", "limma-voom"),
                  formula  = ~ condition,
                  contrast = NULL,
                  shrink   = c("none", "apeglm", "ashr", "normal"),
                  ...) {
  if (!is(study, "metaXpressStudy"))
    stop("'study' must be a metaXpressStudy object")
  method <- match.arg(method)
  shrink <- match.arg(shrink)

  samples <- .samples(study)
  if (!("condition" %in% colnames(samples)))
    stop("colData(study) must contain a 'condition' column")
  if (!inherits(formula, "formula") || !"condition" %in% all.vars(formula))
    stop("'formula' must be a formula containing 'condition'")
  missing_vars <- setdiff(all.vars(formula), colnames(samples))
  if (length(missing_vars) > 0)
    stop("Formula variable(s) not found in colData(study): ",
         paste(missing_vars, collapse = ", "))
  if (shrink != "none" && method != "DESeq2")
    stop("'shrink' is only available for method = 'DESeq2'")

  design_info <- .prepare_design(samples, formula, contrast)
  if (any(design_info$n < 2))
    stop("Each compared condition needs >= 2 samples (got ",
         paste(names(design_info$n), design_info$n, sep = "=", collapse = ", "),
         ")")

  message("  Running ", method, " on study: ", study@accession, " (",
          design_info$case, " vs ", design_info$control, ")")

  de_result <- switch(method,
    DESeq2       = .de_deseq2(study, design_info, shrink = shrink, ...),
    edgeR        = .de_edger(study, design_info, ...),
    `limma-voom` = .de_limma_voom(study, design_info, ...)
  )

  de_result$n_case    <- unname(design_info$n[design_info$case])
  de_result$n_control <- unname(design_info$n[design_info$control])
  de_result$method    <- method
  attr(de_result, "contrast") <- c(case    = design_info$case,
                                   control = design_info$control)
  study@de_result <- de_result
  study
}

#' Run differential expression on all studies
#'
#' Applies \code{\link{mx_de}} to each study in a list, optionally in
#' parallel using \pkg{BiocParallel}.
#'
#' @param studies A named list of \code{\linkS4class{metaXpressStudy}} objects.
#' @param method Character scalar. DE method passed to \code{\link{mx_de}}.
#'   Default: \code{"DESeq2"}.
#' @param formula A design formula (see \code{\link{mx_de}}), or a list of
#'   formulas with one entry per study when covariates differ between studies.
#'   Default: \code{~ condition}.
#' @param contrast \code{NULL}, a length-2 character vector
#'   \code{c(case, control)} applied to every study, or a list with one such
#'   vector per study (labels often differ between GEO series).
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object
#'   controlling parallelization. Defaults to
#'   \code{\link[BiocParallel]{SerialParam}()} (single-core). Use
#'   \code{\link[BiocParallel]{MulticoreParam}()} for multi-core execution.
#' @param ... Additional arguments passed to \code{\link{mx_de}}.
#'
#' @return The input list of studies, each with \code{de_result} slot filled.
#'
#' @examples
#' \dontrun{
#'   studies <- mx_de_all(studies, method = "DESeq2")
#'   # Parallel:
#'   studies <- mx_de_all(studies, method = "DESeq2",
#'                         BPPARAM = BiocParallel::MulticoreParam(4))
#' }
#'
#' @seealso \code{\link{mx_de}}, \code{\link{mx_de_summary}}
#'
#' @importFrom BiocParallel bpmapply SerialParam
#' @export
mx_de_all <- function(studies, method = "DESeq2", formula = ~ condition,
                      contrast = NULL,
                      BPPARAM = BiocParallel::SerialParam(), ...) {
  if (!is.list(studies) || length(studies) == 0)
    stop("'studies' must be a non-empty list of metaXpressStudy objects")

  message("Running ", method, " on ", length(studies), " studies...")

  if (!is.list(formula)) formula <- rep(list(formula), length(studies))
  if (!is.list(contrast)) contrast <- rep(list(contrast), length(studies))
  if (length(formula) != length(studies) ||
      length(contrast) != length(studies))
    stop("'formula' and 'contrast' must be single values or lists with one ",
         "entry per study")

  out <- BiocParallel::bpmapply(function(s, f, ct) {
    mx_de(s, method = method, formula = f, contrast = ct, ...)
  }, studies, formula, contrast, SIMPLIFY = FALSE, BPPARAM = BPPARAM)
  names(out) <- names(studies)
  out
}

#' Summarise per-study DE results
#'
#' Returns a summary table showing the number of significant DEGs in each
#' study at the default threshold (\code{padj <= 0.05}, \code{|log2FC| >= 1}).
#'
#' @param studies A named list of \code{\linkS4class{metaXpressStudy}} objects
#'   after running \code{\link{mx_de_all}}.
#' @param padj_threshold Numeric. Adjusted p-value threshold. Default: 0.05.
#' @param lfc_threshold Numeric. Absolute log2 fold-change threshold.
#'   Default: 1.
#'
#' @return A \code{data.frame} with columns \code{study}, \code{n_total},
#'   \code{n_up}, \code{n_down}, \code{method}.
#'
#' @examples
#' \dontrun{
#'   summary_tbl <- mx_de_summary(studies)
#' }
#'
#' @seealso \code{\link{mx_de_all}}
#'
#' @export
mx_de_summary <- function(studies, padj_threshold = 0.05, lfc_threshold = 1) {
  if (!is.list(studies))
    stop("'studies' must be a list of metaXpressStudy objects")

  rows <- lapply(names(studies), function(nm) {
    de <- studies[[nm]]@de_result
    if (nrow(de) == 0) {
      warning("Study ", nm, " has no DE results; run mx_de_all() first.")
      return(data.frame(study = nm, n_total = NA, n_up = NA, n_down = NA,
                        method = NA, stringsAsFactors = FALSE))
    }
    sig <- de[de$padj <= padj_threshold &
                abs(de$log2FC) >= lfc_threshold &
                !is.na(de$padj), ]
    data.frame(
      study   = nm,
      n_total = nrow(sig),
      n_up    = sum(sig$log2FC > 0),
      n_down  = sum(sig$log2FC < 0),
      method  = if ("method" %in% colnames(de)) de$method[1] else NA,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

# ============================================================================
# Internal DE method implementations
# ============================================================================

# Resolve case/control levels, relevel condition so control is the reference,
# and build the model matrix shared by edgeR and limma.
.prepare_design <- function(samples, formula, contrast) {
  cond <- if (is.null(contrast)) .make_condition_factor(samples$condition)
          else factor(samples$condition)
  if (is.null(contrast)) {
    control <- levels(cond)[1]
    case    <- levels(cond)[2]
    if (nlevels(cond) > 2)
      warning("'condition' has ", nlevels(cond), " levels; comparing '",
              case, "' vs '", control, "'. Pass 'contrast' to choose.")
  } else {
    if (!is.character(contrast) || length(contrast) != 2)
      stop("'contrast' must be c(case_level, control_level)")
    bad <- setdiff(contrast, levels(cond))
    if (length(bad) > 0)
      stop("contrast level(s) not found in condition: ",
           paste(bad, collapse = ", "), ". Available: ",
           paste(levels(cond), collapse = ", "))
    case    <- contrast[1]
    control <- contrast[2]
  }
  if (is.na(case)) stop("'condition' must have at least 2 levels")

  samples$condition <- stats::relevel(factor(cond), ref = control)
  for (v in setdiff(all.vars(formula), "condition"))
    if (is.character(samples[[v]])) samples[[v]] <- factor(samples[[v]])

  design <- stats::model.matrix(formula, data = samples)
  coef   <- paste0("condition", make.names(case))
  colnames(design) <- make.names(colnames(design))
  if (!coef %in% colnames(design))
    stop("Could not find coefficient '", coef, "' in the design matrix")
  if (qr(design)$rank < ncol(design))
    stop("Design matrix is not full rank; check that covariates in 'formula' ",
         "are not confounded with condition")

  n <- table(samples$condition)[c(case, control)]
  list(samples = samples, formula = formula, design = design, coef = coef,
       case = case,
       control = control, n = stats::setNames(as.integer(n), c(case, control)))
}

# Signed z-statistic from a two-sided p-value and the effect direction
.p_to_z <- function(p, lfc) {
  sign(lfc) * stats::qnorm(pmax(p, .Machine$double.xmin) / 2,
                           lower.tail = FALSE)
}

.de_deseq2 <- function(study, d, shrink = "none", ...) {
  .require_package("DESeq2", "DESeq2 differential expression")
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(.counts(study)),
    colData   = d$samples,
    design    = d$formula
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE, ...)

  res <- DESeq2::results(dds, contrast = c("condition", d$case, d$control),
                         independentFiltering = TRUE)
  out <- data.frame(
    gene_id  = rownames(res),
    log2FC   = res$log2FoldChange,
    lfcSE    = res$lfcSE,
    stat     = res$stat,
    pvalue   = res$pvalue,
    padj     = res$padj,
    baseMean = res$baseMean,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  if (shrink != "none") {
    if (shrink == "apeglm") .require_package("apeglm", "apeglm shrinkage")
    if (shrink == "ashr")   .require_package("ashr", "ashr shrinkage")
    coef_name <- paste0("condition_", make.names(d$case), "_vs_",
                        make.names(d$control))
    shr <- if (shrink == "ashr") {
      DESeq2::lfcShrink(dds, contrast = c("condition", d$case, d$control),
                        res = res, type = "ashr", quiet = TRUE)
    } else {
      DESeq2::lfcShrink(dds, coef = coef_name, type = shrink, quiet = TRUE)
    }
    out$log2FC_shrunk <- shr$log2FoldChange
  }
  out
}

.de_edger <- function(study, d, ...) {
  .require_package("edgeR", "edgeR differential expression")
  counts <- round(.counts(study))
  dge    <- edgeR::DGEList(counts = counts, group = d$samples$condition)

  # Filter low-expression genes (standard edgeR practice)
  keep <- edgeR::filterByExpr(dge, design = d$design)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- edgeR::calcNormFactors(dge)
  dge  <- edgeR::estimateDisp(dge, design = d$design)

  fit  <- edgeR::glmQLFit(dge, design = d$design, ...)
  qlft <- edgeR::glmQLFTest(fit, coef = d$coef)
  tt   <- edgeR::topTags(qlft, n = Inf, sort.by = "none")$table

  # GLM coefficients are on the natural-log scale, so the QL test is
  # converted to a z-score and the log2-scale SE is back-calculated.
  z <- .p_to_z(tt$PValue, tt$logFC)
  se <- abs(tt$logFC) / abs(z)
  se[!is.finite(se) | se == 0] <- NA_real_

  data.frame(
    gene_id  = rownames(tt),
    log2FC   = tt$logFC,
    lfcSE    = se,
    stat     = z,
    pvalue   = tt$PValue,
    padj     = tt$FDR,
    baseMean = rowMeans(edgeR::cpm(dge, log = FALSE)) *
      mean(dge$samples$lib.size) / 1e6,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

.de_limma_voom <- function(study, d, ...) {
  .require_package("limma", "limma-voom differential expression")
  .require_package("edgeR", "limma-voom differential expression")
  dge  <- edgeR::DGEList(counts = round(.counts(study)))

  # Filter low-expression genes (standard limma-voom practice)
  keep <- edgeR::filterByExpr(dge, design = d$design,
                              group = d$samples$condition)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- edgeR::calcNormFactors(dge)

  v   <- limma::voom(dge, design = d$design)
  fit <- limma::lmFit(v, design = d$design, ...)
  fit <- limma::eBayes(fit)
  tt  <- limma::topTable(fit, coef = d$coef, number = Inf, sort.by = "none")

  # Moderated SE of the log2FC: sqrt(s2.post) * stdev.unscaled
  se <- sqrt(fit$s2.post) * fit$stdev.unscaled[, d$coef]

  data.frame(
    gene_id  = rownames(tt),
    log2FC   = tt$logFC,
    lfcSE    = unname(se[rownames(tt)]),
    stat     = .p_to_z(tt$P.Value, tt$logFC),
    pvalue   = tt$P.Value,
    padj     = tt$adj.P.Val,
    baseMean = 2^tt$AveExpr * mean(dge$samples$lib.size) / 1e6,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}
