#' @include AllClasses.R
NULL

# ============================================================================
# Module 3: Per-Study Differential Expression
# ============================================================================

#' Run differential expression analysis on a single study
#'
#' Performs differential expression analysis on one
#' \code{metaXpressStudy} using DESeq2, edgeR, or limma-voom. Stores the
#' full result (all genes) in the \code{de_result} slot of the returned
#' study object.
#'
#' @param study A \code{\linkS4class{metaXpressStudy}} object with raw counts.
#' @param method Character scalar. DE method. One of \code{"DESeq2"}
#'   (default), \code{"edgeR"}, or \code{"limma-voom"}.
#' @param formula A formula specifying the design. Default: \code{~ condition}.
#'   The \code{condition} column in \code{study@metadata} is used as the
#'   contrast variable; the last level is the reference.
#' @param ... Additional arguments passed to the underlying DE function.
#'
#' @return The input \code{study} with the \code{de_result} slot populated.
#'   \code{de_result} is a \code{data.frame} with columns: \code{gene_id},
#'   \code{log2FC}, \code{pvalue}, \code{padj}, \code{baseMean}, \code{method}.
#'
#' @details
#' \strong{CRITICAL:} Always use \code{padj} (not raw \code{pvalue}) for
#' significance filtering downstream. The default significance threshold is
#' \code{padj <= 0.05} AND \code{|log2FC| >= 1}.
#'
#' @references
#' Ritchie et al. (2015) Nucleic Acids Research. \doi{10.1093/nar/gkv007}
#'
#' @examples
#' \dontrun{
#'   study <- mx_de(study, method = "DESeq2", formula = ~ condition)
#'   head(study@de_result)
#' }
#'
#' @seealso \code{\link{mx_de_all}}, \code{\link{mx_de_summary}}
#'
#' @importFrom methods is
#' @export
mx_de <- function(study,
                   method  = c("DESeq2", "edgeR", "limma-voom"),
                   formula = ~ condition,
                   ...) {
  if (!is(study, "metaXpressStudy"))
    stop("'study' must be a metaXpressStudy object")
  method <- match.arg(method)

  if (!("condition" %in% colnames(study@metadata)))
    stop("study@metadata must contain a 'condition' column")

  message("  Running ", method, " on study: ", study@accession)

  de_result <- switch(method,
    DESeq2      = .de_deseq2(study, formula = formula, ...),
    edgeR       = .de_edger(study, formula = formula, ...),
    `limma-voom` = .de_limma_voom(study, formula = formula, ...)
  )

  de_result$method <- method
  study@de_result  <- de_result
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
#' @param formula A formula specifying the design. Default: \code{~ condition}.
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
#' @importFrom BiocParallel bplapply SerialParam
#' @export
mx_de_all <- function(studies, method = "DESeq2", formula = ~ condition,
                       BPPARAM = BiocParallel::SerialParam(), ...) {
  if (!is.list(studies) || length(studies) == 0)
    stop("'studies' must be a non-empty list of metaXpressStudy objects")

  message("Running ", method, " on ", length(studies), " studies...")

  BiocParallel::bplapply(studies, function(s) {
    mx_de(s, method = method, formula = formula, ...)
  }, BPPARAM = BPPARAM)
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

.de_deseq2 <- function(study, formula, ...) {
  counts   <- round(study@counts)
  metadata <- study@metadata

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData   = metadata,
    design    = formula
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE, ...)

  cond_levels <- levels(factor(metadata$condition))
  contrast    <- c("condition", cond_levels[2], cond_levels[1])

  res <- DESeq2::results(dds, contrast = contrast,
                          independentFiltering = TRUE)
  res_df <- as.data.frame(res)

  data.frame(
    gene_id  = rownames(res_df),
    log2FC   = res_df$log2FoldChange,
    pvalue   = res_df$pvalue,
    padj     = res_df$padj,
    baseMean = res_df$baseMean,
    lfcSE    = res_df$lfcSE,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

.de_edger <- function(study, formula, ...) {
  counts   <- round(study@counts)
  metadata <- study@metadata
  cond     <- factor(metadata$condition)

  design <- stats::model.matrix(~ cond)
  dge    <- edgeR::DGEList(counts = counts, group = cond)
  dge    <- edgeR::calcNormFactors(dge)
  dge    <- edgeR::estimateDisp(dge, design = design)

  fit  <- edgeR::glmQLFit(dge, design = design)
  qlft <- edgeR::glmQLFTest(fit, coef = ncol(design))

  tt <- edgeR::topTags(qlft, n = Inf, sort.by = "none")$table

  data.frame(
    gene_id  = rownames(tt),
    log2FC   = tt$logFC,
    pvalue   = tt$PValue,
    padj     = tt$FDR,
    baseMean = rowMeans(counts),
    lfcSE    = NA_real_,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

.de_limma_voom <- function(study, formula, ...) {
  counts   <- study@counts
  metadata <- study@metadata
  cond     <- factor(metadata$condition)

  design <- stats::model.matrix(~ cond)
  dge    <- edgeR::DGEList(counts = round(counts))
  dge    <- edgeR::calcNormFactors(dge)
  v      <- limma::voom(dge, design = design)
  fit    <- limma::lmFit(v, design = design)
  fit    <- limma::eBayes(fit)

  tt <- limma::topTable(fit, coef = ncol(design), n = Inf, sort.by = "none")

  data.frame(
    gene_id  = rownames(tt),
    log2FC   = tt$logFC,
    pvalue   = tt$P.Value,
    padj     = tt$adj.P.Val,
    baseMean = rowMeans(counts),
    lfcSE    = tt$logFC / tt$t,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}
