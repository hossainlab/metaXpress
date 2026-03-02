#' @import methods
NULL

# ============================================================================
# metaXpressStudy — single RNA-seq study
# ============================================================================

#' The metaXpressStudy class
#'
#' Represents a single bulk RNA-seq study, holding raw counts, sample metadata,
#' accession information, and QC results. Populated by \code{\link{mx_fetch_geo}},
#' \code{\link{mx_fetch_sra}}, or \code{\link{mx_load_local}}.
#'
#' @slot counts A numeric matrix of raw counts (genes x samples). Rownames must
#'   be gene identifiers; colnames must be sample identifiers.
#' @slot metadata A \code{data.frame} of sample metadata. Must contain columns
#'   \code{condition} (case/control labels) and \code{sample_id}.
#' @slot accession A length-1 character giving the GEO/SRA accession ID
#'   (e.g., \code{"GSE12345"}) or \code{"local"} for user-supplied data.
#' @slot organism A length-1 character giving the species
#'   (e.g., \code{"Homo sapiens"}).
#' @slot qc_score A length-1 numeric in [0, 10]. Filled by
#'   \code{\link{mx_qc_study}}; \code{NA} until then.
#' @slot de_result A \code{data.frame} of per-study differential expression
#'   results. Empty until \code{\link{mx_de}} is called. Required columns when
#'   populated: \code{gene_id}, \code{log2FC}, \code{pvalue}, \code{padj},
#'   \code{baseMean}.
#'
#' @seealso \code{\link{mx_fetch_geo}}, \code{\link{mx_qc_study}},
#'   \code{\link{mx_de}}
#'
#' @exportClass metaXpressStudy
setClass(
  "metaXpressStudy",
  representation(
    counts    = "matrix",
    metadata  = "data.frame",
    accession = "character",
    organism  = "character",
    qc_score  = "numeric",
    de_result = "data.frame"
  ),
  prototype(
    counts    = matrix(0L, nrow = 0, ncol = 0),
    metadata  = data.frame(),
    accession = character(0),
    organism  = "Homo sapiens",
    qc_score  = NA_real_,
    de_result = data.frame()
  ),
  validity = function(object) {
    msgs <- character(0)

    if (length(object@accession) > 0) {
      if (!("condition" %in% colnames(object@metadata)))
        msgs <- c(msgs, "metadata must contain a 'condition' column")

      if (nrow(object@metadata) > 0 &&
          ncol(object@counts) != nrow(object@metadata))
        msgs <- c(msgs,
                  "ncol(counts) must equal nrow(metadata)")

      if (!is.na(object@qc_score) &&
          (object@qc_score < 0 || object@qc_score > 10))
        msgs <- c(msgs, "qc_score must be between 0 and 10")
    }

    if (length(msgs)) msgs else TRUE
  }
)

# ============================================================================
# metaXpressResult — meta-analysis output
# ============================================================================

#' The metaXpressResult class
#'
#' Holds the output of a meta-analysis performed by \code{\link{mx_meta}},
#' including gene-level statistics, heterogeneity estimates, and (optionally)
#' pathway enrichment results.
#'
#' @slot meta_table A \code{data.frame} with one row per gene. Required
#'   columns: \code{gene_id}, \code{meta_log2FC}, \code{meta_pvalue},
#'   \code{meta_padj}, \code{i_squared}, \code{q_stat}, \code{n_studies},
#'   \code{direction_consistency}.
#' @slot method A length-1 character naming the meta-analysis method used.
#'   One of \code{"fisher"}, \code{"stouffer"}, \code{"inverse_normal"},
#'   \code{"fixed_effects"}, \code{"random_effects"}, \code{"awmeta"}.
#' @slot n_studies A length-1 integer giving the total number of studies
#'   integrated.
#' @slot heterogeneity A \code{data.frame} with per-gene heterogeneity
#'   statistics: \code{gene_id}, \code{Q}, \code{df}, \code{I_sq},
#'   \code{tau_sq}, \code{p_heterogeneity}.
#' @slot pathway_result A \code{data.frame} of pathway meta-analysis results.
#'   Empty until \code{\link{mx_pathway_meta}} is called.
#'
#' @seealso \code{\link{mx_meta}}, \code{\link{mx_heterogeneity}},
#'   \code{\link{mx_pathway_meta}}
#'
#' @exportClass metaXpressResult
setClass(
  "metaXpressResult",
  representation(
    meta_table     = "data.frame",
    method         = "character",
    n_studies      = "integer",
    heterogeneity  = "data.frame",
    pathway_result = "data.frame"
  ),
  prototype(
    meta_table     = data.frame(),
    method         = character(0),
    n_studies      = 0L,
    heterogeneity  = data.frame(),
    pathway_result = data.frame()
  ),
  validity = function(object) {
    msgs <- character(0)

    if (nrow(object@meta_table) > 0) {
      required_cols <- c("gene_id", "meta_log2FC", "meta_pvalue",
                         "meta_padj", "i_squared", "n_studies")
      missing_cols <- setdiff(required_cols, colnames(object@meta_table))
      if (length(missing_cols))
        msgs <- c(msgs, paste("meta_table missing columns:",
                              paste(missing_cols, collapse = ", ")))
    }

    if (length(msgs)) msgs else TRUE
  }
)
