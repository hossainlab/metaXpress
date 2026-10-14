#' @import methods
NULL

# ============================================================================
# metaXpressStudy — single RNA-seq study
# ============================================================================

#' The metaXpressStudy class
#'
#' Represents a single bulk RNA-seq study. \code{metaXpressStudy} extends
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}}, so raw counts live
#' in the \code{"counts"} assay, sample metadata in \code{colData()}, and
#' optional gene annotation in \code{rowData()}. All standard
#' SummarizedExperiment accessors and subsetting (\code{study[genes, samples]},
#' \code{study$condition}) work directly. Populated by
#' \code{\link{mx_fetch_geo}}, \code{\link{mx_fetch_sra}},
#' \code{\link{mx_load_local}}, or \code{\link{mx_from_se}}.
#'
#' The \code{colData} must contain a \code{condition} column (case/control
#' labels) once an accession is set.
#'
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
#' @seealso \code{\link{metaXpressStudy}}, \code{\link{mx_from_se}},
#'   \code{\link{mx_fetch_geo}}, \code{\link{mx_qc_study}}, \code{\link{mx_de}}
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @exportClass metaXpressStudy
setClass(
  "metaXpressStudy",
  contains = "SummarizedExperiment",
  representation(
    accession = "character",
    organism  = "character",
    qc_score  = "numeric",
    de_result = "data.frame"
  ),
  prototype(
    accession = character(0),
    organism  = "Homo sapiens",
    qc_score  = NA_real_,
    de_result = data.frame()
  ),
  validity = function(object) {
    msgs <- character(0)

    if (length(object@accession) > 0) {
      if (!("condition" %in% colnames(SummarizedExperiment::colData(object))))
        msgs <- c(msgs, "colData must contain a 'condition' column")

      if (!is.na(object@qc_score) &&
          (object@qc_score < 0 || object@qc_score > 10))
        msgs <- c(msgs, "qc_score must be between 0 and 10")
    }

    if (length(msgs)) msgs else TRUE
  }
)

#' Construct a metaXpressStudy
#'
#' Builds a \code{\linkS4class{metaXpressStudy}} from a count matrix and a
#' sample metadata table.
#'
#' @param counts A numeric matrix of raw counts (genes x samples). Rownames
#'   must be gene identifiers; colnames must be sample identifiers.
#' @param metadata A \code{data.frame} of sample metadata, one row per column
#'   of \code{counts}. Must contain a \code{condition} column when
#'   \code{accession} is set.
#' @param accession A length-1 character accession ID, or \code{"local"}.
#' @param organism A length-1 character species name.
#' @param qc_score A length-1 numeric QC score, \code{NA} if not yet run.
#' @param de_result A \code{data.frame} of DE results (usually empty).
#' @param gene_info Optional \code{data.frame} or \code{DataFrame} of gene
#'   annotation, one row per row of \code{counts}; stored as \code{rowData}.
#'
#' @return A \code{\linkS4class{metaXpressStudy}} object.
#'
#' @examples
#' counts <- matrix(rpois(24, 100), nrow = 4,
#'                  dimnames = list(paste0("GENE", 1:4), paste0("S", 1:6)))
#' meta <- data.frame(sample_id = paste0("S", 1:6),
#'                    condition = rep(c("control", "case"), each = 3))
#' study <- metaXpressStudy(counts, meta, accession = "local")
#' study
#'
#' @export
metaXpressStudy <- function(counts    = matrix(0L, nrow = 0, ncol = 0),
                            metadata  = data.frame(),
                            accession = character(0),
                            organism  = "Homo sapiens",
                            qc_score  = NA_real_,
                            de_result = data.frame(),
                            gene_info = NULL) {
  counts <- as.matrix(counts)
  col_data <- if (nrow(metadata) == 0 && ncol(counts) > 0) {
    S4Vectors::DataFrame(row.names = colnames(counts))
  } else {
    metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
    rownames(metadata) <- colnames(counts)
    S4Vectors::DataFrame(metadata, check.names = FALSE)
  }

  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = counts),
    colData = col_data,
    rowData = gene_info
  )

  new("metaXpressStudy", se,
      accession = accession,
      organism  = organism,
      qc_score  = qc_score,
      de_result = de_result)
}

#' Convert a SummarizedExperiment to a metaXpressStudy
#'
#' Wraps any \code{SummarizedExperiment} (including \code{DESeqDataSet} and
#' \code{RangedSummarizedExperiment} objects, e.g. from recount3 or tximeta)
#' as a \code{\linkS4class{metaXpressStudy}} so it can enter the metaXpress
#' pipeline. \code{as(se, "metaXpressStudy")} calls this with defaults.
#'
#' @param se A \code{SummarizedExperiment} or subclass.
#' @param assay_name Name of the assay holding raw counts. Defaults to
#'   \code{"counts"}, falling back to the first assay if absent.
#' @param condition_col Name of the \code{colData} column holding case/control
#'   labels; copied to \code{condition} if different.
#' @param accession Accession ID. Defaults to \code{metadata(se)$accession}
#'   if present, otherwise \code{"local"}.
#' @param organism Species name. Defaults to \code{metadata(se)$organism} if
#'   present, otherwise \code{"Homo sapiens"}.
#'
#' @return A \code{\linkS4class{metaXpressStudy}} object.
#'
#' @examples
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays  = list(counts = matrix(rpois(24, 100), nrow = 4,
#'     dimnames = list(paste0("GENE", 1:4), paste0("S", 1:6)))),
#'   colData = data.frame(dex = rep(c("untrt", "trt"), each = 3))
#' )
#' study <- mx_from_se(se, condition_col = "dex")
#'
#' @seealso \code{\link{metaXpressStudy}}
#' @export
mx_from_se <- function(se,
                       assay_name    = "counts",
                       condition_col = "condition",
                       accession     = NULL,
                       organism      = NULL) {
  if (!is(se, "SummarizedExperiment"))
    stop("'se' must be a SummarizedExperiment object")

  if (!assay_name %in% SummarizedExperiment::assayNames(se)) {
    if (length(SummarizedExperiment::assays(se)) == 0)
      stop("'se' has no assays")
    message("Assay '", assay_name, "' not found; using the first assay.")
    assay_name <- 1L
  }

  col_data <- as.data.frame(SummarizedExperiment::colData(se))
  if (!condition_col %in% colnames(col_data))
    stop("colData(se) has no column '", condition_col, "'. ",
         "Set 'condition_col' to the column holding case/control labels.")
  col_data$condition <- as.character(col_data[[condition_col]])

  meta <- S4Vectors::metadata(se)
  if (is.null(accession))
    accession <- if (!is.null(meta$accession)) meta$accession else "local"
  if (is.null(organism))
    organism <- if (!is.null(meta$organism)) meta$organism else "Homo sapiens"

  metaXpressStudy(
    counts    = as.matrix(SummarizedExperiment::assay(se, assay_name)),
    metadata  = col_data,
    accession = accession,
    organism  = organism,
    gene_info = SummarizedExperiment::rowData(se)
  )
}

setAs("SummarizedExperiment", "metaXpressStudy", function(from) {
  mx_from_se(from)
})

setAs("metaXpressStudy", "SummarizedExperiment", function(from) {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays   = SummarizedExperiment::assays(from),
    colData  = SummarizedExperiment::colData(from),
    rowData  = SummarizedExperiment::rowData(from),
    metadata = S4Vectors::metadata(from)
  )
  S4Vectors::metadata(se)[c("accession", "organism", "qc_score")] <-
    list(from@accession, from@organism, from@qc_score)
  if (nrow(from@de_result) > 0)
    S4Vectors::metadata(se)$de_result <- from@de_result
  se
})

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
#'   One of \code{"random_effects"}, \code{"fixed_effects"}, \code{"fisher"},
#'   \code{"stouffer"}, \code{"inverse_normal"}, \code{"minp"},
#'   \code{"maxp"}, \code{"rop"}.
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
