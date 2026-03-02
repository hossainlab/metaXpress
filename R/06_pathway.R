#' @include AllClasses.R
NULL

# ============================================================================
# Module 6: Pathway Meta-Analysis
# ============================================================================

#' Perform pathway meta-analysis on meta-analysis results
#'
#' Runs over-representation analysis (ORA) or gene set enrichment analysis
#' (GSEA) using the meta-analysis gene rankings, then combines pathway results
#' across studies.
#'
#' @param meta_result A \code{\linkS4class{metaXpressResult}} object from
#'   \code{\link{mx_meta}}.
#' @param db Character scalar. Gene set database. One of \code{"Hallmarks"}
#'   (default), \code{"KEGG"}, \code{"Reactome"}, or \code{"GO_BP"}.
#' @param method Character scalar. Enrichment method. One of \code{"ORA"}
#'   (default) or \code{"GSEA"}.
#' @param padj_threshold Numeric. Adjusted p-value threshold for defining the
#'   significant gene list (ORA only). Default: \code{0.05}.
#' @param lfc_threshold Numeric. Log2 fold-change threshold (ORA only).
#'   Default: \code{1}.
#' @param ... Additional arguments passed to \pkg{clusterProfiler}.
#'
#' @return A \code{\linkS4class{metaXpressResult}} object with the
#'   \code{pathway_result} slot populated. The \code{pathway_result}
#'   \code{data.frame} contains: \code{pathway_id}, \code{pathway_name},
#'   \code{n_genes}, \code{pvalue}, \code{padj}, \code{gene_ratio}.
#'
#' @references
#' Zeng et al. (2018) Genes. \doi{10.1101/444604}
#'
#' @examples
#' \dontrun{
#'   result <- mx_pathway_meta(meta_result, db = "Hallmarks")
#'   head(result@pathway_result[order(result@pathway_result$padj), ])
#' }
#'
#' @seealso \code{\link{mx_pathway_consensus}}, \code{\link{mx_pathway_dedup}},
#'   \code{\link{mx_pathway_heatmap}}
#'
#' @importFrom clusterProfiler enricher GSEA
#' @importFrom msigdbr msigdbr
#' @importFrom methods is
#' @export
mx_pathway_meta <- function(meta_result, db = c("Hallmarks", "KEGG",
                                                  "Reactome", "GO_BP"),
                              method = c("ORA", "GSEA"),
                              padj_threshold = 0.05,
                              lfc_threshold  = 1,
                              ...) {
  if (!is(meta_result, "metaXpressResult"))
    stop("'meta_result' must be a metaXpressResult object")
  db     <- match.arg(db)
  method <- match.arg(method)

  mt <- meta_result@meta_table

  gene_sets <- .get_gene_sets(db)

  pathway_result <- switch(method,
    ORA  = .pathway_ora(mt, gene_sets, padj_threshold, lfc_threshold),
    GSEA = .pathway_gsea(mt, gene_sets)
  )

  meta_result@pathway_result <- pathway_result
  meta_result
}

#' Find consensus pathways across studies
#'
#' Identifies pathways that are significantly enriched in the majority of
#' individual studies, indicating robust cross-study signal.
#'
#' @param pathway_results A list of \code{data.frame} objects from per-study
#'   enrichment analyses.
#' @param min_fraction Numeric. Minimum fraction of studies in which a pathway
#'   must be significant. Default: \code{0.5}.
#' @param padj_threshold Numeric. Significance threshold per study.
#'   Default: \code{0.05}.
#'
#' @return A \code{data.frame} of consensus pathways with columns
#'   \code{pathway_id}, \code{pathway_name}, \code{n_significant_studies},
#'   \code{fraction_studies}.
#'
#' @examples
#' \dontrun{
#'   consensus <- mx_pathway_consensus(pathway_results)
#' }
#'
#' @seealso \code{\link{mx_pathway_meta}}, \code{\link{mx_pathway_dedup}}
#'
#' @export
mx_pathway_consensus <- function(pathway_results, min_fraction = 0.5,
                                  padj_threshold = 0.05) {
  if (!is.list(pathway_results))
    stop("'pathway_results' must be a list of pathway enrichment data.frames")

  k          <- length(pathway_results)
  all_paths  <- unique(unlist(lapply(pathway_results,
                                      function(r) r$pathway_id)))

  n_sig <- vapply(all_paths, function(pid) {
    sum(vapply(pathway_results, function(r) {
      row <- r[r$pathway_id == pid, ]
      if (nrow(row) == 0) return(FALSE)
      !is.na(row$padj[1]) && row$padj[1] <= padj_threshold
    }, logical(1)))
  }, integer(1))

  fraction <- n_sig / k
  keep     <- fraction >= min_fraction

  data.frame(
    pathway_id          = all_paths[keep],
    n_significant_studies = n_sig[keep],
    fraction_studies    = fraction[keep],
    stringsAsFactors    = FALSE,
    row.names           = NULL
  )
}

#' Remove redundant pathways
#'
#' Reduces pathway result sets by removing pathways with high gene set overlap,
#' retaining the most significant representative from each cluster of similar
#' pathways.
#'
#' @param pathway_results A \code{data.frame} from \code{\link{mx_pathway_meta}}
#'   or a list thereof.
#' @param jaccard_threshold Numeric. Jaccard similarity threshold above which
#'   two pathways are considered redundant. Default: \code{0.5}.
#'
#' @return A \code{data.frame} with redundant pathways removed.
#'
#' @references
#' Zeng et al. (2018) Genes. \doi{10.1101/444604}
#'
#' @examples
#' \dontrun{
#'   result <- mx_pathway_dedup(result@pathway_result)
#' }
#'
#' @export
mx_pathway_dedup <- function(pathway_results, jaccard_threshold = 0.5) {
  stop("mx_pathway_dedup() is not yet implemented.")
}

#' Cross-study pathway heatmap
#'
#' Creates a heatmap showing enrichment scores or -log10(padj) for top
#' pathways across all studies.
#'
#' @param pathway_results A list of per-study pathway enrichment
#'   \code{data.frame} objects, or a single \code{data.frame} from
#'   \code{\link{mx_pathway_meta}}.
#' @param top_n Integer. Number of top pathways to display. Default: \code{30}.
#' @param value Character scalar. Value to display: \code{"padj"} (default)
#'   or \code{"NES"} (for GSEA results).
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' \dontrun{
#'   mx_pathway_heatmap(result@pathway_result, top_n = 20)
#' }
#'
#' @seealso \code{\link{mx_pathway_meta}}
#'
#' @export
mx_pathway_heatmap <- function(pathway_results, top_n = 30,
                                value = c("padj", "NES")) {
  stop("mx_pathway_heatmap() is not yet implemented.")
}

# ============================================================================
# Internal helpers
# ============================================================================

.get_gene_sets <- function(db) {
  category_map <- list(
    Hallmarks = list(category = "H",  subcategory = NULL),
    KEGG      = list(category = "C2", subcategory = "CP:KEGG"),
    Reactome  = list(category = "C2", subcategory = "CP:REACTOME"),
    GO_BP     = list(category = "C5", subcategory = "GO:BP")
  )

  spec <- category_map[[db]]
  gs   <- msigdbr::msigdbr(species = "Homo sapiens",
                             category   = spec$category,
                             subcategory = spec$subcategory)

  split(gs$gene_symbol, gs$gs_name)
}

.pathway_ora <- function(meta_table, gene_sets, padj_threshold, lfc_threshold) {
  sig_genes <- meta_table$gene_id[
    !is.na(meta_table$meta_padj) &
      meta_table$meta_padj  <= padj_threshold &
      abs(meta_table$meta_log2FC) >= lfc_threshold
  ]

  universe <- meta_table$gene_id

  if (length(sig_genes) == 0) {
    warning("No significant genes found at padj <= ", padj_threshold,
            " and |log2FC| >= ", lfc_threshold, ". ",
            "Returning empty pathway result.")
    return(data.frame())
  }

  gs_df <- data.frame(
    gs_name     = rep(names(gene_sets), lengths(gene_sets)),
    gene_symbol = unlist(gene_sets, use.names = FALSE),
    stringsAsFactors = FALSE
  )

  res <- clusterProfiler::enricher(
    gene         = sig_genes,
    universe     = universe,
    TERM2GENE    = gs_df,
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )

  if (is.null(res) || nrow(as.data.frame(res)) == 0)
    return(data.frame())

  res_df <- as.data.frame(res)
  data.frame(
    pathway_id   = res_df$ID,
    pathway_name = res_df$Description,
    n_genes      = res_df$Count,
    pvalue       = res_df$pvalue,
    padj         = res_df$p.adjust,
    gene_ratio   = res_df$GeneRatio,
    stringsAsFactors = FALSE
  )
}

.pathway_gsea <- function(meta_table, gene_sets) {
  ranked <- stats::setNames(meta_table$meta_log2FC, meta_table$gene_id)
  ranked <- sort(ranked[!is.na(ranked)], decreasing = TRUE)

  gs_df <- data.frame(
    gs_name     = rep(names(gene_sets), lengths(gene_sets)),
    gene_symbol = unlist(gene_sets, use.names = FALSE),
    stringsAsFactors = FALSE
  )

  res <- tryCatch(
    clusterProfiler::GSEA(
      geneList     = ranked,
      TERM2GENE    = gs_df,
      pvalueCutoff = 1,
      verbose      = FALSE
    ),
    error = function(e) {
      warning("GSEA failed: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(res)) return(data.frame())
  res_df <- as.data.frame(res)
  data.frame(
    pathway_id   = res_df$ID,
    pathway_name = res_df$Description,
    n_genes      = res_df$setSize,
    pvalue       = res_df$pvalue,
    padj         = res_df$p.adjust,
    NES          = res_df$NES,
    stringsAsFactors = FALSE
  )
}
