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
#' Wu, T. et al. (2021) clusterProfiler 4.0: a universal enrichment tool for
#' interpreting omics data. \emph{The Innovation}, \strong{2}(3), 100141.
#' \doi{10.1016/j.xinn.2021.100141}
#'
#' Liberzon, A. et al. (2015) The Molecular Signatures Database (MSigDB)
#' hallmark gene set collection. \emph{Cell Systems}, \strong{1}(6),
#' 417--425. \doi{10.1016/j.cels.2015.12.004}
#'
#' Subramanian, A. et al. (2005) Gene set enrichment analysis: a
#' knowledge-based approach for interpreting genome-wide expression profiles.
#' \emph{Proceedings of the National Academy of Sciences}, \strong{102}(43),
#' 15545--15550. \doi{10.1073/pnas.0506580102}
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
#' Gu, Z. & Hübschmann, D. (2023) simplifyEnrichment: a Bioconductor package
#' for clustering and visualizing functional enrichment results.
#' \emph{Genomics, Proteomics & Bioinformatics}, \strong{21}(1), 190--202.
#' \doi{10.1016/j.gpb.2022.04.008}
#'
#' @examples
#' \dontrun{
#'   result <- mx_pathway_dedup(result@pathway_result)
#' }
#'
#' @export
mx_pathway_dedup <- function(pathway_results, jaccard_threshold = 0.5) {
  if (is.list(pathway_results) && !is.data.frame(pathway_results)) {
    return(lapply(pathway_results, mx_pathway_dedup, jaccard_threshold = jaccard_threshold))
  }
  
  if (!is.data.frame(pathway_results) || nrow(pathway_results) == 0) {
    return(pathway_results)
  }
  
  if (!"geneID" %in% colnames(pathway_results)) {
    warning("No 'geneID' column found in pathway results. Cannot compute Jaccard similarity. Returning original.")
    return(pathway_results)
  }
  
  # Parse gene lists
  gene_lists <- strsplit(pathway_results$geneID, "/")
  n <- length(gene_lists)
  
  if (n <= 1) return(pathway_results)
  
  # Compute pairwise Jaccard index
  jaccard_mat <- matrix(0, nrow = n, ncol = n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      g1 <- gene_lists[[i]]
      g2 <- gene_lists[[j]]
      intersect_len <- length(intersect(g1, g2))
      union_len <- length(union(g1, g2))
      j_idx <- intersect_len / union_len
      jaccard_mat[i, j] <- j_idx
      jaccard_mat[j, i] <- j_idx
    }
  }
  
  # Convert to distance matrix (1 - Jaccard)
  dist_mat <- as.dist(1 - jaccard_mat)
  
  # Hierarchical clustering
  hc <- hclust(dist_mat, method = "average")
  
  # Cut tree at height corresponding to 1 - jaccard_threshold
  clusters <- cutree(hc, h = 1 - jaccard_threshold)
  
  # Select the most significant pathway per cluster
  pathway_results$cluster <- clusters
  
  # Ensure padj is not NA, replace NA with 1 for sorting purposes
  sort_padj <- pathway_results$padj
  sort_padj[is.na(sort_padj)] <- 1
  
  keep_idx <- integer(0)
  for (cl in unique(clusters)) {
    cl_idx <- which(clusters == cl)
    # Get index of minimum padj
    best <- cl_idx[which.min(sort_padj[cl_idx])]
    keep_idx <- c(keep_idx, best)
  }
  
  dedup <- pathway_results[sort(keep_idx), ]
  dedup$cluster <- NULL
  dedup
}

#' Cross-study pathway heatmap
#'
#' Creates a heatmap showing enrichment scores or -log10(padj) for top
#' pathways across all studies.
#'
#' @param pathway_results A list of per-study pathway enrichment
#'   \code{data.frame} objects.
#' @param top_n Integer. Number of top pathways to display. Default: \code{30}.
#' @param value Character scalar. Value to display: \code{"padj"} (default)
#'   or \code{"NES"} (for GSEA results).
#'
#' @return A \code{ComplexHeatmap::Heatmap} object.
#'
#' @examples
#' \dontrun{
#'   mx_pathway_heatmap(per_study_pathways, top_n = 20)
#' }
#'
#' @seealso \code{\link{mx_pathway_meta}}
#'
#' @export
mx_pathway_heatmap <- function(pathway_results, top_n = 30,
                                value = c("padj", "NES")) {
  if (!is.list(pathway_results) || is.data.frame(pathway_results)) {
    stop("'pathway_results' must be a list of data.frames from per-study enrichment")
  }
  
  value <- match.arg(value)
  
  # Extract top pathways based on significance across all studies
  all_paths <- do.call(rbind, lapply(pathway_results, function(x) {
    if(nrow(x) == 0) return(NULL)
    x[, c("pathway_id", "pathway_name", "padj", if(value == "NES" && "NES" %in% colnames(x)) "NES" else NULL)]
  }))
  
  if(is.null(all_paths) || nrow(all_paths) == 0) {
    stop("No pathways to plot.")
  }
  
  # Get top N by lowest median padj
  path_medians <- aggregate(padj ~ pathway_id, data = all_paths, FUN = median, na.rm = TRUE)
  top_paths <- head(path_medians$pathway_id[order(path_medians$padj)], top_n)
  
  k <- length(pathway_results)
  mat <- matrix(NA_real_, nrow = length(top_paths), ncol = k, 
                dimnames = list(top_paths, names(pathway_results)))
  
  # Fill matrix
  for(i in seq_len(k)) {
    pr <- pathway_results[[i]]
    if(nrow(pr) == 0) next
    
    idx <- match(top_paths, pr$pathway_id)
    if(value == "padj") {
      mat[, i] <- -log10(pr$padj[idx])
    } else {
      mat[, i] <- pr$NES[idx]
    }
  }
  
  # Get pathway names
  path_names <- all_paths$pathway_name[match(top_paths, all_paths$pathway_id)]
  rownames(mat) <- path_names
  
  # Clean up NAs
  mat[is.na(mat)] <- 0
  
  # Color scale
  if(value == "padj") {
    col_fun <- circlize::colorRamp2(c(0, max(mat, na.rm=TRUE)), c("white", "red"))
    title <- "-log10(padj)"
  } else {
    max_val <- max(abs(mat), na.rm=TRUE)
    col_fun <- circlize::colorRamp2(c(-max_val, 0, max_val), c("blue", "white", "red"))
    title <- "NES"
  }
  
  ComplexHeatmap::Heatmap(
    mat, 
    name = title,
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 8),
    column_title = paste0("Top ", length(top_paths), " Consensus Pathways")
  )
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
    geneID       = res_df$geneID,
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
    geneID       = res_df$core_enrichment,
    stringsAsFactors = FALSE
  )
}
