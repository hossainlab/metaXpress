#' @include AllClasses.R
NULL

# ============================================================================
# Module 6: Pathway Meta-Analysis
# ============================================================================

#' Retrieve gene sets from MSigDB
#'
#' Downloads a gene set collection with \pkg{msigdbr} and returns it as a
#' named list, using the identifier type that matches your gene IDs.
#'
#' @param db Character scalar. One of \code{"Hallmarks"} (default),
#'   \code{"Reactome"}, \code{"KEGG"} (KEGG MEDICUS), \code{"WikiPathways"},
#'   \code{"GO_BP"}.
#' @param species Species name understood by \pkg{msigdbr}. Default:
#'   \code{"Homo sapiens"}.
#' @param id_type Gene identifier: \code{"symbol"} (default), \code{"entrez"}
#'   or \code{"ensembl"}. NCBI GEO counts from \code{\link{mx_fetch_geo}} use
#'   Entrez IDs.
#'
#' @return A named list of character vectors (gene set name -> gene IDs).
#'   The MSigDB version is stored in \code{attr(x, "version")}.
#'
#' @references
#' Liberzon, A. et al. (2015) The Molecular Signatures Database (MSigDB)
#' hallmark gene set collection. \emph{Cell Systems}, \strong{1}(6),
#' 417--425. \doi{10.1016/j.cels.2015.12.004}
#'
#' @examples
#' if (requireNamespace("msigdbr", quietly = TRUE)) {
#'   hallmarks <- mx_gene_sets("Hallmarks")
#'   length(hallmarks)
#'   head(hallmarks$HALLMARK_INFLAMMATORY_RESPONSE)
#' }
#'
#' @seealso \code{\link{mx_pathway_meta}}, \code{\link{mx_pathway_studies}}
#' @export
mx_gene_sets <- function(db = c("Hallmarks", "Reactome", "KEGG",
                                "WikiPathways", "GO_BP"),
                         species = "Homo sapiens",
                         id_type = c("symbol", "entrez", "ensembl")) {
  db      <- match.arg(db)
  id_type <- match.arg(id_type)
  .require_package("msigdbr", "MSigDB gene sets")

  spec <- switch(db,
    Hallmarks    = list("H",  NULL),
    Reactome     = list("C2", "CP:REACTOME"),
    KEGG         = list("C2", "CP:KEGG_MEDICUS"),
    WikiPathways = list("C2", "CP:WIKIPATHWAYS"),
    GO_BP        = list("C5", "GO:BP")
  )
  gs <- tryCatch(
    msigdbr::msigdbr(species = species, collection = spec[[1]],
                     subcollection = spec[[2]]),
    error = function(e)
      msigdbr::msigdbr(species = species, category = spec[[1]],
                       subcategory = spec[[2]]))
  col <- switch(id_type, symbol = "gene_symbol", entrez = "ncbi_gene",
                ensembl = "ensembl_gene")
  ids <- as.character(gs[[col]])
  ok  <- !is.na(ids) & ids != ""
  out <- lapply(split(ids[ok], gs$gs_name[ok]), unique)
  attr(out, "version") <- if ("db_version" %in% names(gs))
    unique(gs$db_version)[1] else NA_character_
  out
}

#' Pathway enrichment of meta-analysis results
#'
#' Runs pre-ranked gene set enrichment analysis (GSEA, via \pkg{fgsea}) or
#' over-representation analysis (ORA, hypergeometric test) on the
#' gene-level meta-analysis.
#'
#' @param meta_result A \code{\linkS4class{metaXpressResult}} from
#'   \code{\link{mx_meta}}.
#' @param db Gene set collection passed to \code{\link{mx_gene_sets}}.
#'   Ignored when \code{gene_sets} is supplied.
#' @param method \code{"GSEA"} (default) or \code{"ORA"}.
#' @param gene_sets Optional named list of gene sets (e.g. from
#'   \code{\link{mx_gene_sets}} or a GMT file) whose IDs match
#'   \code{meta_table$gene_id}.
#' @param padj_threshold,lfc_threshold ORA only: thresholds on
#'   \code{meta_padj} and \code{|meta_log2FC|} defining the gene lists.
#'   Defaults: \code{0.05} and \code{0}.
#' @param min_size,max_size Gene set size limits after restricting to genes
#'   in the analysis. Defaults: \code{15} and \code{500}.
#' @param seed Integer seed for the fgsea permutation sampling, for
#'   reproducibility. The global RNG state is restored afterwards.
#' @param ... Further arguments passed to \code{fgsea::fgsea()}.
#'
#' @return The \code{meta_result} with \code{pathway_result} filled:
#'   \describe{
#'     \item{GSEA}{\code{pathway_id}, \code{set_size}, \code{NES},
#'       \code{pvalue}, \code{padj}, \code{leading_edge}.}
#'     \item{ORA}{\code{pathway_id}, \code{direction} (up/down),
#'       \code{set_size}, \code{n_overlap}, \code{n_list}, \code{fold_enrichment},
#'       \code{pvalue}, \code{padj}, \code{overlap_genes}.}
#'   }
#'
#' @details
#' GSEA ranks genes by the signed meta-analysis z-score
#' \eqn{\mathrm{sign}(\beta)\,\Phi^{-1}(1 - p/2)}, which uses both effect
#' size and its precision; this avoids arbitrary significance cut-offs and is
#' the recommended default. ORA tests up- and down-regulated genes
#' separately against the universe of all genes in the meta-analysis
#' (not the whole genome), as recommended to avoid inflated enrichment.
#'
#' @references
#' Subramanian, A. et al. (2005) Gene set enrichment analysis: a
#' knowledge-based approach for interpreting genome-wide expression profiles.
#' \emph{PNAS}, \strong{102}(43), 15545--15550.
#' \doi{10.1073/pnas.0506580102}
#'
#' Korotkevich, G. et al. (2016) Fast gene set enrichment analysis.
#' \emph{bioRxiv}. \doi{10.1101/060012}
#'
#' Timmons, J.A., Szkop, K.J. & Gallagher, I.J. (2015) Multiple sources of
#' bias confound functional enrichment analysis of global -omics data.
#' \emph{Genome Biology}, \strong{16}, 186. \doi{10.1186/s13059-015-0761-7}
#'
#' @examples
#' if (requireNamespace("fgsea", quietly = TRUE)) {
#'   set.seed(1)
#'   genes <- paste0("G", 1:1000)
#'   sets  <- list(UP_SET = genes[1:40], NULL_SET = genes[501:560])
#'   de <- lapply(1:3, function(s) {
#'     lfc <- rnorm(1000, c(rep(1, 40), rep(0, 960)), 0.3)
#'     p <- 2 * pnorm(-abs(lfc / 0.3))
#'     data.frame(gene_id = genes, log2FC = lfc, lfcSE = 0.3, pvalue = p,
#'                padj = p.adjust(p, "BH"))
#'   })
#'   res <- mx_meta(de)
#'   res <- mx_pathway_meta(res, gene_sets = sets, min_size = 10)
#'   res@pathway_result
#' }
#'
#' @seealso \code{\link{mx_pathway_studies}}, \code{\link{mx_gene_sets}},
#'   \code{\link{mx_pathway_dedup}}, \code{\link{mx_pathway_heatmap}}
#'
#' @importFrom methods is
#' @export
mx_pathway_meta <- function(meta_result,
                            db = c("Hallmarks", "Reactome", "KEGG",
                                   "WikiPathways", "GO_BP"),
                            method = c("GSEA", "ORA"),
                            gene_sets = NULL,
                            padj_threshold = 0.05,
                            lfc_threshold  = 0,
                            min_size = 15,
                            max_size = 500,
                            seed = 1L,
                            ...) {
  if (!is(meta_result, "metaXpressResult"))
    stop("'meta_result' must be a metaXpressResult object")
  method <- match.arg(method)
  mt <- meta_result@meta_table
  gene_sets <- .resolve_gene_sets(gene_sets, match.arg(db), mt$gene_id)

  pathway_result <- switch(method,
    GSEA = .pathway_gsea(.signed_z(mt$meta_pvalue, mt$meta_log2FC, mt$gene_id),
                         gene_sets, min_size, max_size, seed, ...),
    ORA  = .pathway_ora(mt, gene_sets, padj_threshold, lfc_threshold,
                        min_size, max_size)
  )
  attr(pathway_result, "method") <- method
  meta_result@pathway_result <- pathway_result
  meta_result
}

#' Cross-study pathway meta-analysis
#'
#' Runs pre-ranked GSEA separately in every study and combines the
#' pathway-level evidence across studies. This answers whether a pathway is
#' consistently perturbed across cohorts, which is more robust than
#' enrichment of a single pooled gene list when studies differ in platform,
#' depth or gene coverage (Shen & Tseng 2010).
#'
#' @param de_results A named list of per-study DE \code{data.frame}s or
#'   \code{\linkS4class{metaXpressStudy}} objects after
#'   \code{\link{mx_de_all}}.
#' @inheritParams mx_pathway_meta
#' @param weights Per-study weights for the Stouffer combination:
#'   \code{"equal"} (default), \code{"sample_size"} (square root of the
#'   per-study sample size from \code{n_case} / \code{n_control}), or a
#'   numeric vector with one weight per study. Equal weights are the default
#'   because GSEA p-values of a large cohort are already far more extreme,
#'   and sample-size weights would let one cohort dominate.
#'
#' @return A \code{data.frame} with one row per gene set: \code{pathway_id},
#'   \code{n_studies}, \code{mean_NES}, \code{n_up} / \code{n_down}
#'   (studies with positive / negative NES), \code{n_sig} (studies with
#'   per-study \code{padj <= 0.05}), \code{meta_z}, \code{meta_pvalue},
#'   \code{meta_padj}, followed by one \code{NES_<study>} column per study.
#'   Per-study fgsea tables are in \code{attr(x, "per_study")}.
#'
#' @details
#' For each study, genes are ranked by the signed z-score derived from the
#' per-study test (the \code{stat} column from \code{\link{mx_de}} when
#' present). Each pathway's two-sided GSEA p-value is converted to a signed
#' z-score using the sign of its NES, and the z-scores are combined by
#' weighted Stouffer's method, so pathways enriched in opposite directions
#' cancel out. This is the p-value-combination variant of MAPE (Shen & Tseng
#' 2010) with directional combination.
#'
#' @references
#' Shen, K. & Tseng, G.C. (2010) Meta-analysis for pathway enrichment
#' analysis when combining multiple genomic studies. \emph{Bioinformatics},
#' \strong{26}(10), 1316--1323. \doi{10.1093/bioinformatics/btq148}
#'
#' Korotkevich, G. et al. (2016) Fast gene set enrichment analysis.
#' \emph{bioRxiv}. \doi{10.1101/060012}
#'
#' @examples
#' if (requireNamespace("fgsea", quietly = TRUE)) {
#'   set.seed(2)
#'   genes <- paste0("G", 1:1000)
#'   sets  <- list(UP_SET = genes[1:40], NULL_SET = genes[501:560])
#'   de <- lapply(1:3, function(s) {
#'     lfc <- rnorm(1000, c(rep(0.8, 40), rep(0, 960)), 0.3)
#'     p <- 2 * pnorm(-abs(lfc / 0.3))
#'     data.frame(gene_id = genes, log2FC = lfc, pvalue = p,
#'                padj = p.adjust(p, "BH"))
#'   })
#'   names(de) <- c("S1", "S2", "S3")
#'   mx_pathway_studies(de, gene_sets = sets, min_size = 10)
#' }
#'
#' @seealso \code{\link{mx_pathway_meta}}, \code{\link{mx_pathway_heatmap}}
#' @export
mx_pathway_studies <- function(de_results,
                               db = c("Hallmarks", "Reactome", "KEGG",
                                      "WikiPathways", "GO_BP"),
                               gene_sets = NULL,
                               min_size = 15,
                               max_size = 500,
                               weights = "equal",
                               seed = 1L,
                               ...) {
  de_results <- .extract_de_results(de_results)
  .validate_de_list(de_results)
  if (is.null(names(de_results)))
    names(de_results) <- paste0("study_", seq_along(de_results))
  k <- length(de_results)
  all_ids <- unique(unlist(lapply(de_results, `[[`, "gene_id")))
  gene_sets <- .resolve_gene_sets(gene_sets, match.arg(db), all_ids)

  w <- if (is.numeric(weights)) {
    if (length(weights) != k) stop("'weights' must have length ", k)
    as.numeric(weights)
  } else if (identical(weights, "sample_size")) {
    n <- vapply(de_results, function(d)
      if (all(c("n_case", "n_control") %in% names(d)))
        as.numeric(d$n_case[1] + d$n_control[1]) else NA_real_, numeric(1))
    if (anyNA(n)) stop("'sample_size' weights need n_case/n_control columns")
    sqrt(n)
  } else if (identical(weights, "equal")) {
    rep(1, k)
  } else {
    stop("'weights' must be 'equal', 'sample_size' or a numeric vector")
  }

  message("Running per-study GSEA on ", k, " studies...")
  per_study <- lapply(names(de_results), function(nm) {
    d <- de_results[[nm]]
    stat <- if ("stat" %in% names(d) && !all(is.na(d$stat))) d$stat
            else .signed_z(d$pvalue, d$log2FC)
    .pathway_gsea(stats::setNames(stat, d$gene_id), gene_sets,
                  min_size, max_size, seed, ...)
  })
  names(per_study) <- names(de_results)

  pids <- unique(unlist(lapply(per_study, `[[`, "pathway_id")))
  get <- function(col) vapply(per_study, function(r)
    r[[col]][match(pids, r$pathway_id)], numeric(length(pids)))
  nes  <- matrix(get("NES"), nrow = length(pids))
  pv   <- matrix(get("pvalue"), nrow = length(pids))
  padj <- matrix(get("padj"), nrow = length(pids))

  z <- sign(nes) * stats::qnorm(pmax(pv, .Machine$double.xmin) / 2,
                                lower.tail = FALSE)
  wm <- matrix(w, nrow(z), ncol(z), byrow = TRUE)
  wm[is.na(z)] <- NA_real_
  meta_z <- rowSums(wm * z, na.rm = TRUE) / sqrt(rowSums(wm^2, na.rm = TRUE))
  n_st   <- rowSums(!is.na(z))
  meta_p <- 2 * stats::pnorm(-abs(meta_z))
  meta_p[n_st < 2] <- NA_real_

  out <- data.frame(
    pathway_id  = pids,
    n_studies   = as.integer(n_st),
    mean_NES    = rowMeans(nes, na.rm = TRUE),
    n_up        = as.integer(rowSums(nes > 0, na.rm = TRUE)),
    n_down      = as.integer(rowSums(nes < 0, na.rm = TRUE)),
    n_sig       = as.integer(rowSums(padj <= 0.05, na.rm = TRUE)),
    meta_z      = meta_z,
    meta_pvalue = meta_p,
    meta_padj   = stats::p.adjust(meta_p, method = "BH"),
    stringsAsFactors = FALSE
  )
  nes_df <- as.data.frame(nes)
  names(nes_df) <- paste0("NES_", names(de_results))
  out <- cbind(out, nes_df)
  out <- out[order(out$meta_pvalue), ]
  rownames(out) <- NULL
  attr(out, "per_study") <- per_study
  out
}

#' Find consensus pathways across studies (vote counting)
#'
#' Counts in how many per-study enrichment analyses each pathway is
#' significant. Vote counting ignores effect size and is less powerful than
#' \code{\link{mx_pathway_studies}}; use it as a descriptive summary.
#'
#' @param pathway_results A list of per-study enrichment \code{data.frame}s
#'   with columns \code{pathway_id} and \code{padj}, e.g.
#'   \code{attr(mx_pathway_studies(...), "per_study")}.
#' @param min_fraction Numeric. Minimum fraction of studies in which a pathway
#'   must be significant. Default: \code{0.5}.
#' @param padj_threshold Numeric. Significance threshold per study.
#'   Default: \code{0.05}.
#'
#' @return A \code{data.frame} with \code{pathway_id},
#'   \code{n_significant_studies} and \code{fraction_studies}.
#'
#' @examples
#' pr1 <- data.frame(pathway_id = c("A", "B"), padj = c(0.01, 0.2))
#' pr2 <- data.frame(pathway_id = c("A", "B"), padj = c(0.03, 0.01))
#' mx_pathway_consensus(list(pr1, pr2), min_fraction = 1)
#'
#' @seealso \code{\link{mx_pathway_studies}}
#' @export
mx_pathway_consensus <- function(pathway_results, min_fraction = 0.5,
                                 padj_threshold = 0.05) {
  if (!is.list(pathway_results) || is.data.frame(pathway_results))
    stop("'pathway_results' must be a list of pathway enrichment data.frames")

  k         <- length(pathway_results)
  all_paths <- unique(unlist(lapply(pathway_results, `[[`, "pathway_id")))

  n_sig <- vapply(all_paths, function(pid) {
    sum(vapply(pathway_results, function(r) {
      p <- r$padj[match(pid, r$pathway_id)]
      !is.na(p) && p <= padj_threshold
    }, logical(1)))
  }, integer(1))

  fraction <- n_sig / k
  keep     <- fraction >= min_fraction

  data.frame(
    pathway_id            = all_paths[keep],
    n_significant_studies = n_sig[keep],
    fraction_studies      = fraction[keep],
    stringsAsFactors      = FALSE,
    row.names             = NULL
  )
}

#' Remove redundant pathways
#'
#' Greedy redundancy filter: pathways are visited from most to least
#' significant and a pathway is dropped when its gene-set Jaccard similarity
#' with an already retained pathway is at least \code{jaccard_threshold}.
#'
#' @param pathway_results A \code{data.frame} with \code{pathway_id} and a
#'   p-value column (\code{padj} or \code{meta_padj}).
#' @param gene_sets The named list of gene sets used for the enrichment.
#'   When \code{NULL}, results are returned unchanged with a warning.
#' @param jaccard_threshold Numeric. Default: \code{0.5}.
#'
#' @return The filtered \code{data.frame}, ordered by significance.
#'
#' @examples
#' sets <- list(A = letters[1:10], B = letters[1:9], C = letters[15:25])
#' res  <- data.frame(pathway_id = c("A", "B", "C"), padj = c(0.01, 0.02, 0.03))
#' mx_pathway_dedup(res, gene_sets = sets)
#'
#' @export
mx_pathway_dedup <- function(pathway_results, gene_sets = NULL,
                             jaccard_threshold = 0.5) {
  if (!is.data.frame(pathway_results))
    stop("'pathway_results' must be a data.frame from mx_pathway_meta()")
  if (nrow(pathway_results) == 0)
    return(pathway_results)
  if (!"pathway_id" %in% colnames(pathway_results))
    stop("'pathway_results' must contain a 'pathway_id' column")
  if (is.null(gene_sets)) {
    warning("'gene_sets' not supplied; returning results unchanged.")
    return(pathway_results)
  }

  pcol <- intersect(c("meta_padj", "padj", "pvalue"),
                    colnames(pathway_results))[1]
  if (!is.na(pcol))
    pathway_results <- pathway_results[order(pathway_results[[pcol]],
                                             na.last = TRUE), , drop = FALSE]
  pids <- pathway_results$pathway_id
  keep <- rep(TRUE, length(pids))
  kept <- list()
  for (i in seq_along(pids)) {
    gi <- gene_sets[[pids[i]]]
    if (is.null(gi)) next
    for (gj in kept) {
      if (length(intersect(gi, gj)) / length(union(gi, gj)) >=
          jaccard_threshold) {
        keep[i] <- FALSE
        break
      }
    }
    if (keep[i]) kept[[length(kept) + 1]] <- gi
  }
  pathway_results[keep, , drop = FALSE]
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
  value <- match.arg(value)

  # Accept single data.frame or named list
  if (is.data.frame(pathway_results))
    pathway_results <- list(Combined = pathway_results)

  if (!is.list(pathway_results))
    stop("'pathway_results' must be a data.frame or named list of data.frames")

  study_names <- names(pathway_results)
  if (is.null(study_names))
    study_names <- paste0("Study_", seq_along(pathway_results))

  all_pathways <- unique(unlist(lapply(pathway_results,
                                        function(r) r$pathway_id)))
  if (length(all_pathways) == 0)
    stop("No pathways found in 'pathway_results'")

  # Select top_n pathways by mean significance
  if (value == "padj") {
    mean_padj <- vapply(all_pathways, function(pid) {
      vals <- vapply(pathway_results, function(r) {
        row <- r[r$pathway_id == pid, ]
        if (nrow(row) == 0 || is.na(row$padj[1])) 1 else row$padj[1]
      }, numeric(1))
      mean(vals)
    }, numeric(1))
    top_paths <- all_pathways[order(mean_padj)][seq_len(min(top_n,
                                                             length(all_pathways)))]
  } else {
    mean_abs_nes <- vapply(all_pathways, function(pid) {
      vals <- vapply(pathway_results, function(r) {
        row <- r[r$pathway_id == pid, ]
        if (nrow(row) == 0 || !"NES" %in% names(row)) NA_real_ else row$NES[1]
      }, numeric(1))
      mean(abs(vals), na.rm = TRUE)
    }, numeric(1))
    top_paths <- all_pathways[order(-mean_abs_nes)][
      seq_len(min(top_n, length(all_pathways)))]
  }

  # Build matrix
  val_mat <- matrix(NA_real_, nrow = length(top_paths),
                    ncol = length(pathway_results),
                    dimnames = list(top_paths, study_names))

  for (i in seq_along(pathway_results)) {
    r <- pathway_results[[i]]
    for (pid in top_paths) {
      row <- r[r$pathway_id == pid, ]
      if (nrow(row) == 0) next
      val_mat[pid, i] <- if (value == "NES" && "NES" %in% colnames(row))
        row$NES[1]
      else
        -log10(pmax(row$padj[1], .Machine$double.eps))
    }
  }

  # Replace NA with 0 for display
  val_mat[is.na(val_mat)] <- 0

  # Truncate long pathway names for display
  row_labels <- sub("^HALLMARK_|^KEGG_|^REACTOME_|^GOBP_", "",
                    rownames(val_mat))
  row_labels <- substr(row_labels, 1, 40)

  # Build tidy data.frame for ggplot2
  plot_df <- data.frame(
    pathway = rep(row_labels, ncol(val_mat)),
    study   = rep(study_names, each = nrow(val_mat)),
    value   = as.vector(val_mat),
    stringsAsFactors = FALSE
  )
  plot_df$pathway <- factor(plot_df$pathway,
                             levels = rev(row_labels))
  plot_df$study   <- factor(plot_df$study, levels = study_names)

  y_label <- if (value == "NES") "NES" else "-log10(padj)"
  mid_val  <- if (value == "NES") 0 else mean(range(plot_df$value, na.rm = TRUE))

  ggplot2::ggplot(plot_df, ggplot2::aes(x = study, y = pathway,
                                          fill = value)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::scale_fill_gradient2(
      low      = "#377EB8",
      mid      = "white",
      high     = "#E41A1C",
      midpoint = mid_val,
      name     = y_label
    ) +
    ggplot2::labs(title = paste0("Top ", length(top_paths),
                                  " pathways - ", y_label),
                  x = NULL, y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 7)
    )
}

# ============================================================================
# Internal helpers
# ============================================================================

# Signed z-score from two-sided p-value and effect direction
.signed_z <- function(p, lfc, ids = NULL) {
  z <- sign(lfc) * stats::qnorm(pmax(p, .Machine$double.xmin) / 2,
                                lower.tail = FALSE)
  if (!is.null(ids)) names(z) <- ids
  z
}

# Use supplied gene sets or download them, matching the gene ID type
.resolve_gene_sets <- function(gene_sets, db, gene_ids) {
  if (is.null(gene_sets)) {
    id_type <- switch(.detect_id_type(gene_ids),
                      ENTREZID = "entrez", ENSEMBL = "ensembl", "symbol")
    message("Using MSigDB ", db, " gene sets (", id_type, " IDs)")
    gene_sets <- mx_gene_sets(db, id_type = id_type)
  }
  if (!is.list(gene_sets) || is.null(names(gene_sets)))
    stop("'gene_sets' must be a named list of character vectors")
  overlap <- length(intersect(unique(unlist(gene_sets)), gene_ids))
  if (overlap == 0)
    stop("No overlap between gene set members and gene IDs. Check that ",
         "both use the same identifier type (symbol / Entrez / Ensembl).")
  gene_sets
}

.pathway_gsea <- function(stats_vec, gene_sets, min_size, max_size, seed,
                          ...) {
  .require_package("fgsea", "GSEA")
  stats_vec <- stats_vec[is.finite(stats_vec) & !is.na(names(stats_vec))]
  stats_vec <- stats_vec[!duplicated(names(stats_vec))]

  if (!is.null(seed)) {
    has_seed <- exists(".Random.seed", envir = globalenv())
    if (has_seed) old <- get(".Random.seed", envir = globalenv())
    on.exit(if (has_seed) assign(".Random.seed", old, envir = globalenv())
            else rm(".Random.seed", envir = globalenv()))
    set.seed(seed)
  }
  res <- suppressWarnings(fgsea::fgsea(pathways = gene_sets,
                                       stats = stats_vec,
                                       minSize = min_size,
                                       maxSize = max_size, ...))
  if (nrow(res) == 0)
    return(data.frame(pathway_id = character(0), set_size = integer(0),
                      NES = numeric(0), pvalue = numeric(0),
                      padj = numeric(0), leading_edge = character(0)))
  out <- data.frame(
    pathway_id   = res$pathway,
    set_size     = as.integer(res$size),
    NES          = res$NES,
    pvalue       = res$pval,
    padj         = res$padj,
    leading_edge = vapply(res$leadingEdge, paste, character(1),
                          collapse = ";"),
    stringsAsFactors = FALSE
  )
  out[order(out$pvalue), , drop = FALSE]
}

.pathway_ora <- function(meta_table, gene_sets, padj_threshold, lfc_threshold,
                         min_size, max_size) {
  universe <- unique(meta_table$gene_id[!is.na(meta_table$meta_padj)])
  sig <- !is.na(meta_table$meta_padj) &
    meta_table$meta_padj <= padj_threshold &
    abs(meta_table$meta_log2FC) >= lfc_threshold
  lists <- list(
    up   = unique(meta_table$gene_id[sig & meta_table$meta_log2FC > 0]),
    down = unique(meta_table$gene_id[sig & meta_table$meta_log2FC < 0]))

  sets <- lapply(gene_sets, intersect, universe)
  sets <- sets[lengths(sets) >= min_size & lengths(sets) <= max_size]
  if (length(sets) == 0) stop("No gene sets within the size limits")
  N <- length(universe)

  rows <- lapply(names(lists), function(dir) {
    g <- lists[[dir]]
    if (length(g) == 0) return(NULL)
    ov <- vapply(sets, function(s) length(intersect(s, g)), integer(1))
    M  <- lengths(sets)
    p  <- stats::phyper(ov - 1, M, N - M, length(g), lower.tail = FALSE)
    data.frame(
      pathway_id      = names(sets),
      direction       = dir,
      set_size        = as.integer(M),
      n_overlap       = ov,
      n_list          = length(g),
      fold_enrichment = (ov / length(g)) / (M / N),
      pvalue          = p,
      padj            = stats::p.adjust(p, method = "BH"),
      overlap_genes   = vapply(sets, function(s)
        paste(intersect(s, g), collapse = ";"), character(1)),
      stringsAsFactors = FALSE, row.names = NULL)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    warning("No significant genes at meta_padj <= ", padj_threshold,
            "; returning empty ORA result.")
    return(data.frame())
  }
  out[order(out$pvalue), , drop = FALSE]
}
