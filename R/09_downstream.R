#' @include AllClasses.R
NULL

# ============================================================================
# Module 9: Downstream Biology
# ============================================================================
#
# Functions that turn a meta-signature into biological hypotheses:
#   - mx_tf_activity()           upstream transcription factor activity
#   - mx_cell_composition()      cross-study cell-type composition shifts
#   - mx_coexpression_modules()  consensus co-expression modules
#   - mx_drug_connectivity()     drug signatures that reverse the meta-signature
#
# External resources (regulons, marker genes, drug signatures) are passed in
# as plain R objects so every analysis runs offline and with any source.

#' Infer transcription factor activity from the meta-signature
#'
#' Estimates the activity of each transcription factor (TF) from the
#' expression changes of its target genes, using a univariate linear model
#' (ULM): for each TF, the gene-level meta statistic is regressed on the TF's
#' signed target weights, and the t-value of the slope is the activity score.
#' This is the default method of decoupleR and is well suited to CollecTRI or
#' DoRothEA regulons.
#'
#' @param meta_result A \code{\linkS4class{metaXpressResult}} object.
#' @param network A \code{data.frame} of TF-target interactions with columns
#'   \code{source} (TF), \code{target} (gene) and \code{mor} (mode of
#'   regulation: positive for activation, negative for repression). OmniPath
#'   exports (\code{source_genesymbol}, \code{target_genesymbol},
#'   \code{is_stimulation}, \code{is_inhibition}) are converted automatically.
#' @param id_col Column of \code{meta_table} whose values match
#'   \code{network$target}. Default \code{"gene_id"}; use \code{"symbol"} when
#'   gene IDs are Entrez/Ensembl and a symbol column is present.
#' @param stat Gene-level statistic: \code{"z"} (signed meta z-score, default)
#'   or \code{"log2FC"}.
#' @param min_size Minimum number of measured targets per TF. Default 5.
#'
#' @return A \code{data.frame} sorted by p-value with columns \code{source},
#'   \code{score} (positive = TF more active in cases), \code{pvalue},
#'   \code{padj} (BH) and \code{n_targets}.
#'
#' @references
#' Badia-i-Mompel, P. et al. (2022) decoupleR: ensemble of computational
#' methods to infer biological activities from omics data.
#' \emph{Bioinformatics Advances}, \strong{2}(1), vbac016.
#' \doi{10.1093/bioadv/vbac016}
#'
#' Müller-Dott, S. et al. (2023) Expanding the coverage of regulons from
#' high-confidence prior knowledge for accurate estimation of transcription
#' factor activities. \emph{Nucleic Acids Research}, \strong{51}(20),
#' 10934--10949. \doi{10.1093/nar/gkad841}
#'
#' @examples
#' mt <- data.frame(gene_id = paste0("G", 1:200),
#'                  meta_log2FC = c(rep(2, 20), rnorm(180, sd = 0.3)),
#'                  meta_pvalue = c(rep(1e-6, 20), runif(180)),
#'                  meta_padj = 1, i_squared = 0, n_studies = 3L)
#' res <- new("metaXpressResult", meta_table = mt, method = "fisher",
#'            n_studies = 3L)
#' net <- data.frame(source = rep(c("TF_A", "TF_B"), each = 20),
#'                   target = paste0("G", c(1:20, 101:120)), mor = 1)
#' mx_tf_activity(res, net)
#'
#' @seealso \code{\link{mx_drug_connectivity}}, \code{\link{mx_downstream_plot}}
#' @export
mx_tf_activity <- function(meta_result, network, id_col = "gene_id",
                           stat = c("z", "log2FC"), min_size = 5) {
  stat <- match.arg(stat)
  y    <- .gene_statistic(meta_result, id_col, stat)
  net  <- .standardize_network(network)
  net  <- net[net$target %in% names(y), , drop = FALSE]

  tfs <- split(net, net$source)
  tfs <- tfs[vapply(tfs, nrow, integer(1)) >= min_size]
  if (length(tfs) == 0)
    stop("No TF has >= ", min_size, " targets among the measured genes. ",
         "Check that 'id_col' matches the identifiers in 'network$target'.")

  n  <- length(y)
  yc <- y - mean(y)
  syy <- sum(yc^2)

  rows <- lapply(names(tfs), function(tf) {
    x <- numeric(n)
    names(x) <- names(y)
    x[tfs[[tf]]$target] <- tfs[[tf]]$mor
    xc <- x - mean(x)
    r  <- sum(xc * yc) / sqrt(sum(xc^2) * syy)
    t  <- r * sqrt((n - 2) / max(1 - r^2, .Machine$double.eps))
    data.frame(source = tf, score = t,
               pvalue = 2 * stats::pt(-abs(t), df = n - 2),
               n_targets = nrow(tfs[[tf]]), stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, rows)
  out$padj <- stats::p.adjust(out$pvalue, method = "BH")
  out <- out[order(out$pvalue), c("source", "score", "pvalue", "padj",
                                  "n_targets")]
  rownames(out) <- NULL
  .tag(out, "tf_activity")
}

#' Meta-analyse cell-type composition shifts across studies
#'
#' Scores each sample for each cell type as the mean within-study z-score of
#' log-CPM expression across that cell type's marker genes (an MCP-counter
#' style abundance proxy), computes Hedges' g for case vs control in every
#' study, and combines the per-study effects with a DerSimonian-Laird
#' random-effects model. A significant shift indicates that part of the
#' bulk differential expression may reflect changes in tissue composition.
#'
#' @param studies A named list of \code{\linkS4class{metaXpressStudy}} objects.
#' @param markers A named list mapping cell types to marker gene identifiers,
#'   or a \code{data.frame} with columns \code{gene} and \code{cell_type}.
#'   Identifiers must match \code{rownames} of the studies.
#' @param min_markers Minimum number of expressed markers a cell type needs in
#'   a study to be scored there. Default 3.
#'
#' @details Expression is log2-CPM. If a study's \code{colData} has a
#'   \code{lib_size} column it is used as the library size, so studies that
#'   carry only a subset of genes are still normalised to their true depth.
#'
#' @return A \code{data.frame} sorted by p-value with columns
#'   \code{cell_type}, \code{meta_g} (pooled Hedges' g; positive = higher in
#'   cases), \code{pvalue}, \code{padj}, \code{i_squared}, \code{n_studies},
#'   and one \code{g_<study>} column per study. Per-sample scores are stored
#'   in \code{attr(x, "scores")}.
#'
#' @references
#' Becht, E. et al. (2016) Estimating the population abundance of
#' tissue-infiltrating immune and stromal cell populations using gene
#' expression. \emph{Genome Biology}, \strong{17}, 218.
#' \doi{10.1186/s13059-016-1070-5}
#'
#' Hedges, L. V. (1981) Distribution theory for Glass's estimator of effect
#' size and related estimators. \emph{Journal of Educational Statistics},
#' \strong{6}(2), 107--128. \doi{10.3102/10769986006002107}
#'
#' @examples
#' \dontrun{
#'   markers <- list(T_cells = c("CD3D", "CD3E", "CD2"),
#'                   Fibroblasts = c("COL1A1", "COL3A1", "DCN"))
#'   comp <- mx_cell_composition(studies, markers)
#' }
#'
#' @seealso \code{\link{mx_coexpression_modules}},
#'   \code{\link{mx_downstream_plot}}
#' @export
mx_cell_composition <- function(studies, markers, min_markers = 3) {
  studies <- .named_study_list(studies)
  markers <- .standardize_markers(markers)

  scores <- lapply(studies, function(s) {
    z <- .row_z(.log_cpm(s))
    sc <- vapply(markers, function(g) {
      g <- intersect(g, rownames(z))
      if (length(g) < min_markers) return(rep(NA_real_, ncol(z)))
      colMeans(z[g, , drop = FALSE])
    }, numeric(ncol(z)))
    matrix(sc, nrow = ncol(z), dimnames = list(colnames(z), names(markers)))
  })

  es <- lapply(seq_along(studies), function(i) {
    cond <- .make_condition_factor(studies[[i]]$condition)
    apply(scores[[i]], 2, .hedges_g, group = cond)
  })

  g_mat  <- vapply(es, function(e) e["g", ], numeric(length(markers)))
  se_mat <- vapply(es, function(e) e["se", ], numeric(length(markers)))
  g_mat  <- matrix(g_mat, ncol = length(studies),
                   dimnames = list(names(markers), names(studies)))
  se_mat <- matrix(se_mat, ncol = length(studies), dimnames = dimnames(g_mat))

  out <- .pool_effects(g_mat, se_mat, "cell_type")
  attr(out, "scores") <- scores
  .tag(out, "cell_composition")
}

#' Find consensus co-expression modules across studies
#'
#' Builds a gene-gene Pearson correlation matrix on log-CPM expression in
#' each study, combines them into a consensus (by default the element-wise
#' minimum, as in WGCNA consensus networks, so a link must hold in every
#' study), clusters genes hierarchically, and cuts the tree into modules.
#' Each module is summarised by its eigengene (first principal component) in
#' every study; eigengene differences between cases and controls are pooled
#' across studies with a random-effects model.
#'
#' @param studies A named list of \code{\linkS4class{metaXpressStudy}} objects.
#' @param genes Optional character vector of genes to cluster (e.g. meta-DEGs).
#'   Default \code{NULL} uses the \code{n_genes} most variable genes shared by
#'   all studies.
#' @param n_genes Number of most variable genes to use when \code{genes} is
#'   \code{NULL}. Default 2000.
#' @param consensus How to combine study correlation matrices: \code{"min"}
#'   (conservative, default) or \code{"mean"}.
#' @param cut_height Tree cut height on the distance \eqn{(1 - r)/2}. The
#'   default 0.25 groups genes whose average consensus correlation is at least
#'   about 0.5.
#' @param min_module_size Modules smaller than this are merged into module
#'   \code{"M0"} (unassigned). Default 20.
#'
#' @details Like \code{\link{mx_cell_composition}}, a \code{lib_size}
#'   column in \code{colData} overrides the column sums used for log2-CPM.
#' @param meta_result Optional \code{\linkS4class{metaXpressResult}}. If
#'   given, each module reports the fraction of its genes that are significant
#'   meta-DEGs (\code{meta_padj <= 0.05}) up or down.
#'
#' @return A list with
#'   \describe{
#'     \item{modules}{\code{data.frame} of \code{gene_id}, \code{module},
#'       \code{kME} (mean correlation with the module eigengene across
#'       studies).}
#'     \item{module_stats}{\code{data.frame} with one row per module:
#'       \code{module}, \code{n_genes}, \code{meta_g}, \code{pvalue},
#'       \code{padj}, \code{i_squared}, per-study \code{g_<study>},
#'       \code{hub_genes}, and optionally \code{pct_meta_up},
#'       \code{pct_meta_down}.}
#'     \item{eigengenes}{List of per-study sample x module eigengene
#'       matrices.}
#'   }
#'
#' @references
#' Langfelder, P. and Horvath, S. (2007) Eigengene networks for studying the
#' relationships between co-expression modules. \emph{BMC Systems Biology},
#' \strong{1}, 54. \doi{10.1186/1752-0509-1-54}
#'
#' Langfelder, P. and Horvath, S. (2008) WGCNA: an R package for weighted
#' correlation network analysis. \emph{BMC Bioinformatics}, \strong{9}, 559.
#' \doi{10.1186/1471-2105-9-559}
#'
#' @examples
#' \dontrun{
#'   net <- mx_coexpression_modules(studies, n_genes = 1500,
#'                                  meta_result = meta_result)
#'   net$module_stats
#' }
#'
#' @seealso \code{\link{mx_cell_composition}}, \code{\link{mx_downstream_plot}}
#' @export
mx_coexpression_modules <- function(studies, genes = NULL, n_genes = 2000,
                                    consensus = c("min", "mean"),
                                    cut_height = 0.25, min_module_size = 20,
                                    meta_result = NULL) {
  studies   <- .named_study_list(studies)
  consensus <- match.arg(consensus)

  expr   <- lapply(studies, .log_cpm)
  shared <- Reduce(intersect, lapply(expr, rownames))

  if (is.null(genes)) {
    v <- rowMeans(vapply(expr, function(x) apply(x[shared, , drop = FALSE],
                                                 1, stats::var),
                         numeric(length(shared))))
    v <- v[is.finite(v) & v > 0]
    genes <- names(sort(v, decreasing = TRUE))[seq_len(min(n_genes, length(v)))]
  } else {
    genes <- intersect(genes, shared)
  }
  # Genes must vary in every study to have a defined correlation
  genes <- genes[Reduce(`&`, lapply(expr, function(x)
    apply(x[genes, , drop = FALSE], 1, stats::sd) > 0))]
  if (length(genes) < 2 * min_module_size)
    stop("Too few usable genes (", length(genes), ") for module detection.")

  cors <- lapply(expr, function(x) stats::cor(t(x[genes, , drop = FALSE])))
  cons <- if (consensus == "min") Reduce(pmin, cors) else Reduce(`+`, cors) / length(cors)

  tree <- stats::hclust(stats::as.dist((1 - cons) / 2), method = "average")
  cl   <- stats::cutree(tree, h = cut_height)

  sizes <- table(cl)
  keep  <- names(sizes)[sizes >= min_module_size]
  keep  <- keep[order(-sizes[keep])]
  module <- ifelse(as.character(cl) %in% keep,
                   paste0("M", match(as.character(cl), keep)), "M0")
  names(module) <- genes
  mods <- setdiff(paste0("M", seq_along(keep)), "M0")
  if (length(mods) == 0)
    stop("No module reached min_module_size = ", min_module_size,
         ". Try a larger cut_height or a smaller min_module_size.")

  eig <- lapply(expr, function(x) {
    vapply(mods, function(m) .eigengene(x[names(module)[module == m], ,
                                          drop = FALSE]),
           numeric(ncol(x)))
  })

  kme <- vapply(seq_along(expr), function(i) {
    vapply(genes, function(g) {
      m <- module[[g]]
      if (m == "M0") return(NA_real_)
      stats::cor(expr[[i]][g, ], eig[[i]][, m])
    }, numeric(1))
  }, numeric(length(genes)))
  kme <- stats::setNames(rowMeans(matrix(kme, nrow = length(genes))), genes)

  es <- lapply(seq_along(studies), function(i) {
    cond <- .make_condition_factor(studies[[i]]$condition)
    apply(eig[[i]], 2, .hedges_g, group = cond)
  })
  g_mat  <- matrix(vapply(es, function(e) e["g", ], numeric(length(mods))),
                   ncol = length(studies),
                   dimnames = list(mods, names(studies)))
  se_mat <- matrix(vapply(es, function(e) e["se", ], numeric(length(mods))),
                   ncol = length(studies), dimnames = dimnames(g_mat))

  stats_df <- .pool_effects(g_mat, se_mat, "module")
  stats_df$n_genes <- as.integer(table(module)[stats_df$module])
  stats_df$hub_genes <- vapply(stats_df$module, function(m) {
    k <- kme[module == m]
    paste(names(sort(k, decreasing = TRUE))[seq_len(min(5, length(k)))],
          collapse = ", ")
  }, character(1))

  if (!is.null(meta_result)) {
    mt  <- meta_result@meta_table
    sig <- !is.na(mt$meta_padj) & mt$meta_padj <= 0.05
    up   <- mt$gene_id[sig & mt$meta_log2FC > 0]
    down <- mt$gene_id[sig & mt$meta_log2FC < 0]
    stats_df$pct_meta_up <- vapply(stats_df$module, function(m)
      100 * mean(names(module)[module == m] %in% up), numeric(1))
    stats_df$pct_meta_down <- vapply(stats_df$module, function(m)
      100 * mean(names(module)[module == m] %in% down), numeric(1))
  }

  lead <- c("module", "n_genes")
  stats_df <- stats_df[, c(lead, setdiff(colnames(stats_df), lead))]

  modules_df <- data.frame(gene_id = genes, module = unname(module),
                           kME = unname(kme), stringsAsFactors = FALSE)
  modules_df <- modules_df[order(modules_df$module, -modules_df$kME), ]
  rownames(modules_df) <- NULL

  .tag(list(modules = modules_df,
            module_stats = .tag(stats_df, "coexpression"),
            eigengenes = eig),
       "coexpression")
}

#' Rank drugs by how strongly they reverse the meta-signature
#'
#' Scores each drug perturbation signature (sets of genes the drug pushes up
#' and down) against the meta-analysis ranking. The reversal score is a
#' standardised difference of the mean signed meta z-score of the drug's up
#' and down genes; strongly negative values mean the drug moves genes in the
#' opposite direction to the disease and are candidate therapeutics. The
#' Connectivity Map Kolmogorov-Smirnov score (Lamb et al. 2006) is reported
#' alongside for comparison.
#'
#' @param meta_result A \code{\linkS4class{metaXpressResult}} object.
#' @param signatures A named list of drug signatures, each a list with
#'   character vectors \code{up} and \code{down}, as returned by
#'   \code{\link{mx_read_signatures}}.
#' @param id_col Column of \code{meta_table} whose values match the signature
#'   genes. Default \code{"gene_id"}.
#' @param min_size Minimum number of measured genes in each of the up and
#'   down sets. Default 10.
#'
#' @return A \code{data.frame} sorted by reversal score with columns
#'   \code{drug}, \code{reversal_score} (negative = reverses the disease
#'   signature), \code{pvalue}, \code{padj}, \code{cmap_score}, \code{n_up},
#'   \code{n_down}. P-values assume independent genes and are best used for
#'   ranking rather than as calibrated error rates.
#'
#' @references
#' Lamb, J. et al. (2006) The Connectivity Map: using gene-expression
#' signatures to connect small molecules, genes, and disease.
#' \emph{Science}, \strong{313}(5795), 1929--1935.
#' \doi{10.1126/science.1132939}
#'
#' Subramanian, A. et al. (2017) A next generation Connectivity Map: L1000
#' platform and the first 1,000,000 profiles. \emph{Cell}, \strong{171}(6),
#' 1437--1452. \doi{10.1016/j.cell.2017.10.049}
#'
#' @examples
#' mt <- data.frame(gene_id = paste0("G", 1:300),
#'                  meta_log2FC = c(rep(2, 30), rep(-2, 30), rnorm(240, sd = 0.2)),
#'                  meta_pvalue = c(rep(1e-8, 60), runif(240)),
#'                  meta_padj = 1, i_squared = 0, n_studies = 3L)
#' res <- new("metaXpressResult", meta_table = mt, method = "fisher",
#'            n_studies = 3L)
#' sigs <- list(reverser = list(up = paste0("G", 31:60), down = paste0("G", 1:30)),
#'              mimic    = list(up = paste0("G", 1:30), down = paste0("G", 31:60)))
#' mx_drug_connectivity(res, sigs)
#'
#' @seealso \code{\link{mx_read_signatures}}, \code{\link{mx_downstream_plot}}
#' @export
mx_drug_connectivity <- function(meta_result, signatures, id_col = "gene_id",
                                 min_size = 10) {
  z <- .gene_statistic(meta_result, id_col, "z")
  ranking <- names(sort(z, decreasing = TRUE))
  mu <- mean(z)
  sigma2 <- stats::var(z)

  rows <- lapply(names(signatures), function(d) {
    up   <- intersect(signatures[[d]]$up, names(z))
    down <- intersect(signatures[[d]]$down, names(z))
    if (length(up) < min_size || length(down) < min_size) return(NULL)
    diff <- mean(z[up]) - mean(z[down])
    stat <- diff / sqrt(sigma2 * (1 / length(up) + 1 / length(down)))
    ks_up   <- .cmap_ks(match(up, ranking), length(ranking))
    ks_down <- .cmap_ks(match(down, ranking), length(ranking))
    cmap <- if (sign(ks_up) == sign(ks_down)) 0 else ks_up - ks_down
    data.frame(drug = d, reversal_score = stat,
               pvalue = 2 * stats::pnorm(-abs(stat)), cmap_score = cmap,
               n_up = length(up), n_down = length(down),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out))
    stop("No signature has >= ", min_size, " measured genes in both its up ",
         "and down sets. Check that 'id_col' matches the signature genes.")

  out$padj <- stats::p.adjust(out$pvalue, method = "BH")
  out <- out[order(out$reversal_score), c("drug", "reversal_score", "pvalue",
                                          "padj", "cmap_score", "n_up",
                                          "n_down")]
  rownames(out) <- NULL
  .tag(out, "drug_connectivity")
}

#' Read up/down drug signatures from a gene set library
#'
#' Parses a GMT file or an Enrichr gene set library text file and pairs sets
#' named \code{"<drug> up"} and \code{"<drug> down"} (case-insensitive,
#' separated by a space, underscore or hyphen) into signatures for
#' \code{\link{mx_drug_connectivity}}. Works with Enrichr's
#' \code{LINCS_L1000_Chem_Pert_Consensus_Sigs} library.
#'
#' @param path Path to a GMT or Enrichr library file (tab-separated: set name,
#'   description, genes...).
#'
#' @return A named list; each element is a list with character vectors
#'   \code{up} and \code{down}. Drugs lacking either direction are dropped.
#'
#' @examples
#' f <- tempfile(fileext = ".gmt")
#' writeLines(c("drugA up\t\tG1\tG2\tG3", "drugA down\t\tG4\tG5"), f)
#' mx_read_signatures(f)
#'
#' @seealso \code{\link{mx_drug_connectivity}}
#' @export
mx_read_signatures <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  parts <- strsplit(lines, "\t", fixed = TRUE)
  set_names <- vapply(parts, `[`, character(1), 1)
  genes <- lapply(parts, function(p) {
    g <- p[-(1:2)]
    unique(sub(",.*$", "", g[nzchar(g)]))
  })

  pattern   <- "^(.*?)[ _-](up|down)$"
  is_dir    <- grepl(pattern, set_names, ignore.case = TRUE)
  drug      <- sub(pattern, "\\1", set_names[is_dir], ignore.case = TRUE)
  direction <- tolower(sub(pattern, "\\2", set_names[is_dir],
                           ignore.case = TRUE))
  genes     <- genes[is_dir]

  drugs <- unique(drug)
  sigs <- lapply(drugs, function(d) {
    list(up   = unlist(genes[drug == d & direction == "up"]),
         down = unlist(genes[drug == d & direction == "down"]))
  })
  names(sigs) <- drugs
  sigs[vapply(sigs, function(s) length(s$up) > 0 && length(s$down) > 0,
              logical(1))]
}

#' Plot downstream biology results
#'
#' Draws a ranked bar chart for any result of \code{\link{mx_tf_activity}},
#' \code{\link{mx_cell_composition}}, \code{\link{mx_coexpression_modules}}
#' or \code{\link{mx_drug_connectivity}}. Bars show the score (TF activity,
#' pooled Hedges' g or reversal score); faded bars are not significant.
#'
#' @param x A result from one of the downstream functions.
#' @param top Number of entries to show. For TF activity and drug
#'   connectivity the top entries in each direction are shown (\code{top / 2}
#'   each). Default 20.
#' @param padj_threshold Significance threshold for bar shading. Default 0.05.
#' @param title Optional plot title.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' \dontrun{
#'   mx_downstream_plot(mx_tf_activity(meta_result, collectri))
#' }
#'
#' @importFrom ggplot2 geom_col coord_flip scale_fill_manual scale_alpha_manual
#'   geom_hline
#' @export
mx_downstream_plot <- function(x, top = 20, padj_threshold = 0.05,
                               title = NULL) {
  type <- attr(x, "mx_type")
  if (is.null(type))
    stop("'x' must be a result of a metaXpress downstream function")
  if (type == "coexpression" && !is.data.frame(x)) x <- x$module_stats

  spec <- switch(type,
    tf_activity       = list(label = "source", score = "score",
                             y = "TF activity (ULM t-value)",
                             t = "Transcription factor activity",
                             two_sided = TRUE),
    cell_composition  = list(label = "cell_type", score = "meta_g",
                             y = "Pooled Hedges' g (case vs control)",
                             t = "Cell-type composition shift",
                             two_sided = FALSE),
    coexpression      = list(label = "module", score = "meta_g",
                             y = "Pooled eigengene Hedges' g",
                             t = "Consensus co-expression modules",
                             two_sided = FALSE),
    drug_connectivity = list(label = "drug", score = "reversal_score",
                             y = "Reversal score (negative = reverses disease)",
                             t = "Drug signature connectivity",
                             two_sided = TRUE)
  )

  df <- x[!is.na(x[[spec$score]]), , drop = FALSE]
  if (spec$two_sided) {
    k  <- ceiling(top / 2)
    df <- df[order(df[[spec$score]]), , drop = FALSE]
    idx <- unique(c(utils::head(seq_len(nrow(df)), k),
                    utils::tail(seq_len(nrow(df)), k)))
    df <- df[idx, , drop = FALSE]
  } else {
    df <- utils::head(df[order(df$pvalue), , drop = FALSE], top)
  }

  plot_df <- data.frame(
    label = factor(df[[spec$label]],
                   levels = df[[spec$label]][order(df[[spec$score]])]),
    score = df[[spec$score]],
    direction = ifelse(df[[spec$score]] > 0, "Up", "Down"),
    significant = ifelse(!is.na(df$padj) & df$padj <= padj_threshold,
                         "yes", "no")
  )

  ggplot2::ggplot(plot_df, ggplot2::aes(x = label, y = score,
                                        fill = direction,
                                        alpha = significant)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(Up = "#E41A1C", Down = "#377EB8"),
                               guide = "none") +
    ggplot2::scale_alpha_manual(
      values = c(yes = 1, no = 0.35), limits = c("yes", "no"), drop = FALSE,
      name = paste0("padj <= ", padj_threshold),
      guide = ggplot2::guide_legend(override.aes = list(fill = "grey30"))) +
    ggplot2::labs(title = if (is.null(title)) spec$t else title,
                  x = NULL, y = spec$y) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}

# ============================================================================
# Internal helpers
# ============================================================================

.tag <- function(x, type) {
  attr(x, "mx_type") <- type
  x
}

# Named vector of a gene-level statistic, one value per identifier
.gene_statistic <- function(meta_result, id_col, stat) {
  if (!is(meta_result, "metaXpressResult"))
    stop("'meta_result' must be a metaXpressResult object")
  mt <- meta_result@meta_table
  if (!id_col %in% colnames(mt))
    stop("meta_table has no column '", id_col, "'")

  val <- if (stat == "z") {
    p <- pmax(mt$meta_pvalue, .Machine$double.xmin)
    sign(mt$meta_log2FC) * stats::qnorm(p / 2, lower.tail = FALSE)
  } else {
    mt$meta_log2FC
  }
  ids <- as.character(mt[[id_col]])
  ok  <- !is.na(val) & !is.na(ids) & nzchar(ids)
  val <- val[ok]; ids <- ids[ok]
  # Keep the strongest value for identifiers that appear more than once
  o <- order(-abs(val))
  val <- val[o]; ids <- ids[o]
  keep <- !duplicated(ids)
  stats::setNames(val[keep], ids[keep])
}

.standardize_network <- function(network) {
  net <- as.data.frame(network, stringsAsFactors = FALSE)
  if (all(c("source_genesymbol", "target_genesymbol") %in% colnames(net))) {
    mor <- if ("is_inhibition" %in% colnames(net))
      ifelse(as.logical(net$is_inhibition) &
               !as.logical(net$is_stimulation), -1, 1)
    else 1
    net <- data.frame(source = net$source_genesymbol,
                      target = net$target_genesymbol, mor = mor,
                      stringsAsFactors = FALSE)
  }
  if (!all(c("source", "target") %in% colnames(net)))
    stop("'network' needs columns 'source' and 'target' (and optionally 'mor')")
  if (!"mor" %in% colnames(net)) net$mor <- 1
  net <- net[!is.na(net$source) & !is.na(net$target), c("source", "target", "mor")]
  net[!duplicated(net[, c("source", "target")]), ]
}

.standardize_markers <- function(markers) {
  if (is.data.frame(markers)) {
    if (!all(c("gene", "cell_type") %in% colnames(markers)))
      stop("'markers' data.frame needs columns 'gene' and 'cell_type'")
    markers <- split(as.character(markers$gene), markers$cell_type)
  }
  if (!is.list(markers) || is.null(names(markers)))
    stop("'markers' must be a named list or a data.frame")
  lapply(markers, function(g) unique(as.character(g)))
}

.named_study_list <- function(studies) {
  if (!is.list(studies) || length(studies) < 2 ||
      !all(vapply(studies, is, logical(1), "metaXpressStudy")))
    stop("'studies' must be a list of at least two metaXpressStudy objects")
  if (is.null(names(studies)))
    names(studies) <- vapply(studies, function(s) s@accession, character(1))
  studies
}

# log2-CPM of a study. A 'lib_size' colData column, if present, overrides the
# column sums so that gene-subset studies are normalised to their true depth.
.log_cpm <- function(study) {
  counts <- .counts(study)
  lib <- if ("lib_size" %in% colnames(SummarizedExperiment::colData(study)))
    study$lib_size else colSums(counts)
  log2(t(t(counts) / lib * 1e6) + 1)
}

.row_z <- function(x) {
  s <- apply(x, 1, stats::sd)
  x <- x[s > 0, , drop = FALSE]
  (x - rowMeans(x)) / s[s > 0]
}

# Hedges' g (case - control) and its standard error for one numeric vector
.hedges_g <- function(x, group) {
  lv <- levels(group)
  a <- x[group == lv[1]]; b <- x[group == lv[2]]
  a <- a[!is.na(a)]; b <- b[!is.na(b)]
  n1 <- length(a); n2 <- length(b)
  if (n1 < 2 || n2 < 2) return(c(g = NA_real_, se = NA_real_))
  sp <- sqrt(((n1 - 1) * stats::var(a) + (n2 - 1) * stats::var(b)) /
               (n1 + n2 - 2))
  if (!is.finite(sp) || sp == 0) return(c(g = NA_real_, se = NA_real_))
  j <- 1 - 3 / (4 * (n1 + n2) - 9)
  g <- j * (mean(b) - mean(a)) / sp
  c(g = g, se = sqrt((n1 + n2) / (n1 * n2) + g^2 / (2 * (n1 + n2))))
}

# Random-effects pooling of a feature x study matrix of effects
.pool_effects <- function(g_mat, se_mat, label) {
  k <- rowSums(!is.na(g_mat) & !is.na(se_mat))
  if (any(k < 2))
    message("Dropping ", sum(k < 2), " ", label, "(s) scored in fewer than ",
            "two studies: ", paste(rownames(g_mat)[k < 2], collapse = ", "))
  g_mat  <- g_mat[k >= 2, , drop = FALSE]
  se_mat <- se_mat[k >= 2, , drop = FALSE]
  se_mat[!is.finite(se_mat) | se_mat <= 0] <- NA_real_
  het <- .heterogeneity(g_mat, se_mat, tau_method = "DL")
  re  <- .meta_re(g_mat, se_mat, het$tau2, hksj = FALSE, alpha = 0.05)
  out <- data.frame(
    label     = rownames(g_mat),
    meta_g    = re$est,
    pvalue    = re$p,
    padj      = stats::p.adjust(re$p, method = "BH"),
    i_squared = het$I_sq,
    n_studies = rowSums(!is.na(g_mat) & !is.na(se_mat)),
    stringsAsFactors = FALSE
  )
  names(out)[1] <- label
  per_study <- as.data.frame(g_mat)
  names(per_study) <- paste0("g_", colnames(g_mat))
  out <- cbind(out, per_study)
  out <- out[order(out$pvalue), ]
  rownames(out) <- NULL
  out
}

# Module eigengene: first PC of row-scaled expression, signed to follow the
# module's average expression
.eigengene <- function(x) {
  z  <- t(scale(t(x)))
  z[!is.finite(z)] <- 0
  pc <- stats::prcomp(t(z), center = FALSE, rank. = 1)$x[, 1]
  if (stats::cor(pc, colMeans(z)) < 0) -pc else pc
}

# Connectivity Map Kolmogorov-Smirnov enrichment of positions in a ranking
.cmap_ks <- function(pos, n) {
  pos <- sort(pos[!is.na(pos)])
  t <- length(pos)
  if (t == 0) return(0)
  j <- seq_len(t)
  a <- max(j / t - pos / n)
  b <- max(pos / n - (j - 1) / t)
  if (a > b) a else -b
}
