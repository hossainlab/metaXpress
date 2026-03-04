#' @include AllClasses.R
NULL

# ============================================================================
# Module 7: Visualization
# ============================================================================

#' Meta-analysis volcano plot
#'
#' Creates a volcano plot from a \code{metaXpressResult} object, displaying
#' meta-analysis effect sizes (log2 fold-change) against statistical
#' significance (-log10 meta_padj). The top differentially expressed genes
#' are labelled. Volcano plots are a standard visualization for differential
#' expression results (Li 2012).
#'
#' @param meta_result A \code{\linkS4class{metaXpressResult}} object.
#' @param padj_threshold Numeric. Adjusted p-value threshold for colouring
#'   significant genes. Default: \code{0.05}.
#' @param lfc_threshold Numeric. Log2 fold-change threshold. Default: \code{1}.
#' @param label_top Integer. Number of top significant genes to label.
#'   Default: \code{10}.
#' @param title Character scalar. Plot title. Default: \code{NULL} (auto).
#'
#' @return A \code{ggplot2} object.
#'
#' @references
#' Li, W. (2012) Volcano plots in analyzing differential expressions with
#' mRNA microarrays. \emph{Journal of Bioinformatics and Computational
#' Biology}, \strong{10}(6), 1231003. \doi{10.1142/S0219720012310038}
#'
#' @examples
#' \dontrun{
#'   p <- mx_volcano(meta_result)
#'   print(p)
#' }
#'
#' @seealso \code{\link{mx_forest}}, \code{\link{mx_heatmap}}
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_hline
#'   scale_colour_manual labs theme_bw theme element_text
#' @importFrom ggrepel geom_label_repel
#' @export
mx_volcano <- function(meta_result, padj_threshold = 0.05,
                        lfc_threshold = 1, label_top = 10,
                        title = NULL) {
  if (!is(meta_result, "metaXpressResult"))
    stop("'meta_result' must be a metaXpressResult object")

  mt <- meta_result@meta_table
  if (nrow(mt) == 0)
    stop("meta_table is empty; run mx_meta() first")

  mt$neg_log10_padj <- -log10(mt$meta_padj)
  mt$sig <- ifelse(
    !is.na(mt$meta_padj) & mt$meta_padj <= padj_threshold &
      abs(mt$meta_log2FC) >= lfc_threshold,
    ifelse(mt$meta_log2FC > 0, "Up", "Down"),
    "NS"
  )

  label_df <- mt[mt$sig != "NS", ]
  if (nrow(label_df) > label_top)
    label_df <- label_df[order(label_df$meta_padj)[seq_len(label_top)], ]

  colours <- c(Up = "#E41A1C", Down = "#377EB8", NS = "grey70")

  auto_title <- if (is.null(title))
    paste0("Meta-analysis volcano (", meta_result@method, ", n=",
           meta_result@n_studies, " studies)")
  else
    title

  ggplot2::ggplot(mt, ggplot2::aes(x = meta_log2FC, y = neg_log10_padj,
                                     colour = sig)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
                        linetype = "dashed", colour = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(padj_threshold),
                        linetype = "dashed", colour = "grey40") +
    ggplot2::scale_colour_manual(values = colours, name = "Direction") +
    ggrepel::geom_label_repel(
      data = label_df,
      ggplot2::aes(label = gene_id),
      size = 3, max.overlaps = 20, show.legend = FALSE
    ) +
    ggplot2::labs(title = auto_title,
                  x = "Meta log\u2082 fold change",
                  y = "-log\u2081\u2080 (meta_padj)") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}

#' Forest plot for a single gene
#'
#' Displays per-study log2 fold-change estimates with 95\% confidence intervals
#' and a pooled estimate diamond. I-squared heterogeneity is annotated. Forest
#' plots are the standard visualization for meta-analysis results, showing both
#' individual study contributions and the pooled summary (Lewis & Clarke 2001).
#'
#' @param gene Character scalar. Gene identifier to plot.
#' @param de_results A named list of per-study DE \code{data.frame}s (from
#'   \code{\link{mx_de_all}} result slots), or a list of
#'   \code{\linkS4class{metaXpressStudy}} objects.
#' @param meta_result Optional. A \code{\linkS4class{metaXpressResult}} object.
#'   When provided, the pooled estimate diamond and I-squared are taken from
#'   the meta-analysis result.
#'
#' @return A \code{ggplot2} object.
#'
#' @references
#' Lewis, S. & Clarke, M. (2001) Forest plots: trying to see the wood and the
#' trees. \emph{BMJ}, \strong{322}(7300), 1479--1480.
#' \doi{10.1136/bmj.322.7300.1479}
#'
#' @examples
#' \dontrun{
#'   p <- mx_forest("TP53", de_results, meta_result)
#'   print(p)
#' }
#'
#' @seealso \code{\link{mx_volcano}}, \code{\link{mx_heterogeneity_plot}}
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbarh geom_vline
#'   geom_polygon annotate scale_y_discrete labs theme_bw theme
#' @export
mx_forest <- function(gene, de_results, meta_result = NULL) {
  if (!is.character(gene) || length(gene) != 1)
    stop("'gene' must be a single character string")

  if (is.list(de_results) && all(vapply(de_results, is, logical(1),
                                         "metaXpressStudy")))
    de_results <- lapply(de_results, function(s) s@de_result)

  study_names <- names(de_results)
  if (is.null(study_names))
    study_names <- paste0("Study_", seq_along(de_results))

  rows <- lapply(seq_along(de_results), function(i) {
    d   <- de_results[[i]]
    row <- d[d$gene_id == gene, , drop = FALSE]
    if (nrow(row) == 0)
      return(data.frame(study = study_names[i], log2FC = NA_real_,
                        se = NA_real_, stringsAsFactors = FALSE))
    se <- if ("lfcSE" %in% colnames(row) && !is.na(row$lfcSE[1]))
      row$lfcSE[1]
    else
      abs(row$log2FC[1]) / qnorm(1 - row$pvalue[1] / 2 + .Machine$double.eps)
    data.frame(study = study_names[i], log2FC = row$log2FC[1], se = se,
               stringsAsFactors = FALSE)
  })

  plot_df <- do.call(rbind, rows)
  plot_df <- plot_df[!is.na(plot_df$log2FC), ]

  if (nrow(plot_df) == 0)
    stop("Gene '", gene, "' not found in any study.")

  plot_df$ci_lo <- plot_df$log2FC - 1.96 * plot_df$se
  plot_df$ci_hi <- plot_df$log2FC + 1.96 * plot_df$se
  plot_df$study <- factor(plot_df$study, levels = rev(plot_df$study))

  pooled_lfc  <- NA_real_
  i_sq_label  <- ""

  if (!is.null(meta_result) && is(meta_result, "metaXpressResult")) {
    idx <- which(meta_result@meta_table$gene_id == gene)
    if (length(idx) == 1) {
      pooled_lfc <- meta_result@meta_table$meta_log2FC[idx]
      i_sq <- meta_result@meta_table$i_squared[idx]
      if (!is.na(i_sq))
        i_sq_label <- sprintf("I\u00b2 = %.1f%%", i_sq)
    }
  }

  p <- ggplot2::ggplot(plot_df,
                        ggplot2::aes(x = log2FC, y = study)) +
    ggplot2::geom_point(size = 3, colour = "#377EB8") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = ci_lo, xmax = ci_hi),
      width = 0.2, colour = "#377EB8", orientation = "y"
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey40") +
    ggplot2::labs(title = paste0("Forest plot: ", gene),
                  x = "log\u2082 fold change (95% CI)",
                  y = NULL,
                  caption = i_sq_label) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  if (!is.na(pooled_lfc)) {
    half_h <- 0.4
    ymax   <- nrow(plot_df)
    diamond_df <- data.frame(
      x = c(pooled_lfc - 0.1, pooled_lfc, pooled_lfc + 0.1, pooled_lfc),
      y = c(0, half_h, 0, -half_h) + ymax + 1.2
    )
    p <- p +
      ggplot2::geom_polygon(data = diamond_df,
                             ggplot2::aes(x = x, y = y),
                             fill = "#E41A1C", inherit.aes = FALSE)
  }

  p
}

#' Heatmap of top DEGs across studies
#'
#' Displays log2 fold-change values for the top differentially expressed genes
#' across all studies as a heatmap.
#'
#' @param meta_result A \code{\linkS4class{metaXpressResult}} object.
#' @param studies A named list of \code{\linkS4class{metaXpressStudy}} objects
#'   (used to extract per-study log2FC if needed).
#' @param top_n Integer. Number of top genes (by \code{meta_padj}) to display.
#'   Default: \code{50}.
#'
#' @return A \code{ComplexHeatmap::Heatmap} object.
#'
#' @examples
#' \dontrun{
#'   mx_heatmap(meta_result, studies, top_n = 30)
#' }
#'
#' @seealso \code{\link{mx_volcano}}, \code{\link{mx_upset}}
#'
#' @export
mx_heatmap <- function(meta_result, studies, top_n = 50) {
  if (!is(meta_result, "metaXpressResult"))
    stop("'meta_result' must be a metaXpressResult object")

  mt <- meta_result@meta_table
  mt_ordered <- mt[order(mt$meta_padj, na.last = TRUE), ]
  top_genes  <- head(mt_ordered$gene_id, top_n)

  de_list <- lapply(studies, function(s) s@de_result)
  k       <- length(de_list)

  lfc_mat <- matrix(NA_real_, nrow = length(top_genes), ncol = k,
                     dimnames = list(top_genes, names(de_list)))

  for (i in seq_len(k)) {
    d   <- de_list[[i]]
    idx <- match(top_genes, d$gene_id)
    lfc_mat[, i] <- d$log2FC[idx]
  }

  max_abs <- max(abs(lfc_mat), na.rm = TRUE)
  col_fun <- circlize::colorRamp2(
    c(-max_abs, 0, max_abs),
    c("#377EB8", "white", "#E41A1C")
  )

  ComplexHeatmap::Heatmap(
    lfc_mat,
    name            = "log2FC",
    col             = col_fun,
    cluster_rows    = TRUE,
    cluster_columns = FALSE,
    show_row_names  = TRUE,
    show_column_names = TRUE,
    row_names_gp    = grid::gpar(fontsize = 8),
    column_title    = paste0("Top ", top_n, " DEGs (by meta_padj)")
  )
}

#' UpSet plot of DEG overlap across studies
#'
#' Shows the intersection structure of significant DEGs across studies using
#' an UpSet plot.
#'
#' @param de_results A named list of per-study DE \code{data.frame}s.
#' @param padj_threshold Numeric. Adjusted p-value threshold. Default: \code{0.05}.
#' @param lfc_threshold Numeric. Log2 fold-change threshold. Default: \code{1}.
#'
#' @return A \code{ggplot2} object (requires \pkg{ggupset} or
#'   \code{ComplexHeatmap::UpSet}).
#'
#' @examples
#' \dontrun{
#'   mx_upset(de_results)
#' }
#'
#' @export
mx_upset <- function(de_results, padj_threshold = 0.05, lfc_threshold = 1) {
  stop("mx_upset() is not yet implemented.")
}

#' Study characteristics overview plot
#'
#' Creates a multi-panel summary of study-level QC metrics including
#' QC scores, sample sizes, sequencing depths, and gene detection rates.
#'
#' @param studies A named list of \code{\linkS4class{metaXpressStudy}} objects.
#'
#' @return A \code{ggplot2} object (or list of plots).
#'
#' @examples
#' \dontrun{
#'   mx_study_overview(studies)
#' }
#'
#' @export
mx_study_overview <- function(studies) {
  stop("mx_study_overview() is not yet implemented.")
}

#' I-squared distribution plot
#'
#' Displays the distribution of I-squared heterogeneity values across all
#' genes in the meta-analysis, with reference lines at 25\% (low
#' heterogeneity) and 75\% (high heterogeneity). These thresholds follow the
#' classification proposed by Higgins et al. (2003): I-squared values of
#' 25\%, 50\%, and 75\% correspond to low, moderate, and high heterogeneity,
#' respectively.
#'
#' @param meta_result A \code{\linkS4class{metaXpressResult}} object.
#'
#' @return A \code{ggplot2} object.
#'
#' @references
#' Higgins, J.P.T. et al. (2003) Measuring inconsistency in meta-analyses.
#' \emph{BMJ}, \strong{327}(7414), 557--560. \doi{10.1136/bmj.327.7414.557}
#'
#' @examples
#' \dontrun{
#'   mx_heterogeneity_plot(meta_result)
#' }
#'
#' @seealso \code{\link{mx_heterogeneity}}, \code{\link{mx_forest}}
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline labs theme_bw
#' @export
mx_heterogeneity_plot <- function(meta_result) {
  if (!is(meta_result, "metaXpressResult"))
    stop("'meta_result' must be a metaXpressResult object")

  i_sq_vals <- meta_result@meta_table$i_squared
  df <- data.frame(i_squared = i_sq_vals[!is.na(i_sq_vals)])

  ggplot2::ggplot(df, ggplot2::aes(x = i_squared)) +
    ggplot2::geom_histogram(binwidth = 5, fill = "#377EB8", colour = "white",
                             alpha = 0.8) +
    ggplot2::geom_vline(xintercept = c(25, 75), linetype = "dashed",
                        colour = c("#E6AB02", "#E41A1C")) +
    ggplot2::labs(title = "I\u00b2 heterogeneity distribution",
                  x     = "I\u00b2 (%)",
                  y     = "Number of genes",
                  caption = "Dashed lines: 25% (low) and 75% (high) thresholds") +
    ggplot2::theme_bw()
}
