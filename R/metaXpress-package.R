#' metaXpress: End-to-End Bulk RNA-seq Meta-Analysis
#'
#' metaXpress provides a unified, reproducible pipeline for integrating
#' multiple bulk RNA-seq studies — from raw GEO data retrieval through
#' meta-analysis statistics to pathway-level interpretation and reporting.
#'
#' @section Main workflow:
#' \enumerate{
#'   \item \code{\link{mx_fetch_geo}} / \code{\link{mx_load_local}} — data ingestion
#'   \item \code{\link{mx_qc_study}} / \code{\link{mx_filter_studies}} — quality control
#'   \item \code{\link{mx_reannotate}} / \code{\link{mx_remove_batch}} — harmonization
#'   \item \code{\link{mx_de_all}} — per-study differential expression
#'   \item \code{\link{mx_meta}} — meta-analysis statistics
#'   \item \code{\link{mx_pathway_meta}} — pathway enrichment
#'   \item \code{\link{mx_volcano}} / \code{\link{mx_forest}} — visualization
#'   \item \code{\link{mx_report}} — reproducible reporting
#' }
#'
#' @importFrom stats median na.omit
#' @importFrom utils head
"_PACKAGE"

# Suppress R CMD check NOTEs for ggplot2 non-standard evaluation column names
utils::globalVariables(c(
  "log2FC", "study", "ci_lo", "ci_hi", "x", "y",
  "i_squared", "meta_log2FC", "neg_log10_padj", "sig", "gene_id"
))
