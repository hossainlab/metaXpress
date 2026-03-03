#' metaXpress: End-to-End Bulk RNA-seq Meta-Analysis
#'
#' @description
#' \pkg{metaXpress} provides a unified, reproducible pipeline for integrating
#' multiple bulk RNA-seq studies into a single meta-analysis. Starting from raw
#' GEO accession IDs, the package handles every step: data ingestion and QC,
#' gene-ID harmonisation, library-type correction and batch removal, per-study
#' differential expression, meta-analysis statistics (p-value combination and
#' effect-size models), missing-gene imputation, pathway-level enrichment, and
#' production-quality reporting.
#'
#' The pipeline is built around two S4 classes that flow through every step:
#' \code{\linkS4class{metaXpressStudy}} (one object per study, output of
#' ingestion) and \code{\linkS4class{metaXpressResult}} (the combined
#' meta-analysis result).  All public functions are prefixed \code{mx_} and
#' parallelised via \pkg{BiocParallel}.
#'
#' @section S4 Classes:
#' \describe{
#'   \item{\code{\linkS4class{metaXpressStudy}}}{
#'     Holds a single RNA-seq study: raw count matrix, sample metadata,
#'     accession ID, organism, QC score, and per-study DE result.
#'     Created by \code{\link{mx_fetch_geo}}, \code{\link{mx_fetch_sra}},
#'     or \code{\link{mx_load_local}}.}
#'   \item{\code{\linkS4class{metaXpressResult}}}{
#'     Holds the combined meta-analysis output: gene-level statistics
#'     (meta log2FC, combined p-value, FDR, I\eqn{^2}), heterogeneity
#'     table, and (optionally) pathway enrichment results.
#'     Created by \code{\link{mx_meta}}.}
#' }
#'
#' @section Module 1 — Data Ingestion and QC:
#' \describe{
#'   \item{\code{\link{mx_fetch_geo}}}{
#'     Download count matrices and sample metadata from GEO for one or more
#'     accession IDs; applies the 10-point QC checklist automatically.}
#'   \item{\code{\link{mx_fetch_sra}}}{
#'     Download raw RNA-seq data for SRA project IDs (requires
#'     \pkg{sratools}).}
#'   \item{\code{\link{mx_load_local}}}{
#'     Construct \code{metaXpressStudy} objects from local count-matrix files
#'     or \code{SummarizedExperiment} objects.}
#'   \item{\code{\link{mx_qc_study}}}{
#'     Score a study against the 10-point QC checklist derived from
#'     Heberle et al. (2025).}
#'   \item{\code{\link{mx_cluster_samples}}}{
#'     Automatically assign samples to case/control groups using metadata
#'     clustering.}
#'   \item{\code{\link{mx_filter_studies}}}{
#'     Retain only studies whose \code{qc_score} meets a minimum threshold.}
#' }
#'
#' @section Module 2 — Normalisation and Harmonisation:
#' \describe{
#'   \item{\code{\link{mx_reannotate}}}{
#'     Convert gene identifiers to a unified namespace (SYMBOL, ENSEMBL, or
#'     ENTREZID) using \pkg{AnnotationDbi}.}
#'   \item{\code{\link{mx_normalize}}}{
#'     Normalise within-study counts by TMM, VST, CPM, TPM, or quantile
#'     methods.}
#'   \item{\code{\link{mx_correct_library_type}}}{
#'     Adjust for systematic bias between polyA-selected and rRNA-depleted
#'     libraries.}
#'   \item{\code{\link{mx_remove_batch}}}{
#'     Remove cross-study batch effects via ComBat-seq, ComBat, limma, or
#'     Harmony.}
#'   \item{\code{\link{mx_align_genes}}}{
#'     Restrict all studies to their common gene universe.}
#' }
#'
#' @section Module 3 — Per-Study Differential Expression:
#' \describe{
#'   \item{\code{\link{mx_de}}}{
#'     Run DE on a single \code{metaXpressStudy} using DESeq2, edgeR, or
#'     limma-voom.}
#'   \item{\code{\link{mx_de_all}}}{
#'     Apply \code{\link{mx_de}} to every study in a list, in parallel via
#'     \pkg{BiocParallel}.}
#'   \item{\code{\link{mx_de_summary}}}{
#'     Tabulate significant DEG counts across studies at user-specified
#'     thresholds.}
#' }
#'
#' @section Module 4 — Meta-Analysis Statistics:
#' \describe{
#'   \item{\code{\link{mx_meta}}}{
#'     Combine per-study DE results using Fisher, Stouffer, inverse-normal,
#'     fixed-effects, random-effects (DerSimonian-Laird), or AWmeta adaptive
#'     weighting.  Returns a \code{\linkS4class{metaXpressResult}}.}
#'   \item{\code{\link{mx_heterogeneity}}}{
#'     Compute per-gene Cochran's Q, I\eqn{^2}, \eqn{\tau^2}, and
#'     heterogeneity p-values.}
#'   \item{\code{\link{mx_sensitivity}}}{
#'     Leave-one-out sensitivity analysis: rerun meta-analysis excluding each
#'     study in turn to assess result stability.}
#'   \item{\code{\link{mx_study_overview}}}{
#'     Summarise each study's contribution to the meta-analysis.}
#' }
#'
#' @section Module 5 — Missing Gene Handling:
#' \describe{
#'   \item{\code{\link{mx_missing_summary}}}{
#'     Report gene-by-study coverage as a binary presence/absence matrix.}
#'   \item{\code{\link{mx_filter_coverage}}}{
#'     Retain only genes detected in at least \code{min_studies} studies.}
#'   \item{\code{\link{mx_impute}}}{
#'     Impute missing per-study statistics (zero, mean, median, or
#'     Bayesian shrinkage).}
#' }
#'
#' @section Module 6 — Pathway Meta-Analysis:
#' \describe{
#'   \item{\code{\link{mx_pathway_meta}}}{
#'     Run cross-study ORA or GSEA against MSigDB gene sets (Hallmarks,
#'     KEGG, Reactome, GO) using \pkg{clusterProfiler}.}
#'   \item{\code{\link{mx_pathway_consensus}}}{
#'     Identify pathways enriched across a majority of individual studies.}
#'   \item{\code{\link{mx_pathway_dedup}}}{
#'     Remove redundant pathways via kappa-coefficient overlap clustering.}
#'   \item{\code{\link{mx_pathway_heatmap}}}{
#'     Heatmap of enrichment scores or \eqn{-\log_{10}}(padj) for top
#'     pathways.}
#' }
#'
#' @section Module 7 — Visualisation:
#' \describe{
#'   \item{\code{\link{mx_volcano}}}{
#'     Volcano plot of meta log2FC vs. combined \eqn{-\log_{10}}(padj) with
#'     optional gene labels.}
#'   \item{\code{\link{mx_forest}}}{
#'     Forest plot showing per-study log2FC ± 95\% CI and the pooled
#'     meta-estimate for a single gene.}
#'   \item{\code{\link{mx_heatmap}}}{
#'     Heatmap of log2FC values for the top \emph{n} meta-significant genes
#'     across all studies.}
#'   \item{\code{\link{mx_heterogeneity_plot}}}{
#'     Histogram of the I\eqn{^2} distribution across all tested genes.}
#'   \item{\code{\link{mx_upset}}}{
#'     UpSet plot of DEG overlap across studies.}
#'   \item{\code{\link{mx_study_overview}}}{
#'     Multi-panel QC summary: scores, sample sizes, sequencing depths, and
#'     gene detection rates.}
#' }
#'
#' @section Module 8 — Reporting and Export:
#' \describe{
#'   \item{\code{\link{mx_report}}}{
#'     Render a fully reproducible HTML or PDF analysis report from a
#'     parameterised R Markdown template.}
#'   \item{\code{\link{mx_export}}}{
#'     Write the gene-level meta-analysis table to CSV, Excel, or RDS.}
#'   \item{\code{\link{mx_session_info}}}{
#'     Capture R session information, timestamp, and package version for
#'     reproducibility records.}
#' }
#'
#' @section Example data:
#' \describe{
#'   \item{\code{\link{metaXpress_example}}}{
#'     A minimal simulated dataset: three \code{metaXpressStudy} objects
#'     (500/500/480 genes, 6/6/8 samples) for use in examples and tests.
#'     Load with \code{data(metaXpress_example)}.}
#' }
#'
#' @section Getting started:
#' The fastest way to explore the package is the quickstart vignette:
#' \code{vignette("quickstart", package = "metaXpress")}.
#' For the complete GEO-to-report workflow see
#' \code{vignette("full_workflow", package = "metaXpress")}.
#'
#' @references
#' Heberle, H. et al. (2025) A comprehensive framework for quality control
#' and meta-analysis of bulk RNA-seq data. \emph{Alzheimer's & Dementia},
#' \strong{21}(1), e70025. \doi{10.1002/alz.70025}
#'
#' Fisher, R.A. (1932) \emph{Statistical Methods for Research Workers}.
#' 4th edn. Edinburgh: Oliver & Boyd.
#'
#' DerSimonian, R. & Laird, N. (1986) Meta-analysis in clinical trials.
#' \emph{Controlled Clinical Trials}, \strong{7}(3), 177--188.
#' \doi{10.1016/0197-2456(86)90046-2}
#'
#' Cochran, W.G. (1954) The combination of estimates from different
#' experiments. \emph{Biometrics}, \strong{10}(1), 101--129.
#' \doi{10.2307/3001666}
#'
#' Higgins, J.P.T. & Thompson, S.G. (2002) Quantifying heterogeneity in a
#' meta-analysis. \emph{Statistics in Medicine}, \strong{21}(11),
#' 1539--1558. \doi{10.1002/sim.1186}
#'
#' Love, M.I., Huber, W. & Anders, S. (2014) Moderated estimation of fold
#' change and dispersion for RNA-seq data with DESeq2. \emph{Genome
#' Biology}, \strong{15}(12), 550. \doi{10.1186/s13059-014-0550-8}
#'
#' Robinson, M.D., McCarthy, D.J. & Smyth, G.K. (2010) edgeR: a
#' Bioconductor package for differential expression analysis of digital
#' gene expression data. \emph{Bioinformatics}, \strong{26}(1), 139--140.
#' \doi{10.1093/bioinformatics/btp616}
#'
#' Ritchie, M.E. et al. (2015) limma powers differential expression
#' analyses for RNA-sequencing and microarray studies. \emph{Nucleic Acids
#' Research}, \strong{43}(7), e47. \doi{10.1093/nar/gkv007}
#'
#' Johnson, W.E., Li, C. & Rabinovic, A. (2007) Adjusting batch effects in
#' microarray expression data using empirical Bayes methods.
#' \emph{Biostatistics}, \strong{8}(1), 118--127.
#' \doi{10.1093/biostatistics/kxj037}
#'
#' @author Md. Jubayer Hossain \email{contact.jubayerhossain@@gmail.com}
#'
#' @importFrom stats median na.omit
#' @importFrom utils head
"_PACKAGE"

# Suppress R CMD check NOTEs for ggplot2 non-standard evaluation column names
utils::globalVariables(c(
  "log2FC", "study", "ci_lo", "ci_hi", "x", "y",
  "i_squared", "meta_log2FC", "neg_log10_padj", "sig", "gene_id"
))
