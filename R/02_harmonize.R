#' @include AllClasses.R
NULL

# ============================================================================
# Module 2: Normalization & Harmonization
# ============================================================================

#' Reannotate gene IDs to a common namespace
#'
#' Converts gene identifiers in one or more studies to a unified namespace
#' (Ensembl gene ID, gene symbol, or Entrez ID) using \pkg{AnnotationDbi}.
#' Genes that cannot be mapped are dropped with a warning. Harmonising gene
#' identifiers across studies is a prerequisite for any cross-study comparison.
#'
#' @param studies A named list of \code{\linkS4class{metaXpressStudy}} objects.
#' @param org Character scalar. Target organism. Currently only
#'   \code{"Homo sapiens"} is supported. Default: \code{"Homo sapiens"}.
#' @param target_id Character scalar. Target gene identifier type.
#'   One of \code{"SYMBOL"} (default), \code{"ENSEMBL"}, or
#'   \code{"ENTREZID"}.
#'
#' @return The input list of studies with \code{counts} rownames translated
#'   to the target namespace.
#'
#' @references
#' Pagès, H. et al. (2024) AnnotationDbi: manipulation of SQLite-based
#' annotations in Bioconductor. R package version 1.66.0.
#' \url{https://bioconductor.org/packages/AnnotationDbi}
#'
#' @examples
#' \dontrun{
#'   studies <- mx_reannotate(studies, org = "Homo sapiens",
#'                             target_id = "SYMBOL")
#' }
#'
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object
#'   controlling parallelization. Default: \code{BiocParallel::bpparam()}.
#'   Set to \code{BiocParallel::SerialParam()} for sequential execution.
#'
#' @seealso \code{\link{mx_normalize}}, \code{\link{mx_align_genes}}
#'
#' @export
mx_reannotate <- function(studies, org = "Homo sapiens",
                           target_id = c("SYMBOL", "ENSEMBL", "ENTREZID"),
                           BPPARAM = BiocParallel::bpparam()) {
  if (!is.list(studies))
    stop("'studies' must be a list of metaXpressStudy objects")
  target_id <- match.arg(target_id)

  org_db <- .get_org_db(org)

  message("Reannotating gene IDs to ", target_id, " for ",
          length(studies), " studies...")

  BiocParallel::bplapply(studies, function(s) {
    .reannotate_study(s, org_db = org_db, target_id = target_id)
  }, BPPARAM = BPPARAM)
}

#' Normalize counts within a single study
#'
#' Applies one of five normalization methods to the raw count matrix of a
#' \code{metaXpressStudy} object.
#'
#' @param study A \code{\linkS4class{metaXpressStudy}} object with raw counts.
#' @param method Character scalar. Normalization method. One of
#'   \code{"TMM"} (default), \code{"VST"}, \code{"CPM"}, \code{"TPM"},
#'   or \code{"quantile"}.
#'
#' @return The input \code{study} with the normalized values added as a new
#'   assay (\code{"logCPM_TMM"}, \code{"vst"}, \code{"cpm"}, \code{"tpm"}
#'   or \code{"log2_quantile"}). The raw \code{"counts"} assay is never
#'   modified, because \code{\link{mx_de}} requires raw counts.
#'
#' @details
#' Normalized values are for visualization, clustering, co-expression and
#' cell-composition analyses. Differential expression always works from raw
#' counts with each engine's own normalization.
#' \itemize{
#'   \item \strong{TMM}: log2 CPM with trimmed-mean-of-M-values library
#'     sizes (Robinson & Oshlack 2010), via \pkg{edgeR} (prior count 2).
#'   \item \strong{VST}: Variance-stabilizing transformation (Anders & Huber
#'     2010), implemented in \pkg{DESeq2}. Produces approximately
#'     homoscedastic values suitable for clustering and visualization.
#'   \item \strong{CPM}: Counts per million via \pkg{edgeR}. Simple library-size
#'     scaling without variance stabilization.
#'   \item \strong{TPM}: Transcripts per million (Li et al. 2010). Requires
#'     gene lengths in \code{metadata}. Preferred when comparing expression
#'     levels across genes.
#'   \item \strong{quantile}: Quantile normalization (Bolstad et al. 2003),
#'     implemented in \pkg{limma}. Forces identical distributions across
#'     samples; operates on log2(CPM + 1) values.
#' }
#'
#' @references
#' Robinson, M.D. & Oshlack, A. (2010) A scaling normalization method for
#' differential expression analysis of RNA-seq data. \emph{Genome Biology},
#' \strong{11}(3), R25. \doi{10.1186/gb-2010-11-3-r25}
#'
#' Anders, S. & Huber, W. (2010) Differential expression analysis for sequence
#' count data. \emph{Genome Biology}, \strong{11}(10), R106.
#' \doi{10.1186/gb-2010-11-10-r106}
#'
#' Bolstad, B.M. et al. (2003) A comparison of normalization methods for high
#' density oligonucleotide array data based on variance and bias.
#' \emph{Bioinformatics}, \strong{19}(2), 185--193.
#' \doi{10.1093/bioinformatics/19.2.185}
#'
#' @examples
#' counts <- matrix(rpois(600, 50), nrow = 100,
#'                  dimnames = list(paste0("G", 1:100), paste0("S", 1:6)))
#' meta <- data.frame(condition = rep(c("control", "case"), each = 3))
#' study <- metaXpressStudy(counts, meta, accession = "local")
#' study <- mx_normalize(study, method = "TMM")
#' SummarizedExperiment::assayNames(study)
#'
#' @seealso \code{\link{mx_reannotate}}, \code{\link{mx_remove_batch}}
#'
#' @export
mx_normalize <- function(study,
                          method = c("TMM", "VST", "CPM", "TPM", "quantile")) {
  if (!is(study, "metaXpressStudy"))
    stop("'study' must be a metaXpressStudy object")
  method <- match.arg(method)

  message("  Normalizing study ", study@accession, " using ", method, "...")

  normalized <- switch(method,
    TMM      = .normalize_tmm(study),
    VST      = .normalize_vst(study),
    CPM      = .normalize_cpm(study),
    TPM      = .normalize_tpm(study),
    quantile = .normalize_quantile(study)
  )
  assay_name <- switch(method, TMM = "logCPM_TMM", VST = "vst", CPM = "cpm",
                       TPM = "tpm", quantile = "log2_quantile")
  dimnames(normalized) <- dimnames(study)
  SummarizedExperiment::assay(study, assay_name) <- normalized
  study
}

#' Correct for library type bias (polyA vs rRNA-depleted)
#'
#' Adjusts for systematic expression differences between polyA-selected and
#' rRNA-depleted RNA-seq libraries with a per-gene ratio correction in the
#' spirit of Bush et al. (2017). The ratio is estimated from
#' \emph{control} samples only, so disease effects do not leak into it.
#'
#' Per-gene multiplicative scaling cancels out of within-study log fold
#' changes, so this correction is \strong{not needed} for the standard
#' per-study DE + meta-analysis workflow. It matters only when expression
#' levels are compared across studies (pooled heatmaps, co-expression,
#' cell-composition scores). The corrected values are stored in a new
#' \code{"cpm_libtype_corrected"} assay; raw counts are untouched.
#'
#' @param studies A named list of \code{\linkS4class{metaXpressStudy}} objects.
#'   Studies must have a \code{library_type} column in \code{metadata} with
#'   values \code{"polyA"} or \code{"rRNA_depleted"}.
#'
#' @return The input list of studies, restricted to their common genes, with
#'   an added \code{"cpm_libtype_corrected"} assay in rRNA-depleted studies
#'   (and plain CPM in polyA studies).
#'
#' @references
#' Bush, S.J. et al. (2017) Integration of quantitated expression estimates
#' from polyA-selected and rRNA-depleted RNA-seq libraries. \emph{BMC
#' Bioinformatics}, \strong{18}, 301. \doi{10.1186/s12859-017-1714-9}
#'
#' @examples
#' \dontrun{
#'   studies <- mx_correct_library_type(studies)
#' }
#'
#' @export
mx_correct_library_type <- function(studies) {
  if (!is.list(studies))
    stop("'studies' must be a list of metaXpressStudy objects")

  lib_types <- vapply(studies, function(s) {
    if ("library_type" %in% colnames(SummarizedExperiment::colData(s)))
      unique(s$library_type)[1]
    else
      NA_character_
  }, character(1))

  n_unknown <- sum(is.na(lib_types))
  if (n_unknown > 0)
    warning(n_unknown, " study/studies lack 'library_type' metadata; ",
            "library type correction skipped for these.")

  if (length(unique(na.omit(lib_types))) < 2) {
    message("All studies use the same library type; no correction needed.")
    return(studies)
  }

  # Bush et al. (2017) ratio-based correction:
  # Scale rRNA-depleted libraries to match polyA by a per-gene scaling factor
  # computed from genes stably expressed across both library types.
  polya_idx    <- which(lib_types == "polyA")
  rrna_idx     <- which(lib_types == "rRNA_depleted")

  if (length(polya_idx) == 0 || length(rrna_idx) == 0) {
    message("Only one library type found in studies with known type; ",
            "no correction applied.")
    return(studies)
  }

  # Align genes first
  aligned <- mx_align_genes(studies)

  .require_package("edgeR", "library-type correction")
  cpm_of <- function(s) edgeR::cpm(edgeR::normLibSizes(
    edgeR::DGEList(counts = .counts(s))), log = FALSE)
  control_mean <- function(s) {
    ctrl <- .make_condition_factor(s$condition)
    rowMeans(cpm_of(s)[, ctrl == levels(ctrl)[1], drop = FALSE])
  }

  # Mean control-sample CPM per gene for each library type
  mean_polya <- rowMeans(vapply(aligned[polya_idx], control_mean,
                                numeric(nrow(aligned[[1]]))))
  mean_rrna  <- rowMeans(vapply(aligned[rrna_idx], control_mean,
                                numeric(nrow(aligned[[1]]))))
  ratio <- mean_polya / pmax(mean_rrna, 0.01)

  for (i in seq_along(aligned)) {
    vals <- cpm_of(aligned[[i]])
    if (i %in% rrna_idx) vals <- sweep(vals, 1, ratio, "*")
    SummarizedExperiment::assay(aligned[[i]], "cpm_libtype_corrected") <- vals
  }

  message(sprintf("Library type correction applied to %d rRNA-depleted studies.",
                  length(rrna_idx)))
  aligned
}

#' Remove between-study batch effects for pooled expression analyses
#'
#' Combines all studies into one expression matrix and removes study-level
#' batch effects while preserving the biological \code{condition} effect.
#'
#' @param studies A named list of \code{\linkS4class{metaXpressStudy}} objects.
#' @param method Character scalar. \code{"ComBat-seq"} (default; negative
#'   binomial model on counts), \code{"ComBat"} (empirical Bayes on
#'   log2 CPM) or \code{"limma"} (\code{removeBatchEffect} on log2 CPM).
#' @param covariate Name of a \code{colData} column whose effect must be
#'   preserved. Default: \code{"condition"}. Use \code{NULL} only if you
#'   understand that condition effects confounded with study will be removed.
#'
#' @return The input list of studies, restricted to their common genes, with
#'   the corrected values in a new assay: \code{"combat_seq_counts"} for
#'   ComBat-seq or \code{"log2cpm_batch_corrected"} otherwise. The raw
#'   \code{"counts"} assay is untouched.
#'
#' @details
#' \strong{Batch correction is not required, and should not be used, before
#' per-study DE}: \code{\link{mx_de}} analyses each study separately, so
#' between-study batch effects never enter the per-study fold changes, and
#' feeding batch-corrected values into count models invalidates their
#' variance assumptions. Use this function only for pooled visualizations
#' (PCA, heatmaps), consensus co-expression or other analyses on a merged
#' matrix. Batch is defined as study of origin.
#'
#' @references
#' Johnson, W.E., Li, C. & Rabinovic, A. (2007) Adjusting batch effects in
#' microarray expression data using empirical Bayes methods.
#' \emph{Biostatistics}, \strong{8}(1), 118--127.
#' \doi{10.1093/biostatistics/kxj037}
#'
#' Zhang, Y., Parmigiani, G. & Johnson, W.E. (2020) ComBat-seq: batch effect
#' adjustment for RNA-seq count data. \emph{NAR Genomics and
#' Bioinformatics}, \strong{2}(3), lqaa078. \doi{10.1093/nargab/lqaa078}
#'
#' Nygaard, V., Rødland, E.A. & Hovig, E. (2016) Methods that remove batch
#' effects while retaining group differences may lead to exaggerated
#' confidence in downstream analyses. \emph{Biostatistics}, \strong{17}(1),
#' 29--39. \doi{10.1093/biostatistics/kxv027}
#'
#' @examples
#' \dontrun{
#'   pooled <- mx_remove_batch(studies, method = "limma")
#'   SummarizedExperiment::assay(pooled[[1]], "log2cpm_batch_corrected")
#' }
#'
#' @seealso \code{\link{mx_normalize}}
#'
#' @export
mx_remove_batch <- function(studies,
                            method = c("ComBat-seq", "ComBat", "limma"),
                            covariate = "condition") {
  if (!is.list(studies))
    stop("'studies' must be a list of metaXpressStudy objects")
  if (identical(method, "harmony"))
    stop("method = 'harmony' was removed: harmony corrects low-dimensional ",
         "embeddings and cannot return a gene-level matrix.")
  method <- match.arg(method)

  message("Removing between-study batch effects using ", method, "...")

  aligned <- mx_align_genes(studies)
  combined_counts <- do.call(cbind, lapply(aligned, .counts))
  study_ids <- if (is.null(names(aligned))) seq_along(aligned)
               else names(aligned)
  batch <- factor(rep(study_ids, times = vapply(aligned, ncol, integer(1))))

  group <- NULL
  if (!is.null(covariate)) {
    group <- unlist(lapply(aligned, function(s) {
      if (!covariate %in% colnames(SummarizedExperiment::colData(s)))
        stop("Covariate '", covariate, "' missing in study ", s@accession)
      as.character(.make_condition_factor(s[[covariate]]))
    }), use.names = FALSE)
    group <- factor(group)
    if (any(rowSums(table(batch, group) > 0) < 2))
      warning("Some studies contain only one '", covariate, "' level; ",
              "batch and biology are partly confounded.")
  }

  log2cpm <- function(x) edgeR::cpm(edgeR::DGEList(counts = x), log = TRUE,
                                    prior.count = 2)

  corrected_combined <- switch(method,
    `ComBat-seq` = {
      .require_package("sva", "ComBat-seq")
      sva::ComBat_seq(combined_counts, batch = batch, group = group)
    },
    ComBat = {
      .require_package("sva", "ComBat")
      mod <- if (is.null(group)) NULL else stats::model.matrix(~ group)
      sva::ComBat(log2cpm(combined_counts), batch = batch, mod = mod)
    },
    limma = {
      .require_package("limma", "removeBatchEffect")
      design <- if (is.null(group)) matrix(1, length(batch), 1)
                else stats::model.matrix(~ group)
      limma::removeBatchEffect(log2cpm(combined_counts), batch = batch,
                               design = design)
    }
  )
  assay_name <- if (method == "ComBat-seq") "combat_seq_counts"
                else "log2cpm_batch_corrected"

  col_idx <- 0
  for (i in seq_along(aligned)) {
    n <- ncol(aligned[[i]])
    vals <- corrected_combined[, (col_idx + 1):(col_idx + n), drop = FALSE]
    dimnames(vals) <- dimnames(aligned[[i]])
    SummarizedExperiment::assay(aligned[[i]], assay_name) <- vals
    col_idx <- col_idx + n
  }

  aligned
}

#' Restrict all studies to a common set of genes
#'
#' Finds the intersection of gene IDs across all studies and returns studies
#' whose count matrices are restricted to that common set. Gene IDs are taken
#' from the rownames of each study's count matrix.
#'
#' @param studies A named list of \code{\linkS4class{metaXpressStudy}} objects.
#'
#' @return The input list of studies, each with \code{counts} restricted to the
#'   intersection of all gene sets.
#'
#' @examples
#' \dontrun{
#'   studies <- mx_align_genes(studies)
#' }
#'
#' @seealso \code{\link{mx_reannotate}}
#'
#' @export
mx_align_genes <- function(studies) {
  if (!is.list(studies) || length(studies) == 0)
    stop("'studies' must be a non-empty list of metaXpressStudy objects")

  gene_sets <- lapply(studies, rownames)
  common_genes <- Reduce(intersect, gene_sets)

  if (length(common_genes) == 0)
    stop("No common genes found across studies. ",
         "Run mx_reannotate() first to standardize gene ID namespaces.")

  pct <- length(common_genes) / length(gene_sets[[1]]) * 100
  message(sprintf("Aligning to %d common genes (%.1f%% of first study).",
                   length(common_genes), pct))

  lapply(studies, function(s) s[common_genes, ])
}

# ============================================================================
# Internal helpers
# ============================================================================

.get_org_db <- function(org) {
  switch(org,
    "Homo sapiens" = {
      if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
        stop("Package 'org.Hs.eg.db' is required for human gene annotation")
      org.Hs.eg.db::org.Hs.eg.db
    },
    stop("Organism '", org, "' is not yet supported. ",
         "Currently supported: 'Homo sapiens'")
  )
}

.reannotate_study <- function(study, org_db, target_id) {
  current_ids  <- rownames(study)
  keytype      <- .detect_id_type(current_ids)

  if (keytype == target_id)
    return(study)

  mapped <- AnnotationDbi::mapIds(
    org_db,
    keys    = current_ids,
    column  = target_id,
    keytype = keytype,
    multiVals = "first"
  )

  valid <- !is.na(mapped)
  n_lost <- sum(!valid)
  if (n_lost > 0)
    warning(sprintf("Study %s: %d genes (%.1f%%) could not be mapped to %s and were dropped.",
                    study@accession, n_lost,
                    n_lost / length(valid) * 100, target_id))

  study <- study[valid, ]
  rownames(study) <- mapped[valid]
  study
}

.detect_id_type <- function(ids) {
  sample_ids <- head(ids, 20)
  if (all(grepl("^ENSG[0-9]", sample_ids)))     return("ENSEMBL")
  if (all(grepl("^[0-9]+$", sample_ids)))        return("ENTREZID")
  return("SYMBOL")
}

.normalize_tmm <- function(study) {
  dge <- edgeR::DGEList(counts = .counts(study))
  dge <- edgeR::normLibSizes(dge, method = "TMM")
  edgeR::cpm(dge, log = TRUE, prior.count = 2)
}

.normalize_vst <- function(study) {
  dds  <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(.counts(study)),
    colData   = .samples(study),
    design    = ~ condition
  )
  as.matrix(SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE)))
}

.normalize_cpm <- function(study) {
  dge <- edgeR::DGEList(counts = .counts(study))
  edgeR::cpm(dge, log = FALSE)
}

.normalize_tpm <- function(study) {
  gene_info <- SummarizedExperiment::rowData(study)
  if (!"gene_length" %in% colnames(gene_info))
    stop("TPM normalization requires a 'gene_length' column in rowData(study)")

  gene_lengths <- gene_info$gene_length

  # RPK: reads per kilobase
  rpk <- sweep(.counts(study), 1, gene_lengths / 1000, "/")
  # Scale each sample to per-million
  col_sums <- colSums(rpk) / 1e6
  sweep(rpk, 2, col_sums, "/")
}

.normalize_quantile <- function(study) {
  log_cpm <- log2(edgeR::cpm(edgeR::DGEList(counts = .counts(study))) + 1)
  limma::normalizeBetweenArrays(log_cpm, method = "quantile")
}
