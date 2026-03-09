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
#' identifiers across studies is a prerequisite for any cross-study comparison
#' (Stokholm et al. 2024).
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
#' Stokholm, A., Rabaglino, M.B. & Kadarmideen, H.N. (2024) A protocol for
#' cross-platform meta-analysis of microarray and RNA-seq data.
#' \emph{Current Protocols}, \strong{4}(12), e70046.
#' \doi{10.1002/cpz1.70046}
#'
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
#' @return The input \code{study} with \code{counts} replaced by normalized
#'   values.
#'
#' @details
#' \itemize{
#'   \item \strong{TMM}: Trimmed mean of M-values (Robinson & Oshlack 2010),
#'     implemented in \pkg{edgeR}. Recommended for between-sample normalization
#'     of count data.
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
#'     samples; operates on log1p-transformed counts.
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
#' \dontrun{
#'   study <- mx_normalize(study, method = "TMM")
#' }
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

  study@counts <- normalized
  study
}

#' Correct for library type bias (polyA vs rRNA-depleted)
#'
#' Adjusts for systematic expression differences between polyA-selected and
#' rRNA-depleted RNA-seq libraries using the ratio-based correction described
#' by Bush et al. (2017).
#'
#' @param studies A named list of \code{\linkS4class{metaXpressStudy}} objects.
#'   Studies must have a \code{library_type} column in \code{metadata} with
#'   values \code{"polyA"} or \code{"rRNA_depleted"}.
#'
#' @return The input list of studies with library type bias corrected.
#'
#' @references
#' Bush, S.J. et al. (2017) Cross-species inference of long non-coding RNAs
#' greatly expands the ruminant transcriptome. \emph{Genetics Selection
#' Evolution}, \strong{50}(1), 20. \doi{10.1186/s12859-017-1714-9}
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
    if ("library_type" %in% colnames(s@metadata))
      unique(s@metadata$library_type)[1]
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

  stop("mx_correct_library_type() implementation in progress.")
}

#' Remove batch effects across studies
#'
#' Applies one of four batch effect removal methods to harmonize expression
#' levels across studies before meta-analysis.
#'
#' @param studies A named list of \code{\linkS4class{metaXpressStudy}} objects.
#' @param method Character scalar. Batch removal method. One of
#'   \code{"ComBat-seq"} (default), \code{"ComBat"}, \code{"limma"},
#'   or \code{"harmony"}.
#'
#' @return The input list of studies with batch effects removed from
#'   \code{counts}.
#'
#' @details
#' Batch is defined as study of origin. \code{"ComBat-seq"} operates on raw
#' counts and is preferred for count data. \code{"ComBat"} and
#' \code{"limma"} operate on log-transformed values.
#'
#' @references
#' Johnson, W.E., Li, C. & Rabinovic, A. (2007) Adjusting batch effects in
#' microarray expression data using empirical Bayes methods.
#' \emph{Biostatistics}, \strong{8}(1), 118--127.
#' \doi{10.1093/biostatistics/kxj037}
#'
#' Zhang, Y. et al. (2020) ComBat-seq: batch effect adjustment for RNA-seq
#' count data. \emph{NAR Genomics and Bioinformatics}, \strong{2}(3),
#' lqaa078. \doi{10.1093/nargab/lqaa078}
#'
#' Ritchie, M.E. et al. (2015) limma powers differential expression analyses
#' for RNA-sequencing and microarray studies. \emph{Nucleic Acids Research},
#' \strong{43}(7), e47. \doi{10.1093/nar/gkv007}
#'
#' @examples
#' \dontrun{
#'   studies <- mx_remove_batch(studies, method = "ComBat-seq")
#' }
#'
#' @seealso \code{\link{mx_normalize}}
#'
#' @export
mx_remove_batch <- function(studies,
                             method = c("ComBat-seq", "ComBat",
                                        "limma", "harmony")) {
  if (!is.list(studies))
    stop("'studies' must be a list of metaXpressStudy objects")
  method <- match.arg(method)

  message("Removing batch effects using ", method, "...")

  aligned <- mx_align_genes(studies)
  combined_counts <- do.call(cbind, lapply(aligned, function(s) s@counts))
  batch <- rep(seq_along(aligned),
               times = vapply(aligned, function(s) ncol(s@counts), integer(1)))

  corrected_combined <- switch(method,
    `ComBat-seq` = sva::ComBat_seq(combined_counts, batch = batch),
    ComBat       = {
      log_counts <- log1p(combined_counts)
      sva::ComBat(log_counts, batch = batch)
    },
    limma        = limma::removeBatchEffect(log1p(combined_counts),
                                            batch = batch),
    harmony      = stop("harmony batch correction is not yet implemented")
  )

  col_idx <- 0
  for (i in seq_along(aligned)) {
    n <- ncol(aligned[[i]]@counts)
    aligned[[i]]@counts <- corrected_combined[, (col_idx + 1):(col_idx + n),
                                               drop = FALSE]
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

  gene_sets <- lapply(studies, function(s) rownames(s@counts))
  common_genes <- Reduce(intersect, gene_sets)

  if (length(common_genes) == 0)
    stop("No common genes found across studies. ",
         "Run mx_reannotate() first to standardize gene ID namespaces.")

  pct <- length(common_genes) / length(gene_sets[[1]]) * 100
  message(sprintf("Aligning to %d common genes (%.1f%% of first study).",
                   length(common_genes), pct))

  lapply(studies, function(s) {
    s@counts <- s@counts[common_genes, , drop = FALSE]
    s
  })
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
  current_ids  <- rownames(study@counts)
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

  study@counts <- study@counts[valid, , drop = FALSE]
  rownames(study@counts) <- mapped[valid]
  study
}

.detect_id_type <- function(ids) {
  sample_ids <- head(ids, 20)
  if (all(grepl("^ENSG[0-9]", sample_ids)))     return("ENSEMBL")
  if (all(grepl("^[0-9]+$", sample_ids)))        return("ENTREZID")
  return("SYMBOL")
}

.normalize_tmm <- function(study) {
  dge <- edgeR::DGEList(counts = study@counts)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  edgeR::cpm(dge, log = FALSE)
}

.normalize_vst <- function(study) {
  cond <- factor(study@metadata$condition)
  dds  <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(study@counts),
    colData   = study@metadata,
    design    = ~ condition
  )
  as.matrix(SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE)))
}

.normalize_cpm <- function(study) {
  dge <- edgeR::DGEList(counts = study@counts)
  edgeR::cpm(dge, log = FALSE)
}

.normalize_tpm <- function(study) {
  if (!"gene_length" %in% colnames(study@metadata))
    stop("TPM normalization requires a 'gene_length' column in metadata")
  stop("TPM normalization is not yet implemented.")
}

.normalize_quantile <- function(study) {
  log_counts <- log1p(study@counts)
  limma::normalizeBetweenArrays(log_counts, method = "quantile")
}
