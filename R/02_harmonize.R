#' @include AllClasses.R
NULL

# ============================================================================
# Module 2: Normalization & Harmonization
# ============================================================================

#' Reannotate gene IDs to a common namespace
#'
#' Converts the gene identifiers of one or more studies to a single
#' namespace (gene symbol, Ensembl gene ID or Entrez ID) so that studies can
#' be compared gene by gene.
#'
#' Processing steps for each study:
#' \enumerate{
#'   \item Ensembl version suffixes are removed (\code{ENSG00000141510.17}
#'     becomes \code{ENSG00000141510}), as are Ensembl \code{_PAR_Y}
#'     pseudo-autosomal duplicates.
#'   \item The source identifier type is detected from \emph{all} IDs by
#'     majority (Ensembl, Entrez or symbol), or set with \code{from_id}.
#'   \item IDs are mapped with the organism's \code{org.*.eg.db} package, or
#'     with a user-supplied \code{mapping} table (any organism, offline).
#'   \item Unmapped IDs are dropped and reported.
#'   \item When several source IDs map to the same target ID, they are
#'     collapsed with \code{duplicates}: \code{"max_mean"} keeps the most
#'     highly expressed source row (default, as in WGCNA's
#'     \code{collapseRows}); \code{"sum"} adds counts (appropriate when the
#'     source rows are disjoint features of one gene); \code{"drop"}
#'     removes all ambiguous genes.
#' }
#' The original identifier of every retained row is stored in
#' \code{rowData(study)$source_id}.
#'
#' @param studies A named list of \code{\linkS4class{metaXpressStudy}}
#'   objects, or a single study.
#' @param org Organism: \code{"Homo sapiens"} (default; \pkg{org.Hs.eg.db}),
#'   \code{"Mus musculus"} (\pkg{org.Mm.eg.db}) or \code{"Rattus norvegicus"}
#'   (\pkg{org.Rn.eg.db}). Ignored when \code{mapping} is given.
#' @param target_id Target identifier: \code{"SYMBOL"} (default),
#'   \code{"ENSEMBL"} or \code{"ENTREZID"}.
#' @param from_id Source identifier type, or \code{NULL} (default) to detect
#'   it per study.
#' @param mapping Optional \code{data.frame} with columns \code{from} and
#'   \code{to} giving the ID mapping to use instead of an \code{org.*.eg.db}
#'   package.
#' @param duplicates How to collapse several source IDs mapping to one
#'   target: \code{"max_mean"} (default), \code{"sum"} or \code{"drop"}.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}}. Default:
#'   \code{BiocParallel::SerialParam()}.
#'
#' @return The studies with translated rownames (same type as the input: a
#'   list or a single study). A per-study mapping summary is printed.
#'
#' @references
#' Pages, H. et al. AnnotationDbi: manipulation of SQLite-based annotations
#' in Bioconductor. \url{https://bioconductor.org/packages/AnnotationDbi}
#'
#' Miller, J.A. et al. (2011) Strategies for aggregating gene expression data:
#' the collapseRows R function. \emph{BMC Bioinformatics}, \strong{12}, 322.
#' \doi{10.1186/1471-2105-12-322}
#'
#' @examples
#' counts <- matrix(c(10, 20, 30, 5, 1, 2, 3, 4, 7, 7, 7, 7),
#'                  nrow = 3, byrow = TRUE,
#'                  dimnames = list(c("ENSG00000141510.17", "ENSG00000012048.23",
#'                                    "ENSG00000999999.1"),
#'                                  paste0("S", 1:4)))
#' meta  <- data.frame(condition = c("control", "control", "case", "case"))
#' study <- metaXpressStudy(counts, meta, accession = "local")
#' map   <- data.frame(from = c("ENSG00000141510", "ENSG00000012048"),
#'                     to   = c("TP53", "BRCA1"))
#' out <- mx_reannotate(study, mapping = map)
#' rownames(out)
#' SummarizedExperiment::rowData(out)$source_id
#'
#' @seealso \code{\link{mx_align_genes}}, \code{\link{mx_gene_sets}}
#'
#' @export
mx_reannotate <- function(studies, org = "Homo sapiens",
                          target_id = c("SYMBOL", "ENSEMBL", "ENTREZID"),
                          from_id = NULL,
                          mapping = NULL,
                          duplicates = c("max_mean", "sum", "drop"),
                          BPPARAM = BiocParallel::SerialParam()) {
  single <- is(studies, "metaXpressStudy")
  if (single) studies <- list(studies)
  .validate_study_list(studies)
  target_id  <- match.arg(target_id)
  duplicates <- match.arg(duplicates)
  if (!is.null(from_id))
    from_id <- match.arg(from_id, c("SYMBOL", "ENSEMBL", "ENTREZID"))

  org_db <- NULL
  if (is.null(mapping)) {
    org_db <- .get_org_db(org)
  } else if (!is.data.frame(mapping) ||
             !all(c("from", "to") %in% names(mapping))) {
    stop("'mapping' must be a data.frame with columns 'from' and 'to'")
  }

  out <- BiocParallel::bplapply(studies, .reannotate_study,
                                org_db = org_db, target_id = target_id,
                                from_id = from_id, mapping = mapping,
                                duplicates = duplicates, BPPARAM = BPPARAM)
  for (s in out) message(attr(s, "reannotation"))
  out <- lapply(out, function(s) { attr(s, "reannotation") <- NULL; s })
  names(out) <- names(studies)
  if (single) out[[1]] else out
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
  pkg <- switch(org,
                "Homo sapiens"      = "org.Hs.eg.db",
                "Mus musculus"      = "org.Mm.eg.db",
                "Rattus norvegicus" = "org.Rn.eg.db",
                stop("Organism '", org, "' has no built-in annotation. ",
                     "Supply 'mapping' (data.frame with from/to columns)."))
  .require_package(pkg, paste("gene annotation for", org))
  .require_package("AnnotationDbi", "gene annotation")
  getExportedValue(pkg, pkg)
}

# Remove Ensembl version suffixes and PAR_Y tags
.clean_gene_ids <- function(ids) {
  ids <- sub("_PAR_Y$", "", ids)
  ens <- grepl("^ENS[A-Z]*G[0-9]+\\.[0-9]+$", ids)
  ids[ens] <- sub("\\.[0-9]+$", "", ids[ens])
  ids
}

# Detect identifier type by majority over all IDs
.detect_id_type <- function(ids) {
  ids <- ids[!is.na(ids)]
  if (length(ids) == 0) return("SYMBOL")
  frac_ens <- mean(grepl("^ENS[A-Z]*G[0-9]+(\\.[0-9]+)?$", ids))
  frac_num <- mean(grepl("^[0-9]+$", ids))
  if (frac_ens > 0.5) return("ENSEMBL")
  if (frac_num > 0.5) return("ENTREZID")
  "SYMBOL"
}

# Collapse rows sharing a target ID. Returns list(counts, keep_source).
.collapse_rows <- function(counts, target, source, method) {
  dup <- target %in% target[duplicated(target)]
  if (!any(dup))
    return(list(counts = counts, target = target, source = source))
  if (method == "drop") {
    return(list(counts = counts[!dup, , drop = FALSE], target = target[!dup],
                source = source[!dup]))
  }
  if (method == "sum") {
    summed <- rowsum(counts, target, reorder = FALSE)
    src <- tapply(source, target, paste, collapse = ";")
    return(list(counts = summed, target = rownames(summed),
                source = unname(src[rownames(summed)])))
  }
  # max_mean: keep the source row with the highest mean count
  ord  <- order(target, -rowMeans(counts))
  keep <- ord[!duplicated(target[ord])]
  keep <- sort(keep)
  list(counts = counts[keep, , drop = FALSE], target = target[keep],
       source = source[keep])
}

.reannotate_study <- function(study, org_db, target_id, from_id = NULL,
                              mapping = NULL, duplicates = "max_mean") {
  orig_ids  <- rownames(study)
  clean_ids <- .clean_gene_ids(orig_ids)
  src_type  <- if (is.null(from_id)) .detect_id_type(clean_ids) else from_id

  if (!is.null(mapping)) {
    m <- mapping[!is.na(mapping$from) & !is.na(mapping$to), ]
    m$from <- .clean_gene_ids(as.character(m$from))
    mapped <- as.character(m$to)[match(clean_ids, m$from)]
  } else if (src_type == target_id) {
    mapped <- clean_ids
  } else {
    mapped <- suppressMessages(AnnotationDbi::mapIds(
      org_db, keys = clean_ids, column = target_id, keytype = src_type,
      multiVals = "first"))
    mapped <- unname(as.character(mapped))
  }
  mapped[!is.na(mapped) & !nzchar(mapped)] <- NA_character_

  ok <- !is.na(mapped)
  counts <- .counts(study)[ok, , drop = FALSE]
  col <- .collapse_rows(counts, mapped[ok], orig_ids[ok], duplicates)
  n_dup_rows <- sum(ok) - nrow(col$counts)

  rd <- S4Vectors::DataFrame(source_id = col$source, row.names = col$target)
  dimnames(col$counts) <- list(col$target, colnames(study))
  new_study <- metaXpressStudy(
    counts    = col$counts,
    metadata  = .samples(study),
    accession = study@accession,
    organism  = study@organism,
    qc_score  = study@qc_score,
    gene_info = rd)
  S4Vectors::metadata(new_study) <- S4Vectors::metadata(study)
  attr(new_study, "reannotation") <- sprintf(
    "  %s: %s -> %s | %d IDs, %d unmapped (%.1f%%), %d rows collapsed (%s), %d genes",
    study@accession, src_type, target_id, length(orig_ids), sum(!ok),
    100 * mean(!ok), n_dup_rows, duplicates, nrow(col$counts))
  new_study
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
