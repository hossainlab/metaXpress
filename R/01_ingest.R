#' @include AllClasses.R
NULL

# ============================================================================
# Module 1: Data Ingestion & QC
# ============================================================================

#' Fetch bulk RNA-seq studies from GEO
#'
#' Downloads count matrices and sample metadata from the Gene Expression
#' Omnibus (GEO) for one or more accession IDs and returns a list of
#' \code{metaXpressStudy} objects ready for downstream analysis.
#'
#' @param accessions Character vector of GEO accession IDs
#'   (e.g., \code{c("GSE12345", "GSE67890")}).
#' @param count_type Character scalar. One of \code{"raw"} (default) or
#'   \code{"normalized"}. Always prefer \code{"raw"} for meta-analysis.
#' @param cache_dir Character scalar. Directory for caching downloaded files.
#'   Defaults to \code{tempdir()}.
#'
#' @return A named list of \code{\linkS4class{metaXpressStudy}} objects, one
#'   per accession. QC scores are filled by an internal call to
#'   \code{\link{mx_qc_study}}. Studies failing QC receive a warning but are
#'   returned (use \code{\link{mx_filter_studies}} to remove them).
#'
#' @references
#' Davis, S. & Meltzer, P.S. (2007) GEOquery: a bridge between the Gene
#' Expression Omnibus (GEO) and BioConductor. \emph{Bioinformatics},
#' \strong{23}(14), 1846--1847. \doi{10.1093/bioinformatics/btm254}
#'
#' Heberle, H. et al. (2025) A comprehensive framework for quality control
#' and meta-analysis of bulk RNA-seq data. \emph{Alzheimer's & Dementia},
#' \strong{21}(1), e70025. \doi{10.1002/alz.70025}
#'
#' @examples
#' \dontrun{
#'   studies <- mx_fetch_geo(c("GSE12345", "GSE67890"))
#'   studies <- mx_filter_studies(studies, qc_threshold = 7)
#' }
#'
#' @seealso \code{\link{mx_load_local}}, \code{\link{mx_qc_study}},
#'   \code{\link{mx_filter_studies}}
#'
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object
#'   controlling parallelization. Default: \code{BiocParallel::bpparam()}.
#'   Set to \code{BiocParallel::SerialParam()} for sequential execution.
#'
#' @importFrom methods new
#' @export
mx_fetch_geo <- function(accessions, count_type = "raw",
                          cache_dir = tempdir(),
                          BPPARAM = BiocParallel::bpparam()) {
  if (!is.character(accessions) || length(accessions) == 0)
    stop("'accessions' must be a non-empty character vector of GEO IDs")
  count_type <- match.arg(count_type, c("raw", "normalized"))

  message("Fetching ", length(accessions), " studies from GEO...")

  studies <- BiocParallel::bplapply(accessions, function(acc) {
    message("  Downloading: ", acc)
    .fetch_single_geo(acc, count_type = count_type, cache_dir = cache_dir)
  }, BPPARAM = BPPARAM)

  names(studies) <- accessions
  studies <- lapply(studies, mx_qc_study)
  studies
}

#' Fetch raw RNA-seq data from SRA
#'
#' Downloads count matrices for SRA project IDs. Requires \code{sratools}
#' (prefetch + fasterq-dump) to be installed and on the system PATH.
#'
#' @param srp_ids Character vector of SRA project IDs (e.g., \code{"SRP123456"}).
#' @param cache_dir Character scalar. Directory for caching downloaded files.
#'
#' @return A named list of \code{\linkS4class{metaXpressStudy}} objects.
#'
#' @examples
#' \dontrun{
#'   studies <- mx_fetch_sra("SRP123456")
#' }
#'
#' @seealso \code{\link{mx_fetch_geo}}, \code{\link{mx_load_local}}
#'
#' @export
mx_fetch_sra <- function(srp_ids, cache_dir = tempdir()) {
  if (!is.character(srp_ids) || length(srp_ids) == 0)
    stop("'srp_ids' must be a non-empty character vector")

  stop("mx_fetch_sra() is not yet implemented. ",
       "Use mx_load_local() with pre-downloaded count matrices.")
}

#' Load user-supplied count matrices as metaXpressStudy objects
#'
#' Constructs \code{metaXpressStudy} objects from local count matrix files
#' and matching metadata files. Supports CSV, TSV, and RDS formats.
#'
#' @param count_paths Character vector of file paths to count matrices. Each
#'   file must have genes as rows (with rownames) and samples as columns.
#' @param metadata_paths Character vector of file paths to sample metadata
#'   tables. Each must contain at minimum \code{condition} and \code{sample_id}
#'   columns. Length must match \code{count_paths}.
#' @param accessions Character vector of accession IDs (labels), one per study.
#'   Defaults to the basename of \code{count_paths} without extension.
#' @param organism Character scalar. Species for all studies.
#'   Default: \code{"Homo sapiens"}.
#'
#' @return A named list of \code{\linkS4class{metaXpressStudy}} objects.
#'
#' @examples
#' \dontrun{
#'   studies <- mx_load_local(
#'     count_paths    = c("study1_counts.csv", "study2_counts.csv"),
#'     metadata_paths = c("study1_meta.csv",   "study2_meta.csv")
#'   )
#' }
#'
#' @seealso \code{\link{mx_fetch_geo}}, \code{\link{mx_qc_study}}
#'
#' @importFrom methods new
#' @export
mx_load_local <- function(count_paths, metadata_paths,
                           accessions = NULL,
                           organism = "Homo sapiens") {
  if (!is.character(count_paths) || length(count_paths) == 0)
    stop("'count_paths' must be a non-empty character vector")
  if (length(count_paths) != length(metadata_paths))
    stop("'count_paths' and 'metadata_paths' must have the same length")

  if (is.null(accessions))
    accessions <- tools::file_path_sans_ext(basename(count_paths))

  studies <- mapply(function(cp, mp, acc) {
    counts   <- .read_matrix(cp)
    metadata <- .read_metadata(mp)

    # Ensure sample_id column exists
    if (!"sample_id" %in% colnames(metadata)) {
      if (!is.null(rownames(metadata)) &&
          !identical(rownames(metadata), as.character(seq_len(nrow(metadata)))))
        metadata$sample_id <- rownames(metadata)
      else
        metadata$sample_id <- colnames(counts)[seq_len(nrow(metadata))]
    }

    # Align metadata rows to count columns by sample_id
    if (all(metadata$sample_id %in% colnames(counts))) {
      metadata <- metadata[match(colnames(counts), metadata$sample_id), ]
      rownames(metadata) <- metadata$sample_id
    } else if (nrow(metadata) == ncol(counts)) {
      # Positional match if IDs don't align directly
      metadata$sample_id <- colnames(counts)
      rownames(metadata) <- colnames(counts)
    } else {
      stop("Cannot align metadata to count columns for study '", acc,
           "'. Check that sample IDs match.")
    }

    # Warn if counts appear pre-normalized (non-integer)
    if (!all(counts == floor(counts), na.rm = TRUE)) {
      warning("Study '", acc, "': count matrix contains non-integer values. ",
              "This likely indicates pre-normalized data (e.g., TPM, FPKM). ",
              "DESeq2/edgeR require raw integer counts for valid statistical ",
              "inference. Consider obtaining raw counts or excluding this study.")
    }

    new("metaXpressStudy",
        counts    = counts,
        metadata  = metadata,
        accession = acc,
        organism  = organism,
        qc_score  = NA_real_,
        de_result = data.frame())
  }, count_paths, metadata_paths, accessions, SIMPLIFY = FALSE)

  names(studies) <- accessions
  studies <- lapply(studies, mx_qc_study)
  studies
}

#' Score a study against the 10-point QC checklist
#'
#' Applies each of the 10 QC criteria defined by Heberle et al. (2025) and
#' stores the total score in the \code{qc_score} slot. Criteria that cannot
#' be evaluated (e.g., alignment rate not in metadata) are scored \code{NA}
#' and do not penalise the total.
#'
#' The 10 criteria are:
#' \enumerate{
#'   \item \strong{Minimum sample size:} at least 3 samples per condition
#'     group, required for reliable dispersion estimation in negative binomial
#'     models (Schurch et al. 2016).
#'   \item \strong{Sequencing depth:} median library size >= 10 million reads,
#'     ensuring adequate power for detecting low-abundance transcripts
#'     (Sims et al. 2014).
#'   \item \strong{Alignment rate:} mean alignment rate >= 70\%, indicating
#'     acceptable read quality (optional; skipped if not available).
#'   \item \strong{rRNA contamination:} mean rRNA fraction < 10\%, flagging
#'     insufficient ribosomal depletion (optional).
#'   \item \strong{Duplicate rate:} mean duplicate rate < 50\%, flagging
#'     library complexity issues (optional).
#'   \item \strong{Gene detection:} at least 15,000 genes detected in >= 50\%
#'     of samples, confirming adequate transcriptome coverage.
#'   \item \strong{Metadata completeness:} required columns \code{condition}
#'     and \code{sample_id} are present.
#'   \item \strong{Clear case/control:} exactly two condition levels, ensuring
#'     a well-defined contrast for differential expression.
#'   \item \strong{No batch-condition confounding:} batch and condition are not
#'     perfectly aliased, which would make batch correction impossible
#'     (Leek et al. 2010).
#'   \item \strong{Raw counts:} data are non-negative integers with maximum
#'     > 1,000, confirming untransformed count data suitable for count-based
#'     statistical models.
#' }
#'
#' @param study A \code{\linkS4class{metaXpressStudy}} object.
#'
#' @return The input \code{study} with the \code{qc_score} slot filled.
#'   A \code{qc_details} attribute is attached with a per-criterion breakdown.
#'
#' @references
#' Heberle, H. et al. (2025) A comprehensive framework for quality control
#' and meta-analysis of bulk RNA-seq data. \emph{Alzheimer's & Dementia},
#' \strong{21}(1), e70025. \doi{10.1002/alz.70025}
#'
#' Schurch, N.J. et al. (2016) How many biological replicates are needed in an
#' RNA-seq experiment and which differential expression tool should you use?
#' \emph{RNA}, \strong{22}(6), 839--851. \doi{10.1261/rna.053959.115}
#'
#' Leek, J.T. et al. (2010) Tackling the widespread and critical impact of
#' batch effects in high-throughput data. \emph{Nature Reviews Genetics},
#' \strong{11}(10), 733--739. \doi{10.1038/nrg2825}
#'
#' @examples
#' \dontrun{
#'   study <- mx_qc_study(study)
#'   attr(study, "qc_details")
#' }
#'
#' @export
mx_qc_study <- function(study) {
  if (!is(study, "metaXpressStudy"))
    stop("'study' must be a metaXpressStudy object")

  counts   <- study@counts
  metadata <- study@metadata
  score    <- integer(10)
  details  <- character(10)

  # 1. Minimum sample size
  min_n    <- min(table(metadata$condition))
  score[1] <- as.integer(min_n >= 3)
  details[1] <- sprintf("Min samples per group: %d (need >= 3)", min_n)

  # 2. Sequencing depth
  lib_sizes <- colSums(counts)
  med_depth <- median(lib_sizes)
  score[2]  <- as.integer(med_depth >= 10e6)
  details[2] <- sprintf("Median depth: %.1fM reads (need >= 10M)",
                         med_depth / 1e6)

  # 3. Alignment rate (optional)
  if ("alignment_rate" %in% colnames(metadata)) {
    mean_align <- mean(metadata$alignment_rate, na.rm = TRUE)
    score[3]   <- as.integer(mean_align >= 0.70)
    details[3] <- sprintf("Mean alignment: %.1f%% (need >= 70%%)",
                           mean_align * 100)
  } else {
    score[3]   <- NA_integer_
    details[3] <- "Alignment rate: not available in metadata"
  }

  # 4. rRNA contamination (optional)
  if ("rrna_fraction" %in% colnames(metadata)) {
    mean_rrna <- mean(metadata$rrna_fraction, na.rm = TRUE)
    score[4]  <- as.integer(mean_rrna < 0.10)
    details[4] <- sprintf("Mean rRNA: %.1f%% (need < 10%%)",
                           mean_rrna * 100)
  } else {
    score[4]  <- NA_integer_
    details[4] <- "rRNA fraction: not available in metadata"
  }

  # 5. Duplicate rate (optional)
  if ("duplicate_rate" %in% colnames(metadata)) {
    mean_dup <- mean(metadata$duplicate_rate, na.rm = TRUE)
    score[5] <- as.integer(mean_dup < 0.50)
    details[5] <- sprintf("Mean duplicate rate: %.1f%% (need < 50%%)",
                           mean_dup * 100)
  } else {
    score[5]  <- NA_integer_
    details[5] <- "Duplicate rate: not available in metadata"
  }

  # 6. Gene detection rate
  n_detected    <- rowSums(counts > 0)
  genes_in_half <- sum(n_detected >= (ncol(counts) * 0.5))
  score[6]      <- as.integer(genes_in_half >= 15000)
  details[6]    <- sprintf(
    "Genes detected in >=50%% samples: %d (need >= 15000)", genes_in_half)

  # 7. Metadata completeness
  required     <- c("condition", "sample_id")
  has_required <- all(required %in% colnames(metadata))
  score[7]     <- as.integer(has_required)
  details[7]   <- sprintf("Required metadata columns present: %s",
                           ifelse(has_required, "yes", "no"))

  # 8. Clear case/control definition
  n_levels <- length(unique(metadata$condition))
  score[8] <- as.integer(n_levels == 2)
  details[8] <- sprintf("Condition levels: %d (need exactly 2)", n_levels)

  # 9. No batch-condition confounding
  if ("batch" %in% colnames(metadata)) {
    tab        <- table(metadata$condition, metadata$batch)
    confounded <- all(apply(tab, 2, function(x) sum(x > 0) == 1))
    score[9]   <- as.integer(!confounded)
    details[9] <- sprintf("Batch-condition confounding: %s",
                           ifelse(confounded, "YES (problem)", "no"))
  } else {
    score[9]   <- 1L
    details[9] <- "No batch variable present"
  }

  # 10. Raw counts available
  is_int  <- all(counts == floor(counts), na.rm = TRUE)
  is_pos  <- all(counts >= 0, na.rm = TRUE)
  is_raw  <- is_int && is_pos && max(counts, na.rm = TRUE) > 1000
  score[10]   <- as.integer(is_raw)
  details[10] <- sprintf("Raw integer counts: %s",
                          ifelse(is_raw, "yes", "no (likely pre-normalized)"))

  # Score normalized to 0-10 scale, counting only evaluable criteria
  n_evaluable    <- sum(!is.na(score))
  raw_score      <- sum(score, na.rm = TRUE)
  total_score    <- round(raw_score / n_evaluable * 10)
  study@qc_score <- total_score

  attr(study, "qc_details") <- data.frame(
    criterion = seq_len(10),
    score     = score,
    details   = details,
    stringsAsFactors = FALSE
  )

  if (total_score < 7)
    warning(sprintf(
      "Study %s has low QC score (%d/10). Consider excluding with mx_filter_studies().",
      study@accession, total_score))

  study
}

#' Automatically cluster samples by metadata
#'
#' Uses metadata fields to assign samples to case/control groups, mimicking
#' the sampleclusteR approach for automated group assignment from GEO
#' metadata.
#'
#' @param study A \code{\linkS4class{metaXpressStudy}} object.
#'
#' @return The input \code{study} with the \code{condition} column updated
#'   in \code{metadata}.
#'
#' @references
#' Coke, T., Niranjan, M. & Ewing, R.M. (2025) sampleclusteR: automated
#' case/control group assignment from GEO metadata. \emph{bioRxiv}.
#' \doi{10.1101/2025.04.10.648129}
#'
#' @examples
#' \dontrun{
#'   study <- mx_cluster_samples(study)
#' }
#'
#' @export
mx_cluster_samples <- function(study) {
  if (!is(study, "metaXpressStudy"))
    stop("'study' must be a metaXpressStudy object")
  stop("mx_cluster_samples() is not yet implemented.")
}

#' Filter studies below a QC threshold
#'
#' Removes studies whose \code{qc_score} is below the specified threshold.
#' Studies with \code{NA} QC scores (i.e., \code{\link{mx_qc_study}} has not
#' been run) are always removed with a warning.
#'
#' @param studies A named list of \code{\linkS4class{metaXpressStudy}} objects.
#' @param qc_threshold Numeric scalar. Minimum acceptable QC score (inclusive).
#'   Default: \code{7}.
#'
#' @return A filtered list of \code{\linkS4class{metaXpressStudy}} objects.
#'
#' @examples
#' \dontrun{
#'   studies <- mx_fetch_geo(c("GSE12345", "GSE67890"))
#'   studies <- mx_filter_studies(studies, qc_threshold = 7)
#' }
#'
#' @seealso \code{\link{mx_qc_study}}
#'
#' @export
mx_filter_studies <- function(studies, qc_threshold = 7) {
  if (!is.list(studies))
    stop("'studies' must be a list of metaXpressStudy objects")
  if (!is.numeric(qc_threshold) || length(qc_threshold) != 1)
    stop("'qc_threshold' must be a single numeric value")

  scores <- vapply(studies, function(s) s@qc_score, numeric(1))

  na_idx <- is.na(scores)
  if (any(na_idx))
    warning(sum(na_idx), " study/studies have NA QC scores (mx_qc_study ",
            "not run?) and will be removed: ",
            paste(names(studies)[na_idx], collapse = ", "))

  keep <- !na_idx & scores >= qc_threshold
  removed <- sum(!keep)

  if (removed > 0)
    message(removed, " study/studies removed (QC score < ", qc_threshold, "): ",
            paste(names(studies)[!keep], collapse = ", "))

  studies[keep]
}

# ============================================================================
# Internal helpers
# ============================================================================

.fetch_single_geo <- function(accession, count_type, cache_dir) {
  gse <- GEOquery::getGEO(accession, destdir = cache_dir, GSEMatrix = TRUE)
  if (is.list(gse)) gse <- gse[[1]]

  supp_files <- tryCatch(
    GEOquery::getGEOSuppFiles(accession, baseDir = cache_dir),
    error = function(e) NULL
  )

  if (is.null(supp_files)) {
    warning("No supplementary count files found for ", accession,
            ". Returning empty study.")
    return(new("metaXpressStudy",
               accession = accession,
               organism  = as.character(
                 Biobase::annotation(gse)[[1]]),
               qc_score  = NA_real_,
               de_result = data.frame()))
  }

  counts   <- .parse_geo_counts(supp_files, accession)
  metadata <- .parse_geo_metadata(gse)
  metadata <- .align_metadata_to_counts(metadata, counts)

  new("metaXpressStudy",
      counts    = counts,
      metadata  = metadata,
      accession = accession,
      organism  = "Homo sapiens",
      qc_score  = NA_real_,
      de_result = data.frame())
}

.parse_geo_counts <- function(supp_files, accession) {
  count_files <- rownames(supp_files)[grepl(
    "count|counts|htseq|featurecounts|rsem",
    rownames(supp_files), ignore.case = TRUE)]

  if (length(count_files) == 0)
    stop("Could not identify count file for ", accession,
         ". Please use mx_load_local() instead.")

  .read_matrix(count_files[1])
}

.parse_geo_metadata <- function(gse) {
  pd <- Biobase::pData(gse)
  meta <- as.data.frame(pd, stringsAsFactors = FALSE)
  if (!"sample_id" %in% colnames(meta))
    meta$sample_id <- rownames(meta)
  if (!"condition" %in% colnames(meta))
    meta$condition <- NA_character_
  meta
}

.align_metadata_to_counts <- function(metadata, counts) {
  common <- intersect(colnames(counts), metadata$sample_id)
  if (length(common) == 0)
    common <- intersect(colnames(counts), rownames(metadata))
  if (length(common) == 0)
    stop("Cannot match count column names to metadata sample IDs")
  metadata[match(colnames(counts), metadata$sample_id), , drop = FALSE]
}

.read_matrix <- function(path) {
  ext <- tools::file_ext(path)
  mat <- switch(ext,
    rds = readRDS(path),
    csv = utils::read.csv(path, row.names = 1, check.names = FALSE),
    tsv = ,
    txt = utils::read.delim(path, row.names = 1, check.names = FALSE),
    gz  = {
      con <- gzfile(path)
      on.exit(close(con))
      utils::read.delim(con, row.names = 1, check.names = FALSE)
    },
    stop("Unsupported file format: .", ext)
  )
  as.matrix(mat)
}

.read_metadata <- function(path) {
  ext <- tools::file_ext(path)
  meta <- switch(ext,
    rds = readRDS(path),
    csv = utils::read.csv(path, stringsAsFactors = FALSE),
    tsv = ,
    txt = utils::read.delim(path, stringsAsFactors = FALSE),
    stop("Unsupported metadata format: .", ext)
  )

  # Normalize column names: accept common variants
  cn <- colnames(meta)
  cn_lower <- tolower(cn)

  # Map "Sample" / "sample" / "SampleID" -> "sample_id"
  sample_col <- which(cn_lower %in% c("sample", "sampleid", "sample.id",
                                       "sample_name", "samplename", "geo_accession"))
  if (length(sample_col) > 0 && !("sample_id" %in% cn_lower))
    colnames(meta)[sample_col[1]] <- "sample_id"

  # Map "Condition" / "Group" / "Status" -> "condition"
  cond_col <- which(cn_lower %in% c("condition", "group", "status",
                                     "disease", "phenotype", "class"))
  if (length(cond_col) > 0 && !("condition" %in% cn))
    colnames(meta)[cond_col[1]] <- "condition"

  # Standardize condition labels to lowercase for consistency
  if ("condition" %in% colnames(meta)) {
    meta$condition <- tolower(trimws(meta$condition))
  }

  meta
}
