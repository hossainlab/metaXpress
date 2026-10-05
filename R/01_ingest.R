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

  .require_package("GEOquery", "fetching data from GEO")
  .require_package("Biobase", "parsing GEO metadata")

  message("Fetching ", length(accessions), " studies from GEO...")

  # Use serial processing for GEO downloads to avoid connection issues
  studies <- lapply(accessions, function(acc) {
    message("  Downloading: ", acc)
    tryCatch(
      .fetch_single_geo(acc, count_type = count_type, cache_dir = cache_dir),
      error = function(e) {
        warning("Failed to fetch ", acc, ": ", conditionMessage(e))
        NULL
      }
    )
  })

  names(studies) <- accessions

  # Remove NULLs (failed downloads)
  failed <- vapply(studies, is.null, logical(1))
  if (any(failed)) {
    warning(sum(failed), " study/studies could not be downloaded: ",
            paste(accessions[failed], collapse = ", "))
    studies <- studies[!failed]
  }

  if (length(studies) == 0)
    stop("All downloads failed. Check accession IDs and network connection.")

  # Run QC on studies that have non-empty counts
  studies <- lapply(studies, function(s) {
    if (nrow(s@counts) > 0 && ncol(s@counts) > 0 &&
        nrow(s@metadata) > 0 && "condition" %in% colnames(s@metadata) &&
        !all(is.na(s@metadata$condition))) {
      mx_qc_study(s)
    } else {
      warning("Study ", s@accession, " has empty counts or missing condition ",
              "labels. QC skipped. Assign conditions manually: ",
              "study@metadata$condition <- c(...)")
      s
    }
  })

  studies
}

#' Fetch raw RNA-seq data from SRA
#'
#' Downloads pre-quantified count matrices for SRA project IDs using the
#' \\pkg{recount3} package. This avoids computationally expensive local
#' alignment while ensuring uniform processing across studies (Collado-Torres
#' et al. 2017; Wilks et al. 2021).
#'
#' @param srp_ids Character vector of SRA project IDs (e.g., \code{"SRP123456"}).
#' @param cache_dir Character scalar. Directory for caching downloaded files.
#' @param organism Character scalar. Species for the studies. Default:
#'   \code{"human"}. (Note: \pkg{recount3} uses \code{"human"} or \code{"mouse"}).
#'
#' @return A named list of \code{\linkS4class{metaXpressStudy}} objects.
#'
#' @references
#' Wilks, C. et al. (2021) recount3: summaries and queries for large-scale RNA-seq
#' expression and splicing. \emph{Genome Biology}, \strong{22}(1), 323.
#' \doi{10.1186/s13059-021-02533-6}
#'
#' @examples
#' \dontrun{
#'   studies <- mx_fetch_sra("SRP123456")
#' }
#'
#' @seealso \code{\link{mx_fetch_geo}}, \code{\link{mx_load_local}}
#'
#' @export
mx_fetch_sra <- function(srp_ids, cache_dir = tempdir(), organism = "human") {
  if (!is.character(srp_ids) || length(srp_ids) == 0)
    stop("'srp_ids' must be a non-empty character vector")
  
  .require_package("recount3", "fetching data from SRA")
  
  message("Fetching ", length(srp_ids), " projects from recount3...")
  
  studies <- lapply(srp_ids, function(srp) {
    message("  Downloading: ", srp)
    tryCatch({
      # Find project in recount3
      proj_info <- recount3::available_projects(organism = organism)
      proj_idx <- which(proj_info$project == srp)
      
      if (length(proj_idx) == 0) {
        warning("Project ", srp, " not found in recount3 database.")
        return(NULL)
      }
      
      # Download RangedSummarizedExperiment
      rse <- recount3::create_rse(proj_info[proj_idx[1], ])
      
      # Compute read counts from coverage (recount3 standard transformation)
      assay(rse, "counts") <- recount3::transform_counts(rse)
      counts <- as.matrix(assay(rse, "counts"))
      
      # Extract metadata
      meta <- as.data.frame(colData(rse))
      
      # Attempt to find sample_id and condition
      meta$sample_id <- if ("external_id" %in% colnames(meta)) meta$external_id else rownames(meta)
      
      # Use the phenotype/attributes columns to auto-detect condition
      condition_col <- .detect_condition_column(meta)
      if (!is.null(condition_col)) {
        meta$condition <- tolower(trimws(as.character(meta[[condition_col]])))
      } else {
        meta$condition <- NA_character_
        warning("Could not auto-detect condition for SRA project ", srp, ". Assign manually.")
      }
      
      study <- new("metaXpressStudy",
                   counts    = round(counts),
                   metadata  = meta,
                   accession = srp,
                   organism  = organism,
                   qc_score  = NA_real_,
                   de_result = data.frame())
      
      # Run QC if condition exists
      if (!all(is.na(meta$condition))) {
        study <- mx_qc_study(study)
      }
      study
    }, error = function(e) {
      warning("Failed to fetch ", srp, " from recount3: ", conditionMessage(e))
      NULL
    })
  })
  
  names(studies) <- srp_ids
  valid <- !vapply(studies, is.null, logical(1))
  studies[valid]
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
#' Uses text mining on metadata fields to assign samples to case/control groups,
#' mimicking the \\pkg{sampleclusteR} approach for automated group assignment from
#' GEO metadata. It computes string distances between sample descriptions and
#' uses hierarchical clustering to partition them into two groups.
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
    
  .require_package("stringdist", "string distance calculation for clustering")
  
  meta <- study@metadata
  
  # Concatenate all character/factor columns for each sample into a single text document
  text_cols <- vapply(meta, function(x) is.character(x) || is.factor(x), logical(1))
  if (sum(text_cols) == 0) {
    warning("No character columns found in metadata to cluster on.")
    return(study)
  }
  
  docs <- apply(meta[, text_cols, drop=FALSE], 1, function(x) {
    paste(na.omit(as.character(x)), collapse = " ")
  })
  
  # Compute Jaro-Winkler string distance matrix
  dist_mat <- stringdist::stringdistmatrix(docs, docs, method = "jw")
  
  # Hierarchical clustering
  hc <- hclust(as.dist(dist_mat), method = "ward.D2")
  groups <- cutree(hc, k = 2)
  
  # Assign the groups as "group1" and "group2"
  study@metadata$condition <- paste0("group", groups)
  
  message("Auto-clustered ", nrow(meta), " samples into 2 groups based on metadata text.")
  study
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
  .require_package("GEOquery", "fetching data from GEO")
  .require_package("Biobase", "parsing GEO metadata")

  # --- Download GSE metadata ---
  gse <- tryCatch(
    GEOquery::getGEO(accession, destdir = cache_dir, GSEMatrix = TRUE),
    error = function(e) {
      stop("Failed to download metadata for ", accession, ": ",
           conditionMessage(e))
    }
  )
  if (is.list(gse)) gse <- gse[[1]]

  # --- Detect organism from GEO metadata ---
  organism <- .detect_organism(gse)

  # --- Download supplementary (count) files ---
  supp_files <- tryCatch(
    GEOquery::getGEOSuppFiles(accession, baseDir = cache_dir),
    error = function(e) NULL
  )

  # --- Extract count matrix ---
  counts <- NULL

  if (!is.null(supp_files) && nrow(supp_files) > 0) {
    # Try to find a count matrix in supplementary files
    counts <- tryCatch(
      .parse_geo_counts(supp_files, accession, cache_dir),
      error = function(e) {
        warning("Could not parse supplementary files for ", accession,
                ": ", conditionMessage(e))
        NULL
      }
    )
  }

  # Fallback: use the expression matrix from the GSEMatrix
  if (is.null(counts) && count_type == "normalized") {
    message("  Using expression matrix from GSEMatrix for ", accession)
    expr_mat <- Biobase::exprs(gse)
    if (nrow(expr_mat) > 0 && ncol(expr_mat) > 0) {
      counts <- as.matrix(expr_mat)
    }
  }

  if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) {
    warning("No count matrix could be extracted for ", accession,
            ". Use mx_load_local() with manually downloaded count files.")
    # Return a minimal valid object (empty counts, but valid metadata)
    empty_meta <- data.frame(sample_id = character(0),
                             condition = character(0),
                             stringsAsFactors = FALSE)
    return(new("metaXpressStudy",
               counts    = matrix(0L, nrow = 0, ncol = 0),
               metadata  = empty_meta,
               accession = accession,
               organism  = organism,
               qc_score  = NA_real_,
               de_result = data.frame()))
  }

  # --- Parse sample metadata from GEO phenotype data ---
  metadata <- .parse_geo_metadata(gse)

  # --- Align metadata to count columns ---
  metadata <- .align_metadata_to_counts(metadata, counts)

  # --- Warn if counts look pre-normalized ---
  if (!all(counts == floor(counts), na.rm = TRUE)) {
    warning("Study '", accession, "': count matrix contains non-integer values. ",
            "This likely indicates pre-normalized data. ",
            "DESeq2/edgeR require raw integer counts for valid inference.")
  }

  new("metaXpressStudy",
      counts    = counts,
      metadata  = metadata,
      accession = accession,
      organism  = organism,
      qc_score  = NA_real_,
      de_result = data.frame())
}

.parse_geo_counts <- function(supp_files, accession, cache_dir) {
  fnames <- rownames(supp_files)

  # Unpack tar.gz archives if present
  tar_files <- fnames[grepl("\\.tar(\\.gz)?$", fnames, ignore.case = TRUE)]
  if (length(tar_files) > 0) {
    extract_dir <- file.path(cache_dir, accession, "extracted")
    dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
    for (tf in tar_files) {
      utils::untar(tf, exdir = extract_dir)
    }
    extracted <- list.files(extract_dir, full.names = TRUE, recursive = TRUE)
    fnames <- c(fnames, extracted)
  }

  # Look for count matrix files (prioritize raw counts)
  count_patterns <- c(
    "raw[_\\.]?count", "read[_\\.]?count", "gene[_\\.]?count",
    "counts?\\.(?:csv|tsv|txt)",
    "htseq", "featurecounts", "rsem[_\\.]?gene",
    "expected_count", "raw_count"
  )
  pattern <- paste(count_patterns, collapse = "|")
  count_files <- fnames[grepl(pattern, basename(fnames), ignore.case = TRUE)]

  # Exclude obvious non-count files
  count_files <- count_files[!grepl(
    "normalized|fpkm|tpm|rpkm|metadata|sample|annotation|readme|design",
    basename(count_files), ignore.case = TRUE)]

  # If no specific count files, try any tabular supplementary file
  if (length(count_files) == 0) {
    count_files <- fnames[grepl(
      "\\.(csv|tsv|txt|gz)$", fnames, ignore.case = TRUE)]
    count_files <- count_files[!grepl("\\.tar", count_files, ignore.case = TRUE)]
  }

  if (length(count_files) == 0)
    stop("No count matrix files found in supplementary data for ", accession)

  # Determine if these are per-sample files or a single matrix
  # Heuristic: if filenames contain GSM IDs and there are multiple files,
  # they are per-sample files that need merging
  is_per_sample <- length(count_files) > 2 &&
    sum(grepl("GSM[0-9]+", basename(count_files))) > 2

  if (is_per_sample) {
    message("  Detected ", length(count_files),
            " per-sample count files for ", accession, "; merging...")
    mat <- .merge_per_sample_counts(count_files)
  } else {
    # Single matrix file — read directly
    mat <- .read_matrix(count_files[1])
  }

  mat
}

#' Merge per-sample count files into a single matrix
#'
#' Handles various per-sample formats: 2-column (gene + count),
#' HTSeq output, featureCounts output, and multi-column formats
#' with gene symbol and count columns.
#' @noRd
.merge_per_sample_counts <- function(file_paths) {
  # Extract sample names from filenames (remove GSMxxxxx_ prefix if present)
  sample_names <- basename(file_paths)
  sample_names <- sub("\\.(txt|csv|tsv)(\\.gz)?$", "", sample_names,
                       ignore.case = TRUE)
  # Use GSM ID as sample name
  gsm_match <- regmatches(sample_names, regexpr("GSM[0-9]+", sample_names))
  if (length(gsm_match) == length(sample_names) &&
      all(nchar(gsm_match) > 0)) {
    sample_names <- gsm_match
  }

  counts_list <- lapply(seq_along(file_paths), function(i) {
    fp <- file_paths[i]
    tryCatch({
      parsed <- .read_per_sample_count(fp)
      parsed
    }, error = function(e) {
      warning("Could not parse per-sample file: ", basename(fp),
              " (", conditionMessage(e), ")")
      NULL
    })
  })

  # Remove failures
  valid <- !vapply(counts_list, is.null, logical(1))
  counts_list  <- counts_list[valid]
  sample_names <- sample_names[valid]

  if (length(counts_list) == 0)
    stop("Could not read any per-sample count files")

  # Union all genes, fill missing with 0
  all_genes <- Reduce(union, lapply(counts_list, names))
  merged <- matrix(0L, nrow = length(all_genes), ncol = length(counts_list),
                   dimnames = list(all_genes, sample_names))

  for (i in seq_along(counts_list)) {
    v <- counts_list[[i]]
    merged[names(v), i] <- as.integer(v)
  }

  merged
}

#' Read a single per-sample count file and return a named integer vector
#'
#' Handles: 2-column (gene, count), HTSeq (gene, count with summary lines),
#' featureCounts (multi-column with count as last column), and custom formats
#' with gene symbol + count columns.
#' @return A named integer vector: names = gene IDs, values = counts.
#' @noRd
.read_per_sample_count <- function(path) {
  ext <- tolower(tools::file_ext(path))

  # Read the file
  if (ext == "gz") {
    con <- gzfile(path, "rt")
    on.exit(close(con))
    # Peek at first line to check for header
    first_line <- readLines(con, n = 1)
    has_header <- !grepl("^(ENSG|ENSMUSG|[0-9]+\\t)", first_line)
    seek(con, 0)
    d <- utils::read.delim(con, header = has_header, stringsAsFactors = FALSE,
                           check.names = FALSE)
  } else {
    sep <- if (ext == "csv") "," else "\t"
    first_line <- readLines(path, n = 1)
    has_header <- !grepl("^(ENSG|ENSMUSG|[0-9]+[,\\t])", first_line)
    d <- utils::read.delim(path, sep = sep, header = has_header,
                           stringsAsFactors = FALSE, check.names = FALSE)
  }

  if (nrow(d) == 0 || ncol(d) < 2) return(NULL)

  # Remove HTSeq summary lines (__no_feature, __ambiguous, etc.)
  if (is.character(d[[1]])) {
    d <- d[!grepl("^__", d[[1]]), , drop = FALSE]
  }

  cn <- colnames(d)
  cn_lower <- tolower(cn)

  # --- Detect gene ID column ---
  gene_col <- NULL

  # Look for explicit gene ID column names
  gene_names <- c("gene_id", "geneid", "gene", "symbol", "gene_symbol",
                  "gene_name", "final_id", "ensembl_gene_id")
  for (gn in gene_names) {
    idx <- which(cn_lower == gn)
    if (length(idx) > 0) { gene_col <- idx[1]; break }
  }

  # Fallback: first character column

  if (is.null(gene_col)) {
    char_cols <- which(vapply(d, is.character, logical(1)))
    if (length(char_cols) > 0) gene_col <- char_cols[1]
    else gene_col <- 1
  }

  # --- Detect count column ---
  count_col <- NULL

  # Look for explicit count column names
  count_names <- c("count", "counts", "expected_count", "numreads",
                   "rpb", "reads", "raw_count")
  for (cn_name in count_names) {
    idx <- which(cn_lower == cn_name)
    if (length(idx) > 0) { count_col <- idx[1]; break }
  }

  # Fallback: for 2-column files, use column 2
  if (is.null(count_col) && ncol(d) == 2) {
    count_col <- 2
  }

  # Fallback: last numeric column (featureCounts convention)
  if (is.null(count_col)) {
    num_cols <- which(vapply(d, is.numeric, logical(1)))
    if (length(num_cols) > 0) count_col <- num_cols[length(num_cols)]
  }

  if (is.null(count_col))
    stop("Could not identify count column")

  gene_ids <- as.character(d[[gene_col]])
  counts   <- as.numeric(d[[count_col]])

  # If gene IDs look like coordinates (chr1:xxx), use Symbol column instead
  if (all(grepl("^chr", head(gene_ids, 5)))) {
    sym_col <- which(cn_lower %in% c("symbol", "gene_symbol", "gene_name",
                                      "gene"))
    if (length(sym_col) > 0) {
      gene_ids <- as.character(d[[sym_col[1]]])
    }
  }

  # Aggregate by gene if multiple entries per gene
  if (anyDuplicated(gene_ids)) {
    agg <- stats::aggregate(counts, by = list(gene = gene_ids), FUN = sum)
    result <- agg$x
    names(result) <- agg$gene
  } else {
    result <- counts
    names(result) <- gene_ids
  }

  # Remove NA/empty gene names
  result <- result[!is.na(names(result)) & names(result) != ""]
  result
}

.parse_geo_metadata <- function(gse) {
  pd <- Biobase::pData(gse)
  meta <- as.data.frame(pd, stringsAsFactors = FALSE)

  # Set sample_id from GEO accession column or rownames
  if ("geo_accession" %in% colnames(meta)) {
    meta$sample_id <- meta$geo_accession
  } else {
    meta$sample_id <- rownames(meta)
  }

  # --- Detect condition column from GEO phenotype data ---
  # GEO stores sample characteristics as "characteristic_name:ch1" columns
  if (!"condition" %in% colnames(meta)) {
    condition_col <- .detect_condition_column(meta)
    if (!is.null(condition_col)) {
      raw_vals <- meta[[condition_col]]
      # GEO characteristic values often have "key: value" format
      cleaned <- sub("^[^:]+:\\s*", "", raw_vals)
      meta$condition <- tolower(trimws(cleaned))
      message("  Auto-detected condition column: '", condition_col,
              "' -> levels: ", paste(unique(meta$condition), collapse = ", "))
    } else {
      meta$condition <- NA_character_
      warning("Could not auto-detect condition column for study. ",
              "Assign manually: study@metadata$condition <- ...")
    }
  }

  meta
}

#' Detect which GEO phenotype column contains the case/control condition
#' @noRd
.detect_condition_column <- function(meta) {
  cn <- colnames(meta)

  # GEO system columns to always skip
  skip_patterns <- c("geo_accession", "sample_id", "title", "description",
                     "source_name_ch1", "platform_id", "contact",
                     "data_processing", "supplementary_file",
                     "data_row_count", "channel_count", "scan_date",
                     "hyb_protocol", "label_ch", "label_protocol",
                     "extract_protocol", "submission_date",
                     "last_update_date", "organism_ch1",
                     "molecule_ch1", "taxid_ch1", "relation",
                     "^status$")

  .is_skip <- function(col) {
    any(vapply(skip_patterns, function(p) grepl(p, col, ignore.case = TRUE),
               logical(1)))
  }

  # Priority 1: GEO characteristic columns (":ch1" suffix) — most reliable
  ch_cols <- cn[grepl(":ch1$", cn)]

  # Look for characteristic columns that suggest case/control
  condition_keywords <- c("disease", "state", "condition", "group",
                          "phenotype", "treatment", "diagnosis",
                          "histology", "sample.type", "sample_type",
                          "tissue.type", "tissue_type")
  for (kw in condition_keywords) {
    match_idx <- which(grepl(kw, ch_cols, ignore.case = TRUE))
    if (length(match_idx) > 0) {
      col <- ch_cols[match_idx[1]]
      vals <- unique(sub("^[^:]+:\\s*", "", meta[[col]]))
      vals <- vals[!is.na(vals) & vals != ""]
      if (length(vals) == 2) return(col)
    }
  }

  # Priority 2: any :ch1 column with exactly 2 unique values
  for (col in ch_cols) {
    vals <- unique(sub("^[^:]+:\\s*", "", meta[[col]]))
    vals <- vals[!is.na(vals) & vals != ""]
    if (length(vals) == 2) return(col)
  }

  # Priority 3: explicit non-:ch1 column names (user-provided metadata)
  explicit_names <- c("condition", "group", "disease_state", "disease",
                      "phenotype", "treatment", "diagnosis",
                      "sample_type")
  for (nm in explicit_names) {
    match_idx <- which(tolower(cn) == nm)
    if (length(match_idx) > 0) {
      col <- cn[match_idx[1]]
      vals <- unique(meta[[col]])
      vals <- vals[!is.na(vals) & vals != ""]
      if (length(vals) == 2) return(col)
    }
  }

  # Priority 4: any non-system column with exactly 2 unique values
  for (col in cn) {
    if (.is_skip(col)) next
    if (col %in% ch_cols) next  # already checked
    vals <- unique(meta[[col]])
    vals <- vals[!is.na(vals) & vals != ""]
    if (length(vals) == 2) return(col)
  }

  NULL
}

#' Detect organism from GEO ExpressionSet
#' @noRd
.detect_organism <- function(gse) {
  # Try the organism_ch1 column in phenotype data
  pd <- Biobase::pData(gse)
  if ("organism_ch1" %in% colnames(pd)) {
    org <- unique(pd$organism_ch1)
    org <- org[!is.na(org) & org != ""]
    if (length(org) == 1) return(org)
  }

  # Fallback: check experimentData
  ed <- tryCatch(Biobase::experimentData(gse), error = function(e) NULL)
  if (!is.null(ed)) {
    other <- tryCatch(ed@other, error = function(e) list())
    if ("organism" %in% names(other))
      return(other[["organism"]])
  }

  "Homo sapiens"
}

.align_metadata_to_counts <- function(metadata, counts) {
  count_samples <- colnames(counts)

  # Try matching by sample_id
  common <- intersect(count_samples, metadata$sample_id)

  # Fallback: try matching by rownames
  if (length(common) == 0)
    common <- intersect(count_samples, rownames(metadata))

  # Fallback: try matching GEO accessions (GSMxxxxx) in column names
  if (length(common) == 0 && "geo_accession" %in% colnames(metadata)) {
    # Count columns might contain GSM IDs as part of longer names
    gsm_ids <- metadata$geo_accession
    matches <- vapply(count_samples, function(cs) {
      hit <- which(vapply(gsm_ids, function(g) grepl(g, cs, fixed = TRUE),
                          logical(1)))
      if (length(hit) == 1) hit else NA_integer_
    }, integer(1))
    if (sum(!is.na(matches)) >= 2) {
      metadata <- metadata[matches[!is.na(matches)], , drop = FALSE]
      counts_keep <- count_samples[!is.na(matches)]
      metadata$sample_id <- counts_keep
      rownames(metadata) <- counts_keep
      return(metadata)
    }
  }

  if (length(common) == 0) {
    # Last resort: positional matching if dimensions agree
    if (nrow(metadata) == ncol(counts)) {
      warning("Matching metadata to count columns by position (no ID overlap). ",
              "Verify sample order manually.")
      metadata$sample_id <- count_samples
      rownames(metadata) <- count_samples
      return(metadata)
    }
    stop("Cannot match count column names to metadata sample IDs. ",
         "Count columns: ", paste(head(count_samples, 5), collapse = ", "),
         " ... Metadata IDs: ", paste(head(metadata$sample_id, 5), collapse = ", "))
  }

  # Subset counts to matched samples and reorder metadata
  idx <- match(count_samples, metadata$sample_id)
  if (any(is.na(idx))) {
    # Some count columns don't have metadata — keep only matched
    has_meta <- !is.na(idx)
    warning(sum(!has_meta), " sample(s) in count matrix have no metadata match. ",
            "Keeping ", sum(has_meta), " matched samples.")
    idx <- idx[has_meta]
  }
  metadata[idx, , drop = FALSE]
}

.read_matrix <- function(path) {
  ext <- tolower(tools::file_ext(path))

  # For gzipped files, detect the inner format and read directly
  if (ext == "gz") {
    inner_ext <- tolower(tools::file_ext(sub("\\.gz$", "", path,
                                             ignore.case = TRUE)))
    mat <- if (inner_ext == "csv") {
      utils::read.csv(gzfile(path, "rt"), row.names = 1, check.names = FALSE)
    } else {
      utils::read.delim(gzfile(path, "rt"), row.names = 1, check.names = FALSE)
    }
    return(as.matrix(mat))
  }

  mat <- switch(ext,
    rds = readRDS(path),
    csv = utils::read.csv(path, row.names = 1, check.names = FALSE),
    tsv = ,
    txt = utils::read.delim(path, row.names = 1, check.names = FALSE),
    stop("Unsupported file format: .", ext,
         ". Supported: .csv, .tsv, .txt, .rds, .gz")
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
