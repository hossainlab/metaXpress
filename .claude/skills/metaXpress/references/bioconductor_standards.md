# Bioconductor Coding Standards for metaXpress

## Documentation (roxygen2)

Every exported function MUST have complete roxygen2 documentation.

```r
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
#' @return A named list of \code{\linkS4class{metaXpressStudy}} objects,
#'   one per accession. Studies failing QC are returned with a warning
#'   but not removed (use \code{\link{mx_filter_studies}} for filtering).
#'
#' @references
#' Heberle et al. (2025) Alzheimer's & Dementia.
#' \doi{10.1002/alz.70025}
#'
#' @examples
#' # Fetch two GEO studies
#' \dontrun{
#'   studies <- mx_fetch_geo(c("GSE12345", "GSE67890"))
#'   studies <- mx_filter_studies(studies, qc_threshold = 7)
#' }
#'
#' @seealso \code{\link{mx_load_local}}, \code{\link{mx_qc_study}},
#'   \code{\link{mx_filter_studies}}
#'
#' @export
mx_fetch_geo <- function(accessions, count_type = "raw",
                          cache_dir = tempdir()) {
  # implementation
}
```

### Required roxygen2 tags for all exported functions:
- `@param` — every parameter, with type and valid values
- `@return` — what is returned and its class/structure
- `@examples` — at least one runnable example (use `\dontrun{}` for network calls)
- `@export` — marks function as exported
- `@seealso` — link to related functions

### Additional tags when applicable:
- `@references` — cite the paper the method is from (with `\doi{}`)
- `@importFrom` — preferred over `@import` for specific function imports

---

## NAMESPACE Management

Use `@importFrom` in roxygen2 — never `library()` or `require()` inside functions.

```r
# Good — specific imports
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom SummarizedExperiment assay colData
#' @importFrom methods new is validObject

# Bad — never do this inside a function body
library(DESeq2)
require(ggplot2)
```

For functions from other packages used inline, use `::`:
```r
# Good
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, ...)

# Also acceptable if @importFrom is declared
dds <- DESeqDataSetFromMatrix(countData = counts, ...)
```

---

## Input Validation

Every exported function must validate inputs at the top using `stopifnot()` or informative errors.

```r
mx_meta <- function(de_results, method = "random_effects",
                    min_studies = 2, alpha = 0.05) {

  # --- Input validation ---
  if (!is.list(de_results))
    stop("'de_results' must be a list of data.frames from mx_de_all()")

  valid_methods <- c("fisher", "stouffer", "inverse_normal",
                     "fixed_effects", "random_effects", "awmeta")
  method <- match.arg(method, valid_methods)

  if (!is.numeric(min_studies) || min_studies < 2)
    stop("'min_studies' must be an integer >= 2")

  if (method == "random_effects" && length(de_results) < 3)
    warning("Random effects model is unreliable with fewer than 3 studies. ",
            "Consider method = 'fisher' or method = 'stouffer'.")

  # Check required columns in each DE result
  required_cols <- c("gene_id", "log2FC", "pvalue", "padj")
  for (i in seq_along(de_results)) {
    missing_cols <- setdiff(required_cols, colnames(de_results[[i]]))
    if (length(missing_cols) > 0)
      stop(sprintf("de_results[[%d]] is missing columns: %s",
                   i, paste(missing_cols, collapse = ", ")))
  }

  # --- Implementation ---
}
```

---

## Error and Warning Conventions

```r
# Errors — unrecoverable, stop execution
stop("'counts' must be a matrix with rownames as gene IDs")

# Warnings — recoverable, execution continues
warning("Study GSE12345 has QC score < 5; results may be unreliable")

# Messages — informational progress updates
message("Fetching ", length(accessions), " studies from GEO...")
message("  Running DESeq2 on study: ", study@accession)

# Never use cat() or print() for user-facing output
```

---

## S4 Method Definitions

Define generics in `R/AllGenerics.R`, methods in the relevant module file.

```r
# R/AllGenerics.R
setGeneric("mx_de", function(study, method, formula, ...) {
  standardGeneric("mx_de")
})

setGeneric("show", function(object) standardGeneric("show"))

# R/03_de.R — method implementation
setMethod("show", "metaXpressStudy", function(object) {
  cat("metaXpressStudy\n")
  cat("  Accession :", object@accession, "\n")
  cat("  Organism  :", object@organism, "\n")
  cat("  Genes     :", nrow(object@counts), "\n")
  cat("  Samples   :", ncol(object@counts), "\n")
  cat("  QC score  :", object@qc_score, "/ 10\n")
  cat("  DE run    :", nrow(object@de_result) > 0, "\n")
})

setMethod("show", "metaXpressResult", function(object) {
  cat("metaXpressResult\n")
  cat("  Method    :", object@method, "\n")
  cat("  Studies   :", object@n_studies, "\n")
  cat("  Genes     :", nrow(object@meta_table), "\n")
  n_sig <- sum(object@meta_table$meta_padj <= 0.05, na.rm = TRUE)
  cat("  Sig genes :", n_sig, "(meta_padj <= 0.05)\n")
})
```

---

## Parallelization

Use `BiocParallel` for any function that loops over studies.

```r
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object
#'   controlling parallelization. Defaults to
#'   \code{\link[BiocParallel]{SerialParam}()} (single-core).
#'   Use \code{\link[BiocParallel]{MulticoreParam}()} for multi-core.
mx_de_all <- function(studies, method = "DESeq2", formula = ~ condition,
                       BPPARAM = BiocParallel::SerialParam()) {
  BiocParallel::bplapply(studies, function(study) {
    mx_de(study, method = method, formula = formula)
  }, BPPARAM = BPPARAM)
}
```

---

## Unit Testing Standards

Use `testthat` (>= 3.0). One test file per module. Tests must be self-contained.

```r
# tests/testthat/test-meta.R

# Setup: create minimal test data (no network calls in tests)
make_test_de_results <- function(n_genes = 100, n_studies = 3) {
  set.seed(42)
  lapply(seq_len(n_studies), function(i) {
    data.frame(
      gene_id = paste0("GENE", seq_len(n_genes)),
      log2FC  = rnorm(n_genes),
      pvalue  = runif(n_genes),
      padj    = p.adjust(runif(n_genes), method = "BH"),
      baseMean = abs(rnorm(n_genes, mean = 100)),
      stringsAsFactors = FALSE
    )
  })
}

test_that("mx_meta returns metaXpressResult S4 object", {
  de <- make_test_de_results()
  result <- mx_meta(de, method = "random_effects")
  expect_s4_class(result, "metaXpressResult")
})

test_that("meta_table has all required columns", {
  de <- make_test_de_results()
  result <- mx_meta(de, method = "fisher")
  required <- c("gene_id", "meta_log2FC", "meta_pvalue", "meta_padj",
                 "i_squared", "n_studies")
  expect_true(all(required %in% colnames(result@meta_table)))
})

test_that("meta_padj values are in [0, 1]", {
  de <- make_test_de_results()
  result <- mx_meta(de, method = "stouffer")
  padj <- result@meta_table$meta_padj
  expect_true(all(padj >= 0 & padj <= 1, na.rm = TRUE))
})

test_that("mx_meta warns with fewer than 3 studies for random_effects", {
  de <- make_test_de_results(n_studies = 2)
  expect_warning(mx_meta(de, method = "random_effects"),
                 "fewer than 3 studies")
})

test_that("mx_meta errors on invalid method", {
  de <- make_test_de_results()
  expect_error(mx_meta(de, method = "invalid_method"))
})
```

### Test coverage targets:
- Every exported function: >= 3 tests
- Edge cases: empty input, single study, all-NA genes
- Error conditions: wrong input types, missing columns
- Numerical correctness: known-answer tests for statistical functions

---

## Bioconductor Submission Checklist

Before submitting to Bioconductor:

```
Package structure:
  [ ] R CMD check passes with 0 errors, 0 warnings, <= 1 note
  [ ] BiocCheck passes (run BiocManager::install("BiocCheck"); BiocCheck::BiocCheck(".")
  [ ] Version is 0.99.x for initial submission
  [ ] biocViews terms are correct and specific

Documentation:
  [ ] All exported functions have complete roxygen2 docs
  [ ] All @examples are runnable (or wrapped in \dontrun{})
  [ ] Vignette builds without errors (R CMD build .)
  [ ] NEWS file exists with version history

Code quality:
  [ ] No T/F (use TRUE/FALSE)
  [ ] No 1:length(x) (use seq_along(x) or seq_len(length(x)))
  [ ] No library() / require() inside functions
  [ ] No direct slot access with @ outside the package (use accessors)
  [ ] vapply() preferred over sapply()
  [ ] No hard-coded file paths

Data:
  [ ] Example data is small (< 1MB uncompressed)
  [ ] data/ objects are documented in R/data.R
  [ ] LazyData: false in DESCRIPTION (Bioconductor requirement)

Testing:
  [ ] Test coverage >= 80% (check with covr::package_coverage())
  [ ] No tests that require internet access (mock or skip_if_offline())
```
