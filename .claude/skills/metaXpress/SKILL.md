---
name: metaxpress-rnaseq-meta-analysis
description: Develops the metaXpress R package for end-to-end bulk RNA-seq meta-analysis. Use when user asks to build, scaffold, implement, or extend any module of the metaXpress package, including data ingestion from GEO, normalization, per-study differential expression, meta-analysis statistics, missing gene handling, pathway enrichment, visualization, or report generation. Also use when user asks to write tests, vignettes, or DESCRIPTION file for metaXpress.
license: MIT
metadata:
  author: DeepBio Academy
  version: 1.0.0
  category: bioinformatics
  tags: [rnaseq, meta-analysis, bioconductor, R-package, transcriptomics]
---

# metaXpress R Package Development Skill

## Overview

metaXpress is an end-to-end R package for bulk RNA-seq meta-analysis, targeting Bioconductor. It integrates multiple independent studies from GEO/SRA through harmonization, per-study DE analysis, meta-analysis statistics, pathway enrichment, and reproducible reporting.

**Function naming convention:** All exported functions use the `mx_` prefix.
**S4 classes:** `metaXpressStudy`, `metaXpressResult`
**Target:** Bioconductor (primary), CRAN (secondary)
**License:** MIT

Consult `references/architecture.md` for the full module map and `references/statistical_methods.md` for meta-analysis method details before implementing any statistical functions.

---

## Instructions

### CRITICAL: Before Writing Any Code

1. Check which module is being implemented (see Module Map below)
2. Read the relevant section in `references/architecture.md`
3. For statistical functions, read `references/statistical_methods.md`
4. Follow the R/Bioconductor coding standards in `references/bioconductor_standards.md`
5. Always use `padj` (not `pvalue`) for significance filtering in DE results

---

### Module Map

| Module | File | Key Functions |
|--------|------|---------------|
| 1. Data Ingestion | `R/01_ingest.R` | `mx_fetch_geo`, `mx_load_local`, `mx_qc_study`, `mx_filter_studies` |
| 2. Harmonization | `R/02_harmonize.R` | `mx_reannotate`, `mx_normalize`, `mx_remove_batch`, `mx_align_genes` |
| 3. Per-Study DE | `R/03_de.R` | `mx_de`, `mx_de_all` |
| 4. Meta-Analysis | `R/04_meta.R` | `mx_meta`, `mx_heterogeneity`, `mx_sensitivity` |
| 5. Missing Genes | `R/05_missing.R` | `mx_missing_summary`, `mx_impute`, `mx_filter_coverage` |
| 6. Pathway | `R/06_pathway.R` | `mx_pathway_meta`, `mx_pathway_consensus` |
| 7. Visualization | `R/07_visualize.R` | `mx_volcano`, `mx_forest`, `mx_heatmap`, `mx_upset` |
| 8. Reporting | `R/08_report.R` | `mx_report`, `mx_export` |

---

### Step 1: Package Scaffold

When asked to scaffold or initialize the package:

```r
# Run in R console - creates the full package skeleton
usethis::create_package("metaXpress")
usethis::use_mit_license()
usethis::use_bioc_description()
usethis::use_testthat()
usethis::use_vignette("quickstart")
usethis::use_vignette("full_workflow")
usethis::use_git()

# Create module files
module_files <- c(
  "R/01_ingest.R", "R/02_harmonize.R", "R/03_de.R",
  "R/04_meta.R", "R/05_missing.R", "R/06_pathway.R",
  "R/07_visualize.R", "R/08_report.R", "R/utils.R",
  "R/AllClasses.R", "R/AllGenerics.R"
)
invisible(lapply(module_files, file.create))
```

Expected output: A valid R package directory with DESCRIPTION, NAMESPACE, R/, tests/, vignettes/

---

### Step 2: S4 Class Definitions

Always define classes in `R/AllClasses.R` before implementing functions.

```r
# R/AllClasses.R
setClass("metaXpressStudy",
  representation(
    counts    = "matrix",       # Raw count matrix (genes x samples)
    metadata  = "data.frame",   # Sample metadata with 'condition' column
    accession = "character",    # GEO/SRA accession ID
    organism  = "character",    # e.g., "Homo sapiens"
    qc_score  = "numeric",      # 0-10 QC score
    de_result = "data.frame"    # Populated by mx_de()
  ),
  validity = function(object) {
    if (!("condition" %in% colnames(object@metadata)))
      return("metadata must contain a 'condition' column")
    if (ncol(object@counts) != nrow(object@metadata))
      return("counts columns must match metadata rows")
    TRUE
  }
)

setClass("metaXpressResult",
  representation(
    meta_table     = "data.frame",  # Gene-level meta-analysis results
    method         = "character",   # Method used
    n_studies      = "integer",     # Number of studies
    heterogeneity  = "data.frame",  # I-squared and Q-stat per gene
    pathway_result = "data.frame"   # Pathway meta-analysis results
  )
)
```

---

### Step 3: Module Implementation Guidelines

#### Module 1 — Data Ingestion (R/01_ingest.R)

```r
#' Fetch studies from GEO
#' @param accessions Character vector of GEO accession IDs (e.g., "GSE12345")
#' @param count_type One of "raw", "normalized". Always prefer "raw".
#' @return List of metaXpressStudy objects
mx_fetch_geo <- function(accessions, count_type = "raw") {
  # Use GEOquery::getGEO() for metadata
  # Use GEOquery::getGEOSuppFiles() for count matrices
  # Apply mx_qc_study() to each fetched study
  # Return named list of metaXpressStudy objects
}
```

QC scoring (10 criteria, 1 point each — see `references/qc_criteria.md`):
- Minimum n >= 3 per group
- Sequencing depth >= 10M reads
- Alignment rate >= 70%
- rRNA contamination < 10%
- Duplicate rate < 50%
- Gene detection >= 15,000 genes
- Metadata completeness
- Clear case/control definition
- No batch-condition confounding
- Raw counts available

#### Module 3 — Per-Study DE (R/03_de.R)

CRITICAL rules for DE analysis:
- ALWAYS use `padj` (adjusted p-value) for significance, NEVER raw `pvalue`
- Default significance threshold: `padj <= 0.05` AND `|log2FoldChange| >= 1`
- Use inclusive inequalities (`<=`, `>=`)
- Store full results (all genes), not just significant ones

```r
mx_de <- function(study, method = c("DESeq2", "edgeR", "limma-voom"),
                  formula = ~ condition, ...) {
  method <- match.arg(method)
  # Dispatch to internal method-specific function
  # Return standardized data.frame with columns:
  # gene_id, log2FC, pvalue, padj, baseMean, method
}
```

#### Module 4 — Meta-Analysis (R/04_meta.R)

Supported methods (see `references/statistical_methods.md` for formulas):

| method argument | Statistical approach | Reference |
|---|---|---|
| `"fisher"` | Fisher's combined p-value | Rau et al. 2013 |
| `"stouffer"` | Stouffer's Z-score (weighted) | Rau et al. 2013 |
| `"inverse_normal"` | Fused inverse-normal | Prasad & Li 2021 |
| `"fixed_effects"` | Fixed effects (effect size) | Keel & Lindholm-Perry 2022 |
| `"random_effects"` | DerSimonian-Laird random effects | Keel & Lindholm-Perry 2022 |
| `"awmeta"` | Adaptive weighting | Hu et al. 2025 |

```r
mx_meta <- function(de_results, method = "random_effects",
                    min_studies = 2, ...) {
  # de_results: list of data.frames from mx_de_all()
  # Align genes across studies first (call mx_align_genes internally)
  # Apply chosen method
  # Compute I-squared and Q-statistic for heterogeneity
  # Return metaXpressResult object
  # Output columns: gene_id, meta_log2FC, meta_pvalue, meta_padj,
  #                 i_squared, q_stat, n_studies, direction_consistency
}
```

CRITICAL: Meta-analysis p-values must also be adjusted (BH/FDR). Store as `meta_padj`.

#### Module 7 — Visualization (R/07_visualize.R)

- Use `ggplot2` with `theme_bw()` as default theme
- Use colorblind-friendly palettes (`RColorBrewer` "Set2" or `viridis`)
- Forest plots must show per-study effect sizes + CI + pooled estimate + I-squared
- Avoid overlapping labels; use `ggrepel::geom_label_repel()` for gene labels
- Export both `.svg` and `.png` when saving

```r
mx_forest <- function(gene, de_results, meta_result = NULL) {
  # gene: single gene ID string
  # de_results: list of per-study DE results
  # Shows: per-study log2FC with 95% CI, pooled estimate, I-squared
  # Returns: ggplot object
}
```

---

### Step 4: Writing Tests

Place tests in `tests/testthat/`. One test file per module.

```r
# tests/testthat/test-meta.R
test_that("mx_meta returns metaXpressResult with correct columns", {
  # Use built-in example data: data(metaXpress_example)
  result <- mx_meta(example_de_results, method = "random_effects")
  expect_s4_class(result, "metaXpressResult")
  expect_true(all(c("gene_id", "meta_padj", "i_squared") %in%
                    colnames(result@meta_table)))
  expect_true(all(result@meta_table$meta_padj >= 0, na.rm = TRUE))
  expect_true(all(result@meta_table$meta_padj <= 1, na.rm = TRUE))
})

test_that("mx_meta uses padj not pvalue for significance", {
  # Verify the output contains meta_padj column, not just meta_pvalue
  result <- mx_meta(example_de_results, method = "fisher")
  expect_true("meta_padj" %in% colnames(result@meta_table))
})
```

---

### Step 5: DESCRIPTION File

```
Package: metaXpress
Title: End-to-End Bulk RNA-seq Meta-Analysis
Version: 0.99.0
Authors@R: person("DeepBio Academy", role = c("aut", "cre"),
    email = "contact@deepbioacademy.com")
Description: metaXpress provides a unified pipeline for integrating
    multiple bulk RNA-seq studies. It covers data ingestion from GEO/SRA,
    cross-study normalization and batch correction, per-study differential
    expression analysis, meta-analysis statistics (p-value combination and
    effect size models), missing gene handling, pathway meta-analysis, and
    reproducible HTML/PDF reporting.
License: MIT + file LICENSE
Encoding: UTF-8
LazyData: false
Depends: R (>= 4.3.0)
Imports:
    GEOquery,
    DESeq2,
    edgeR,
    limma,
    sva,
    tximport,
    clusterProfiler,
    msigdbr,
    ggplot2,
    ggrepel,
    ComplexHeatmap,
    BiocParallel,
    AnnotationDbi,
    org.Hs.eg.db,
    methods,
    S4Vectors
Suggests:
    testthat (>= 3.0.0),
    knitr,
    rmarkdown,
    BiocStyle,
    openxlsx
VignetteBuilder: knitr
biocViews: RNASeq, DifferentialExpression, GeneExpression,
    Transcriptomics, BatchEffect, Normalization, Pathways,
    MultipleComparison, ReportWriting
BugReports: https://github.com/deepbioacademy/metaXpress/issues
URL: https://github.com/deepbioacademy/metaXpress
```

---

## Examples

### Example 1: Scaffold the full package

User says: "Set up the metaXpress package skeleton"

Actions:
1. Run `usethis::create_package("metaXpress")`
2. Create all module files in `R/`
3. Define S4 classes in `R/AllClasses.R`
4. Write `DESCRIPTION` file
5. Set up `tests/testthat/`

Result: Valid R package directory ready for Bioconductor submission workflow

---

### Example 2: Implement a specific module

User says: "Implement the meta-analysis statistics module"

Actions:
1. Read `references/statistical_methods.md` for method formulas
2. Implement `mx_meta()` with all 6 method options
3. Implement `mx_heterogeneity()` for I-squared computation
4. Implement `mx_sensitivity()` for leave-one-out analysis
5. Write unit tests in `tests/testthat/test-meta.R`

Result: Fully tested `R/04_meta.R` with roxygen2 documentation

---

### Example 3: Add a visualization

User says: "Create a forest plot function"

Actions:
1. Implement `mx_forest()` in `R/07_visualize.R`
2. Use ggplot2 with per-study CI bars + pooled diamond + I-squared annotation
3. Use colorblind-safe palette
4. Add `ggrepel` for gene labels
5. Return ggplot object (do not save by default)

Result: `mx_forest()` function with roxygen2 docs and a test

---

## Common Issues

### Error: "condition column missing from metadata"
Cause: metaXpressStudy validity check failed
Solution: Ensure metadata data.frame has a column named exactly `condition` with case/control labels

### Error: "padj contains NA for all genes"
Cause: DESeq2 independent filtering removed all genes (too few samples)
Solution: Set `independentFiltering = FALSE` in DESeq2 results call, or increase sample size threshold in `mx_qc_study()`

### Error: "Cannot align genes — no common genes across studies"
Cause: Gene ID namespace mismatch (Ensembl vs Symbol vs Entrez)
Solution: Run `mx_reannotate()` before `mx_align_genes()` to standardize to a common namespace

### Warning: "Fewer than 3 studies — meta-analysis results unreliable"
Cause: Too few studies for stable random effects estimation
Solution: Use `method = "fisher"` or `method = "stouffer"` for 2-study analyses; random effects requires >= 3 studies

---

## Performance Notes

- Take time to implement each function thoroughly before moving to the next
- Quality and correctness are more important than speed
- Do not skip validity checks in S4 class definitions
- Always write at least 3 unit tests per exported function
