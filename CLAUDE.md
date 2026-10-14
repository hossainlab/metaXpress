# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
in this repository.

## Project Overview

`metaXpress` is an R package (targeting Bioconductor) that provides an
end-to-end pipeline for bulk RNA-seq meta-analysis — from GEO data
ingestion and QC through per-study DE, meta-analysis statistics, pathway
enrichment, visualization, and report generation. The package is
currently in the **scaffolding/implementation phase** — the
authoritative specification is `metaXpress_package_plan.md`.

## R Package Development Commands

``` r

# Load all functions during development (replaces source())
devtools::load_all()

# Regenerate NAMESPACE and .Rd docs from roxygen2 comments
devtools::document()

# Run all tests
devtools::test()

# Run a single test file
devtools::test_file("tests/testthat/test-01_ingest.R")

# Full R CMD check (run before any commit)
devtools::check()

# Bioconductor-specific checks (before submission)
BiocCheck::BiocCheck()

# Build vignettes
devtools::build_vignettes()

# Initial package scaffold (run once)
usethis::create_package(".")
```

## Architecture

The package is organized into 9 numbered modules in `R/`, each with a
corresponding test file in `tests/testthat/`:

| File | Module | Purpose |
|----|----|----|
| `R/01_ingest.R` | Data Ingestion & QC | GEO/SRA fetch, 10-point QC, sample clustering |
| `R/02_harmonize.R` | Normalization & Harmonization | Gene ID reannotation, library-type correction, batch removal |
| `R/03_de.R` | Per-Study DE | DESeq2/edgeR/limma-voom wrappers |
| `R/04_meta.R` | Meta-Analysis Statistics | P-value combination, effect-size models, heterogeneity |
| `R/05_missing.R` | Missing Gene Handling | Coverage reporting, imputation strategies |
| `R/06_pathway.R` | Pathway Meta-Analysis | Cross-study ORA/GSEA, redundancy reduction |
| `R/07_visualize.R` | Visualization | Volcano, forest, heatmap, UpSet plots |
| `R/08_report.R` | Report Generation | HTML/PDF rendering, CSV/Excel export |
| `R/09_downstream.R` | Downstream Biology | TF activity, cell composition, consensus co-expression, drug reversal |
| `R/utils.R` | Shared Utilities | Helpers used across modules |

## S4 Data Structures

Two core S4 classes flow through the pipeline:

**`metaXpressStudy`** — output of Module 1, input to Modules 2–3.
Extends `SummarizedExperiment`: - `assay(study, "counts")`: raw count
matrix (genes × samples) - `colData(study)`: sample metadata (must
contain `condition`); `study$condition` works - `rowData(study)`:
optional gene annotation (e.g. `gene_length` for TPM) - extra slots:
`accession`, `organism`, `qc_score`, `de_result` - Build with
[`metaXpressStudy()`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy.md);
convert any SE/DESeqDataSet with
[`mx_from_se()`](https://mdjubayerhossain.com/metaXpress/reference/mx_from_se.md)
or `as(se, "metaXpressStudy")` - Never use `study@metadata` — on SE that
is the free-form list, not sample data. Internally use `.counts()`,
`.samples()`, `.set_counts()` (`R/utils.R`)

**`metaXpressResult`** — output of Module 4, input to Modules 5–8: -
`meta_table`: gene-level results (meta-p, meta-FDR, meta-effect size,
I²) - `method`, `n_studies`, `heterogeneity`, `pathway_result`

## Function Naming

All public functions are prefixed `mx_`. Examples:
[`mx_fetch_geo()`](https://mdjubayerhossain.com/metaXpress/reference/mx_fetch_geo.md),
[`mx_de_all()`](https://mdjubayerhossain.com/metaXpress/reference/mx_de_all.md),
[`mx_meta()`](https://mdjubayerhossain.com/metaXpress/reference/mx_meta.md),
[`mx_volcano()`](https://mdjubayerhossain.com/metaXpress/reference/mx_volcano.md).

## Key Design Decisions

- **All methods take and return the standard S4 objects** — functions
  chain cleanly:
  [`mx_fetch_geo()`](https://mdjubayerhossain.com/metaXpress/reference/mx_fetch_geo.md)
  →
  [`mx_filter_studies()`](https://mdjubayerhossain.com/metaXpress/reference/mx_filter_studies.md)
  →
  [`mx_reannotate()`](https://mdjubayerhossain.com/metaXpress/reference/mx_reannotate.md)
  →
  [`mx_de_all()`](https://mdjubayerhossain.com/metaXpress/reference/mx_de_all.md)
  →
  [`mx_meta()`](https://mdjubayerhossain.com/metaXpress/reference/mx_meta.md).
- **Per-study DE is always run independently** on each study before
  meta-analysis (not pooled).
- **Missing genes across studies** are handled explicitly in Module 5
  (imputation or coverage filtering), not silently dropped.
- **Meta-analysis methods** supported: Fisher, Stouffer, inverse-normal
  (fused), fixed effects, random effects (DerSimonian-Laird), AWmeta
  adaptive weighting — selectable via `method` argument.
- **Parallelization** is via `BiocParallel` throughout (not `parallel`
  or `foreach`).
- **Indentation:** 2 spaces (per `.Rproj` settings).

## Dependencies

Required Bioconductor packages: `GEOquery`, `DESeq2`, `edgeR`, `limma`,
`sva`, `tximport`, `clusterProfiler`, `msigdbr`, `ComplexHeatmap`,
`BiocParallel`. Required CRAN packages: `ggplot2`. Suggested: `harmony`,
`AnnotationDbi`, `org.Hs.eg.db`, `rmarkdown`, `openxlsx`.

## Development Roadmap

See `metaXpress_package_plan.md` for the full phased plan. Current phase
targets: S4 class definitions, Module 1 (GEO ingestion + QC), Module 3
(DESeq2 DE wrapper), unit tests for both, and a quickstart vignette.

## Validation Reference

Statistical correctness is validated against reference implementations
in tests:
[`mx_meta()`](https://mdjubayerhossain.com/metaXpress/reference/mx_meta.md)
against
[`metafor::rma()`](https://wviechtb.github.io/metafor/reference/rma.uni.html)
(DL/REML/PM, HKSJ) and
[`mx_de()`](https://mdjubayerhossain.com/metaXpress/reference/mx_de.md)
against native DESeq2 on `airway` (GSE52778). Heberle et al. (2024,
bioRxiv, <doi:10.1101/2024.11.07.622520>) is a systematic review of AD
brain RNA-seq whose 10-category tool grades study *reporting* quality;
the numeric thresholds in
[`mx_qc_study()`](https://mdjubayerhossain.com/metaXpress/reference/mx_qc_study.md)
are metaXpress defaults, not taken from that paper. Never cite a
reference without verifying it exists (PubMed/DOI).
