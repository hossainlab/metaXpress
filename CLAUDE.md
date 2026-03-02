# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working in this repository.

## Project Overview

`metaXpress` is an R package (targeting Bioconductor) that provides an end-to-end pipeline for bulk RNA-seq meta-analysis — from GEO data ingestion and QC through per-study DE, meta-analysis statistics, pathway enrichment, visualization, and report generation. The package is currently in the **scaffolding/implementation phase** — the authoritative specification is `metaXpress_package_plan.md`.

## R Package Development Commands

```r
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

The package is organized into 8 numbered modules in `R/`, each with a corresponding test file in `tests/testthat/`:

| File | Module | Purpose |
|---|---|---|
| `R/01_ingest.R` | Data Ingestion & QC | GEO/SRA fetch, 10-point QC, sample clustering |
| `R/02_harmonize.R` | Normalization & Harmonization | Gene ID reannotation, library-type correction, batch removal |
| `R/03_de.R` | Per-Study DE | DESeq2/edgeR/limma-voom wrappers |
| `R/04_meta.R` | Meta-Analysis Statistics | P-value combination, effect-size models, heterogeneity |
| `R/05_missing.R` | Missing Gene Handling | Coverage reporting, imputation strategies |
| `R/06_pathway.R` | Pathway Meta-Analysis | Cross-study ORA/GSEA, redundancy reduction |
| `R/07_visualize.R` | Visualization | Volcano, forest, heatmap, UpSet plots |
| `R/08_report.R` | Report Generation | HTML/PDF rendering, CSV/Excel export |
| `R/utils.R` | Shared Utilities | Helpers used across modules |

## S4 Data Structures

Two core S4 classes flow through the pipeline:

**`metaXpressStudy`** — output of Module 1, input to Modules 2–3:
- `counts`: raw count matrix (genes × samples)
- `metadata`: sample data.frame
- `accession`, `organism`, `qc_score`, `de_result`

**`metaXpressResult`** — output of Module 4, input to Modules 5–8:
- `meta_table`: gene-level results (meta-p, meta-FDR, meta-effect size, I²)
- `method`, `n_studies`, `heterogeneity`, `pathway_result`

## Function Naming

All public functions are prefixed `mx_`. Examples: `mx_fetch_geo()`, `mx_de_all()`, `mx_meta()`, `mx_volcano()`.

## Key Design Decisions

- **All methods take and return the standard S4 objects** — functions chain cleanly: `mx_fetch_geo()` → `mx_filter_studies()` → `mx_reannotate()` → `mx_de_all()` → `mx_meta()`.
- **Per-study DE is always run independently** on each study before meta-analysis (not pooled).
- **Missing genes across studies** are handled explicitly in Module 5 (imputation or coverage filtering), not silently dropped.
- **Meta-analysis methods** supported: Fisher, Stouffer, inverse-normal (fused), fixed effects, random effects (DerSimonian-Laird), AWmeta adaptive weighting — selectable via `method` argument.
- **Parallelization** is via `BiocParallel` throughout (not `parallel` or `foreach`).
- **Indentation:** 2 spaces (per `.Rproj` settings).

## Dependencies

Required Bioconductor packages: `GEOquery`, `DESeq2`, `edgeR`, `limma`, `sva`, `tximport`, `clusterProfiler`, `msigdbr`, `ComplexHeatmap`, `BiocParallel`.
Required CRAN packages: `ggplot2`.
Suggested: `harmony`, `AnnotationDbi`, `org.Hs.eg.db`, `rmarkdown`, `openxlsx`.

## Development Roadmap

See `metaXpress_package_plan.md` for the full phased plan. Current phase targets: S4 class definitions, Module 1 (GEO ingestion + QC), Module 3 (DESeq2 DE wrapper), unit tests for both, and a quickstart vignette.

## Validation Reference

Use the Alzheimer's dataset from Heberle et al. (2025) as the primary benchmark for correctness. The QC checklist (10 criteria) in `mx_qc_study()` is directly derived from that paper.
