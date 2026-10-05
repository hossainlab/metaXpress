# Package index

## Package Overview

Top-level package documentation and architecture.

- [`metaXpress`](https://hossainlab.github.io/metaXpress/reference/metaXpress-package.md)
  [`metaXpress-package`](https://hossainlab.github.io/metaXpress/reference/metaXpress-package.md)
  : metaXpress: End-to-End Bulk RNA-seq Meta-Analysis

## S4 Data Structures

Core container classes for studies and meta-analysis results.

- [`metaXpressStudy-class`](https://hossainlab.github.io/metaXpress/reference/metaXpressStudy-class.md)
  : The metaXpressStudy class
- [`metaXpressResult-class`](https://hossainlab.github.io/metaXpress/reference/metaXpressResult-class.md)
  : The metaXpressResult class
- [`show(`*`<metaXpressStudy>`*`)`](https://hossainlab.github.io/metaXpress/reference/show-metaXpressStudy-method.md)
  : Display a metaXpressStudy object
- [`show(`*`<metaXpressResult>`*`)`](https://hossainlab.github.io/metaXpress/reference/show-metaXpressResult-method.md)
  : Display a metaXpressResult object

## Module 1: Ingestion & Quality Control

Fetch studies from GEO/SRA, evaluate 10-point QC criteria, and filter
datasets.

- [`mx_fetch_geo()`](https://hossainlab.github.io/metaXpress/reference/mx_fetch_geo.md)
  : Fetch bulk RNA-seq studies from GEO
- [`mx_fetch_sra()`](https://hossainlab.github.io/metaXpress/reference/mx_fetch_sra.md)
  : Fetch raw RNA-seq data from SRA
- [`mx_load_local()`](https://hossainlab.github.io/metaXpress/reference/mx_load_local.md)
  : Load user-supplied count matrices as metaXpressStudy objects
- [`mx_qc_study()`](https://hossainlab.github.io/metaXpress/reference/mx_qc_study.md)
  : Score a study against the 10-point QC checklist
- [`mx_filter_studies()`](https://hossainlab.github.io/metaXpress/reference/mx_filter_studies.md)
  : Filter studies below a QC threshold
- [`mx_cluster_samples()`](https://hossainlab.github.io/metaXpress/reference/mx_cluster_samples.md)
  : Automatically cluster samples by metadata
- [`mx_study_overview()`](https://hossainlab.github.io/metaXpress/reference/mx_study_overview.md)
  : Study characteristics overview plot

## Module 2: Normalization & Cross-Study Harmonization

Reannotate gene IDs, harmonize library types, correct batch effects, and
align gene spaces.

- [`mx_reannotate()`](https://hossainlab.github.io/metaXpress/reference/mx_reannotate.md)
  : Reannotate gene IDs to a common namespace
- [`mx_normalize()`](https://hossainlab.github.io/metaXpress/reference/mx_normalize.md)
  : Normalize counts within a single study
- [`mx_correct_library_type()`](https://hossainlab.github.io/metaXpress/reference/mx_correct_library_type.md)
  : Correct for library type bias (polyA vs rRNA-depleted)
- [`mx_remove_batch()`](https://hossainlab.github.io/metaXpress/reference/mx_remove_batch.md)
  : Remove batch effects across studies
- [`mx_align_genes()`](https://hossainlab.github.io/metaXpress/reference/mx_align_genes.md)
  : Restrict all studies to a common set of genes

## Module 3: Per-Study Differential Expression

Perform independent differential expression using DESeq2, edgeR, or
limma-voom.

- [`mx_de()`](https://hossainlab.github.io/metaXpress/reference/mx_de.md)
  : Run differential expression analysis on a single study
- [`mx_de_all()`](https://hossainlab.github.io/metaXpress/reference/mx_de_all.md)
  : Run differential expression on all studies
- [`mx_de_summary()`](https://hossainlab.github.io/metaXpress/reference/mx_de_summary.md)
  : Summarise per-study DE results

## Module 4: Meta-Analysis Statistics

P-value combination (Fisher, Stouffer), effect-size meta-analysis (REM,
FEM), and heterogeneity.

- [`mx_meta()`](https://hossainlab.github.io/metaXpress/reference/mx_meta.md)
  : Run meta-analysis across per-study DE results
- [`mx_heterogeneity()`](https://hossainlab.github.io/metaXpress/reference/mx_heterogeneity.md)
  : Compute per-gene heterogeneity statistics
- [`mx_sensitivity()`](https://hossainlab.github.io/metaXpress/reference/mx_sensitivity.md)
  : Leave-one-out sensitivity analysis

## Module 5: Missing Gene Handling

Identify missingness patterns across studies and apply imputation or
coverage filtering.

- [`mx_impute()`](https://hossainlab.github.io/metaXpress/reference/mx_impute.md)
  : Impute missing gene statistics across studies
- [`mx_filter_coverage()`](https://hossainlab.github.io/metaXpress/reference/mx_filter_coverage.md)
  : Filter genes by minimum study coverage
- [`mx_missing_summary()`](https://hossainlab.github.io/metaXpress/reference/mx_missing_summary.md)
  : Summarise gene coverage across studies

## Module 6: Pathway Meta-Analysis

Cross-study over-representation (ORA) and GSEA, consensus pathways, and
deduplication.

- [`mx_pathway_meta()`](https://hossainlab.github.io/metaXpress/reference/mx_pathway_meta.md)
  : Perform pathway meta-analysis on meta-analysis results
- [`mx_pathway_consensus()`](https://hossainlab.github.io/metaXpress/reference/mx_pathway_consensus.md)
  : Find consensus pathways across studies
- [`mx_pathway_dedup()`](https://hossainlab.github.io/metaXpress/reference/mx_pathway_dedup.md)
  : Remove redundant pathways
- [`mx_pathway_heatmap()`](https://hossainlab.github.io/metaXpress/reference/mx_pathway_heatmap.md)
  : Cross-study pathway heatmap

## Module 7: Visualization

Publication-grade plots including volcano plots, forest plots, heatmaps,
and UpSet diagrams.

- [`mx_volcano()`](https://hossainlab.github.io/metaXpress/reference/mx_volcano.md)
  : Meta-analysis volcano plot
- [`mx_forest()`](https://hossainlab.github.io/metaXpress/reference/mx_forest.md)
  : Forest plot for a single gene
- [`mx_heatmap()`](https://hossainlab.github.io/metaXpress/reference/mx_heatmap.md)
  : Heatmap of top DEGs across studies
- [`mx_upset()`](https://hossainlab.github.io/metaXpress/reference/mx_upset.md)
  : UpSet plot of DEG overlap across studies
- [`mx_heterogeneity_plot()`](https://hossainlab.github.io/metaXpress/reference/mx_heterogeneity_plot.md)
  : I-squared distribution plot

## Module 8: Reporting & Export

Generate reproducible HTML/PDF reports and export meta-analysis results.

- [`mx_report()`](https://hossainlab.github.io/metaXpress/reference/mx_report.md)
  : Generate a reproducible analysis report
- [`mx_export()`](https://hossainlab.github.io/metaXpress/reference/mx_export.md)
  : Export meta-analysis results to file
- [`mx_session_info()`](https://hossainlab.github.io/metaXpress/reference/mx_session_info.md)
  : Capture session information for reproducibility
- [`mx_run_app()`](https://hossainlab.github.io/metaXpress/reference/mx_run_app.md)
  : Launch the metaXpress Interactive Explorer Shiny Application
