# metaXpress Full Workflow: GEO to Report

## Overview

`metaXpress` is an end-to-end R package designed for multi-study
transcriptomic meta-analysis. It provides a standardized pipeline
covering data ingestion from NCBI GEO and SRA, 10-point quality control
scoring, cross-platform harmonization, per-study differential
expression, six statistical meta-analysis models, missing gene handling,
pathway enrichment, interactive visualization, and reproducible
reporting.

For full documentation and online tutorials, visit the official website
at <https://mdjubayerhossain.com/metaXpress/>.

## Installation from GitHub

`metaXpress` is currently available on GitHub and targets future
submission to Bioconductor.

``` r

# 1. Ensure BiocManager is present for Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# 2. Install remotes
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

# 3. Install metaXpress
remotes::install_github("hossainlab/metaXpress", dependencies = TRUE)
```

Load the package:

``` r

library(metaXpress)
```

------------------------------------------------------------------------

## 1. Data Ingestion & Quality Control

### 1a. Fetch Public Data from GEO

You can ingest multiple public RNA-seq datasets directly using their
NCBI GEO accessions:

``` r

# Fetch three bulk RNA-seq cohorts from GEO
accessions <- c("GSE130688", "GSE136569", "GSE171485")
studies    <- mx_fetch_geo(accessions, count_type = "raw", cache_dir = "geo_cache/")
```

Alternatively, load local count matrices and metadata with
[`mx_load_local()`](https://mdjubayerhossain.com/metaXpress/reference/mx_load_local.md).

### 1b. Inspect 10-Point QC Scoring

[`mx_qc_study()`](https://mdjubayerhossain.com/metaXpress/reference/mx_qc_study.md)
benchmarks each study against 10 criteria (sample size, depth, alignment
rate, duplicate rate, case/control completeness, integer counts, etc.):

``` r

qc_scores <- vapply(studies, function(s) s@qc_score, numeric(1))
print(qc_scores)

# Inspect detailed metrics for study 1
attr(studies[[1]], "qc_details")
```

### 1c. Filter Low-Quality Cohorts

``` r

# Filter out cohorts failing the QC threshold (minimum score: 7/10)
studies <- mx_filter_studies(studies, qc_threshold = 7)
```

------------------------------------------------------------------------

## 2. Cross-Study Harmonization

### 2a. Reannotate Gene Identifiers

Convert diverse gene identifiers (Ensembl, Entrez, RefSeq) to a uniform
namespace:

``` r

studies <- mx_reannotate(studies, org = "Homo sapiens", target_id = "SYMBOL")
```

### 2b. Correct Library Type Bias

If cohorts mix poly(A)-selected and rRNA-depleted protocols, adjust for
library type biases:

``` r

studies <- mx_correct_library_type(studies)
```

### 2c. Remove Batch Effects

``` r

studies <- mx_remove_batch(studies, method = "ComBat-seq")
```

### 2d. Align Common Gene Space

``` r

studies <- mx_align_genes(studies)
```

------------------------------------------------------------------------

## 3. Per-Study Differential Expression

Run independent statistical tests on each cohort using DESeq2, edgeR, or
limma-voom:

``` r

studies <- mx_de_all(
  studies,
  method  = "DESeq2",
  formula = ~ condition,
  BPPARAM = BiocParallel::MulticoreParam(4)
)
```

Inspect significant gene counts:

``` r

mx_de_summary(studies, padj_threshold = 0.05, lfc_threshold = 1.0)
```

------------------------------------------------------------------------

## 4. Handle Missing Genes

When cohorts have varying gene detection rates:

``` r

de_results <- lapply(studies, function(s) s@de_result)

# Inspect gene coverage across cohorts
cov_mat <- mx_missing_summary(de_results)

# Retain genes present in at least 2 cohorts or apply k-NN imputation
de_results <- mx_filter_coverage(de_results, min_studies = 2)
# de_results <- mx_impute(de_results, method = "knn")
```

------------------------------------------------------------------------

## 5. Statistical Meta-Analysis

### 5a. Model Selection Guide

- **2 cohorts:** `"fisher"` (Fisher’s combined probability) or
  `"stouffer"` (Z-score weighting).
- **3+ cohorts, low heterogeneity ($`I^2 < 25\%`$):** `"fixed_effects"`
  (Inverse-variance).
- **3+ cohorts, moderate-to-high heterogeneity ($`I^2 \ge 25\%`$):**
  `"random_effects"` (DerSimonian-Laird, recommended default).
- **Mixed heterogeneity:** `"awmeta"` (Adaptive weighting).

``` r

meta_result <- mx_meta(de_results, method = "random_effects", min_studies = 2)
meta_result
```

### 5b. Quantify Between-Study Heterogeneity

Assess Cochran’s $`Q`$, Higgins $`I^2`$, and $`\tau^2`$ for every gene:

``` r

het <- mx_heterogeneity(de_results)
summary(het$I_sq)
```

### 5c. Leave-One-Out Sensitivity Analysis

Ensure discoveries are not driven by an individual outlier cohort:

``` r

loo <- mx_sensitivity(de_results, method = "random_effects")
vapply(loo, function(r) sum(r@meta_table$meta_padj <= 0.05, na.rm = TRUE), integer(1))
```

------------------------------------------------------------------------

## 6. Pathway Meta-Analysis

Perform cross-study Over-Representation Analysis (ORA) or Gene Set
Enrichment Analysis (GSEA) on Hallmark, KEGG, or GO sets:

``` r

meta_result <- mx_pathway_meta(meta_result, db = "Hallmarks", method = "ORA")
head(meta_result@pathway_result[order(meta_result@pathway_result$padj), ], 10)
```

------------------------------------------------------------------------

## 7. Publication Visualizations

### Volcano Plot

``` r

mx_volcano(meta_result, padj_threshold = 0.05, lfc_threshold = 1.0, label_top = 15)
```

### Forest Plot

``` r

top_gene <- meta_result@meta_table$gene_id[which.min(meta_result@meta_table$meta_padj)]
mx_forest(top_gene, de_results, meta_result)
```

### Heterogeneity Plot

``` r

mx_heterogeneity_plot(meta_result)
```

------------------------------------------------------------------------

## 8. Interactive Explorer & Report Generation

### Launch Interactive Shiny App

Explore results dynamically in your browser:

``` r

metaXpress::mx_run_app()
```

### Generate Reproducible Reports & Export

``` r

mx_report(meta_result, studies, de_results, format = "html", output_dir = "results/")
mx_export(meta_result, format = "csv", output_dir = "results/")
```

------------------------------------------------------------------------

## Session Information

``` r

mx_session_info()$session_info
```
