# metaXpress Quickstart

## Overview

`metaXpress` provides an end-to-end pipeline for integrating multiple
bulk RNA-seq studies into a single meta-analysis result. This vignette
walks through a minimal example using the built-in simulated dataset.

``` r

library(metaXpress)
```

## Ingest Public Data from GEO

Users can fetch public RNA-seq count matrices and metadata directly from
NCBI GEO using accession IDs:

``` r

# Fetch two bulk RNA-seq cohorts from GEO
studies <- mx_fetch_geo(c("GSE53697", "GSE95587"))
lapply(studies, show)
```

## Step 1 — QC Filtering

Studies are already QC-scored in the example data. In a real workflow,
QC scores are computed by
[`mx_fetch_geo()`](https://hossainlab.github.io/metaXpress/reference/mx_fetch_geo.md)
automatically, or manually via
[`mx_qc_study()`](https://hossainlab.github.io/metaXpress/reference/mx_qc_study.md).

``` r

studies <- mx_filter_studies(studies, qc_threshold = 7)
message("Studies passing QC: ", length(studies))
```

## Step 2 — Align Genes

Find the common gene universe across all studies.

``` r

studies <- mx_align_genes(studies)
```

## Step 3 — Per-Study Differential Expression

Run DESeq2 on each study independently.

``` r

studies <- mx_de_all(studies, method = "DESeq2", formula = ~ condition)
mx_de_summary(studies)
```

## Step 4 — Meta-Analysis

Combine per-study DE results using the random effects model.

``` r

de_results <- lapply(studies, function(s) s@de_result)
meta_result <- mx_meta(de_results, method = "random_effects")
meta_result
```

## Step 5 — Visualise

``` r

mx_volcano(meta_result, label_top = 10)
```

``` r

# Forest plot for the top gene
top_gene <- meta_result@meta_table$gene_id[
  which.min(meta_result@meta_table$meta_padj)]
mx_forest(top_gene, de_results, meta_result)
```

``` r

mx_heterogeneity_plot(meta_result)
```

## Step 6 — Export Results

``` r

mx_export(meta_result, format = "csv", output_dir = tempdir())
```

## Session Info

``` r

mx_session_info()$session_info
```
