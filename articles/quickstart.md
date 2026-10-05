# metaXpress Quickstart

## Overview

`metaXpress` provides an end-to-end pipeline for integrating multiple
bulk RNA-seq studies into a unified meta-analysis. This guide walks
through installation from GitHub and a minimal workflow for ingesting
public data from NCBI GEO, performing quality control, differential
expression analysis, meta-analysis, and visualization.

## Installation from GitHub

`metaXpress` is currently available on GitHub and targets future
submission to Bioconductor.

You can install the package directly from GitHub using `remotes` or
`pak`:

``` r

# Install remotes if needed
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

# Install metaXpress directly from GitHub
remotes::install_github("hossainlab/metaXpress")
```

Alternatively, using `pak`:

``` r

# install.packages("pak")
pak::pkg_install("hossainlab/metaXpress")
```

Once installed, load the library:

``` r

library(metaXpress)
```

## Ingest Public Data from GEO

You can fetch public RNA-seq count matrices and curated metadata
directly from NCBI GEO using accession numbers:

``` r

# Fetch two bulk RNA-seq cohorts from GEO
studies <- mx_fetch_geo(c("GSE53697", "GSE95587"))
lapply(studies, show)
```

## Step 1 — QC Filtering

Each study is evaluated against a 10-point quality control standard
(replicate counts, depth, alignment rate, raw integer counts, etc.):

``` r

# Filter studies passing QC threshold (score >= 7)
studies <- mx_filter_studies(studies, qc_threshold = 7)
message("Studies passing QC: ", length(studies))
```

## Step 2 — Align Gene Spaces

Identify the common gene universe across all cohorts:

``` r

studies <- mx_align_genes(studies)
```

## Step 3 — Per-Study Differential Expression

Run negative binomial GLMs on each study independently using DESeq2,
edgeR, or limma-voom:

``` r

studies <- mx_de_all(studies, method = "DESeq2", formula = ~ condition)
mx_de_summary(studies)
```

## Step 4 — Meta-Analysis

Combine per-study effect sizes and standard errors using the
DerSimonian-Laird Random Effects Model:

``` r

de_results <- lapply(studies, function(s) s@de_result)
meta_result <- mx_meta(de_results, method = "random_effects")
meta_result
```

## Step 5 — Visualise

Generate publication-grade plots:

``` r

# Volcano plot of meta-analysis results
mx_volcano(meta_result, label_top = 10)
```

``` r

# Multi-study forest plot for the top gene
top_gene <- meta_result@meta_table$gene_id[
  which.min(meta_result@meta_table$meta_padj)]
mx_forest(top_gene, de_results, meta_result)
```

``` r

# Between-study heterogeneity distribution
mx_heterogeneity_plot(meta_result)
```

## Step 6 — Interactive Shiny Explorer

You can also explore your meta-analysis results interactively in your
browser:

``` r

# Launch local Shiny dashboard
metaXpress::mx_run_app()
```

Or test the online demo at
<https://hossainlab.shinyapps.io/metaXpress-demo/>.

## Step 7 — Export Results

``` r

mx_export(meta_result, format = "csv", output_dir = tempdir())
```

## Session Info

``` r

mx_session_info()$session_info
```
