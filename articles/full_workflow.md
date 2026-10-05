# metaXpress Full Workflow: GEO to Report

## Overview

This vignette demonstrates the complete metaXpress workflow, from
fetching studies from GEO through pathway enrichment and report
generation. All `mx_fetch_*()` calls fetch public data directly from
NCBI GEO. Users can provide accession IDs or local count matrices via
[`mx_load_local()`](https://hossainlab.github.io/metaXpress/reference/mx_load_local.md).

``` r

library(metaXpress)
```

## 1. Data Ingestion

### 1a. Fetch from GEO

``` r

# Fetch three Alzheimer's disease RNA-seq studies
accessions <- c("GSE53697", "GSE95587", "GSE118553")
studies    <- mx_fetch_geo(accessions, count_type = "raw",
                            cache_dir = "geo_cache/")
```

### 1b. Inspect QC Results

``` r

qc_scores <- vapply(studies, function(s) s@qc_score, numeric(1))
print(qc_scores)

# View detailed QC breakdown for the first study
attr(studies[[1]], "qc_details")
```

### 1c. Filter Low-Quality Studies

``` r

studies <- mx_filter_studies(studies, qc_threshold = 7)
```

## 2. Harmonization

### 2a. Reannotate Gene IDs

``` r

studies <- mx_reannotate(studies, org = "Homo sapiens",
                          target_id = "SYMBOL")
```

### 2b. Correct Library Type Bias

If studies mix polyA-selected and rRNA-depleted libraries:

``` r

studies <- mx_correct_library_type(studies)
```

### 2c. Remove Batch Effects

``` r

studies <- mx_remove_batch(studies, method = "ComBat-seq")
```

### 2d. Align to Common Genes

``` r

studies <- mx_align_genes(studies)
```

## 3. Per-Study Differential Expression

``` r

studies <- mx_de_all(studies, method = "DESeq2",
                      formula = ~ condition,
                      BPPARAM = BiocParallel::MulticoreParam(4))
```

``` r

mx_de_summary(studies, padj_threshold = 0.05, lfc_threshold = 1)
```

## 4. Handle Missing Genes

``` r

de_results <- lapply(studies, function(s) s@de_result)

# Check coverage across studies
cov_mat <- mx_missing_summary(de_results)
hist(attr(cov_mat, "coverage_pct"),
     main = "Gene coverage across studies",
     xlab = "% studies with gene detected")

# Keep genes in at least 2 of 3 studies
de_results <- mx_filter_coverage(de_results, min_studies = 2)
```

## 5. Meta-Analysis

### 5a. Choose the Right Method

Use the decision guide: - 2 studies → `"fisher"` or `"stouffer"` - 3+
studies, low heterogeneity → `"fixed_effects"` - 3+ studies, high
heterogeneity → `"random_effects"` (default) - Mixed heterogeneity →
`"awmeta"`

``` r

meta_result <- mx_meta(de_results, method = "random_effects",
                        min_studies = 2)
meta_result
```

### 5b. Heterogeneity Assessment

``` r

het <- mx_heterogeneity(de_results)
summary(het$I_sq)

# Genes with high heterogeneity
high_het <- het[!is.na(het$I_sq) & het$I_sq > 75, ]
nrow(high_het)
```

### 5c. Sensitivity Analysis

``` r

loo <- mx_sensitivity(de_results, method = "random_effects")
# Compare significant gene counts across LOO runs
vapply(loo, function(r) sum(r@meta_table$meta_padj <= 0.05, na.rm = TRUE),
        integer(1))
```

## 6. Pathway Meta-Analysis

``` r

meta_result <- mx_pathway_meta(meta_result, db = "Hallmarks",
                                method = "ORA")
head(meta_result@pathway_result[
  order(meta_result@pathway_result$padj), ], 10)
```

## 7. Visualization

``` r

mx_volcano(meta_result, padj_threshold = 0.05, lfc_threshold = 1,
            label_top = 15)
```

``` r

top_gene <- meta_result@meta_table$gene_id[
  which.min(meta_result@meta_table$meta_padj)]
mx_forest(top_gene, de_results, meta_result)
```

``` r

mx_heatmap(meta_result, studies, top_n = 30)
```

``` r

mx_heterogeneity_plot(meta_result)
```

## 8. Report & Export

``` r

mx_report(meta_result, studies, de_results,
           format     = "html",
           output_dir = "results/")
```

``` r

mx_export(meta_result, format = "excel", output_dir = "results/")
mx_export(meta_result, format = "csv",   output_dir = "results/")
```

## Session Information

``` r

mx_session_info()$session_info
```
