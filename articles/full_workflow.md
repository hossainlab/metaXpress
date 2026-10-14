# Full workflow: from GEO accessions to report

## Overview

This vignette shows the full metaXpress pipeline from GEO accessions to
a report. Each step needs internet access or takes minutes on full data,
so the chunks are not evaluated when the package is built. For executed
analyses of real data, see:

- [`vignette("quickstart")`](https://mdjubayerhossain.com/metaXpress/articles/quickstart.md):
  a five-minute tour.
- [`vignette("case_study")`](https://mdjubayerhossain.com/metaXpress/articles/case_study.md):
  a complete Alzheimer’s disease meta-analysis.

The steps below reproduce that case study from scratch.

``` r

if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
remotes::install_github("hossainlab/metaXpress")
library(metaXpress)
```

------------------------------------------------------------------------

## 1. Data ingestion and quality control

### 1a. Uniformly processed counts from GEO

By default,
[`mx_fetch_geo()`](https://mdjubayerhossain.com/metaXpress/reference/mx_fetch_geo.md)
downloads the **NCBI GEO-generated raw counts**. NCBI aligns the SRA
reads of every eligible human and mouse series with HISAT2 and counts
them with featureCounts against RefSeq. All studies therefore share one
alignment, annotation and count definition, which removes a major source
of technical heterogeneity. Genes are Entrez IDs, and symbols and gene
types are added to `rowData`.

Always name the phenotype column that holds the case/control labels.
Automatic detection is only a fallback.

``` r

accessions <- c("GSE125583", "GSE104704", "GSE53697", "GSE173955")
studies <- mx_fetch_geo(
  accessions,
  condition_col = c(GSE125583 = "diagnosis:ch1",
                    GSE104704 = "study group:ch1",
                    GSE53697  = "disease status:ch1",
                    GSE173955 = "title"),
  cache_dir = "geo_cache")
lapply(studies, function(s) table(s$condition))
S4Vectors::metadata(studies[[1]])$count_source
```

Series without NCBI counts fall back to the author’s supplementary files
(`source = "supplementary"`). Local matrices load with
[`mx_load_local()`](https://mdjubayerhossain.com/metaXpress/reference/mx_load_local.md),
and any `SummarizedExperiment` (e.g. from recount3 or tximeta) converts
with
[`mx_from_se()`](https://mdjubayerhossain.com/metaXpress/reference/mx_from_se.md).

### 1b. Study selection is a scientific step

Before the analysis, record for every study:

- tissue or region;
- disease definition and control type;
- library type;
- **donor independence**. Two GEO series from the same lab or brain bank
  often share donors (GSE95587 and GSE125583 share 89 cases), and
  pooling them violates the independence assumption.

### 1c. Data-level QC

``` r

studies <- lapply(studies, mx_qc_study, min_genes = 15000)
studies[[1]]                          # "x / 10 passed (y failed, z not evaluable)"
S4Vectors::metadata(studies[[1]])$qc  # per-criterion status and details
studies <- mx_filter_studies(studies) # drop studies failing required criteria
```

Each criterion is recorded as pass, fail or not evaluable. Alignment
rate, for example, is usually missing from GEO. `qc_score` counts the
criteria passed and is never rescaled. By default,
[`mx_filter_studies()`](https://mdjubayerhossain.com/metaXpress/reference/mx_filter_studies.md)
removes studies that fail a criterion required for valid DE: at least 3
replicates per group, complete condition labels and raw counts.
Criterion 8 (exactly two groups) fails for series with extra groups
until they are subset in step 2.

------------------------------------------------------------------------

## 2. Preparing each study

Restrict each study to the two groups being compared, choose a gene
namespace and filter out low-count genes:

``` r

prep <- function(s, keep_levels) {
  s <- s[, s$condition %in% keep_levels]
  rd <- SummarizedExperiment::rowData(s)
  s <- s[!is.na(rd$symbol) & rd$gene_type == "protein-coding" &
           !duplicated(rd$symbol), ]
  rownames(s) <- SummarizedExperiment::rowData(s)$symbol
  s[rowSums(SummarizedExperiment::assay(s) >= 10) >=
      min(table(s$condition)), ]
}
contrasts <- list(GSE125583 = c("Alzheimer's disease", "control"),
                  GSE104704 = c("Aged, diseased", "Old"),
                  GSE53697  = c("advanced Alzheimer's Disease", "control"),
                  GSE173955 = c("Alzheimer's Disease",
                                "non-Alzheimer's Disease"))
studies$GSE173955$condition <- sub("_Biological.*$", "",
                                   studies$GSE173955$condition)
studies <- Map(prep, studies, contrasts)
```

For author-supplied count tables with mixed identifiers,
[`mx_reannotate()`](https://mdjubayerhossain.com/metaXpress/reference/mx_reannotate.md)
converts Ensembl, Entrez and symbols to one namespace.

**No batch correction before DE.** Each study is analysed separately, so
between-study batch effects never enter the per-study fold changes.
[`mx_normalize()`](https://mdjubayerhossain.com/metaXpress/reference/mx_normalize.md)
and
[`mx_remove_batch()`](https://mdjubayerhossain.com/metaXpress/reference/mx_remove_batch.md)
(which preserves `condition`) write new assays for pooled PCA, heatmaps
or co-expression. They never overwrite the raw counts.

------------------------------------------------------------------------

## 3. Per-study differential expression

Covariates enter the per-study design formula. The contrast is set
explicitly per study, so every fold change means case / control:

``` r

studies <- mx_de_all(
  studies,
  method   = "DESeq2",               # or "edgeR" (QL) / "limma-voom"
  formula  = list(~ sex + condition, ~ condition, ~ condition, ~ condition),
  contrast = contrasts,
  BPPARAM  = BiocParallel::MulticoreParam(4))
mx_de_summary(studies, padj_threshold = 0.05, lfc_threshold = 0)
de_results <- lapply(studies, function(s) s@de_result)
```

(`sex` must exist in `colData`, e.g.
`studies$GSE125583$sex <- studies$GSE125583$"Sex:ch1"`.)

------------------------------------------------------------------------

## 4. Pre-pooling checks

``` r

mx_study_concordance(de_results)   # negative r => flipped contrast
cov <- mx_missing_summary(de_results)
table(studies_measured = rowSums(cov))
```

[`mx_meta()`](https://mdjubayerhossain.com/metaXpress/reference/mx_meta.md)
analyses each gene over the studies that measured it (`min_studies`,
default 2). Prefer this to imputation:
[`mx_impute()`](https://mdjubayerhossain.com/metaXpress/reference/mx_impute.md)
fills fold changes from the other studies, which is circular.

------------------------------------------------------------------------

## 5. Meta-analysis

Choose the method a priori from the question being asked:

| Question | Method |
|----|----|
| Average effect across heterogeneous cohorts (recommended) | `"random_effects"`, `tau_method = "REML"` |
| Studies are technical replicates | `"fixed_effects"` |
| Only p-values or z-scores are available | `"inverse_normal"` (√n weights) or `"stouffer"` |
| Same-direction effect in at least one study | `"fisher"` (directional) |
| DE in all, or in at least r, studies | `"maxp"`, `"rop"` |

``` r

res <- mx_meta(de_results, method = "random_effects", tau_method = "REML")
mt  <- res@meta_table
head(mt[order(mt$meta_pvalue), c("gene_id", "meta_log2FC", "ci_lower",
                                 "ci_upper", "pi_lower", "pi_upper",
                                 "meta_padj", "i_squared", "tau2")])
```

The HKSJ adjustment (`hksj = TRUE`) gives better per-gene confidence
intervals with few studies. With fewer than about 5 studies, though, its
t-distribution makes genome-wide FDR control unreachable (see the case
study).

``` r

loo <- mx_sensitivity(de_results, method = "random_effects",
                      tau_method = "REML")
rob <- attr(loo, "robustness")
core <- rob$gene_id[rob$full_padj <= 0.05 & rob$frac_loo_sig == 1]
```

------------------------------------------------------------------------

## 6. Pathways

``` r

# GSEA on the signed meta z-score (fgsea + MSigDB)
res <- mx_pathway_meta(res, db = "Hallmarks", method = "GSEA")

# Per-study GSEA combined across studies (directional Stouffer)
ps <- mx_pathway_studies(de_results, db = "Reactome")
mx_pathway_heatmap(attr(ps, "per_study"), top_n = 25, value = "NES")
mx_pathway_dedup(ps, gene_sets = mx_gene_sets("Reactome"))
```

------------------------------------------------------------------------

## 7. Figures

``` r

mx_volcano(res, padj_threshold = 0.05, lfc_threshold = 0.5, label_top = 15)
mx_forest("GFAP", de_results, res)      # CI diamond + prediction interval
mx_heterogeneity_plot(res)
mx_upset(de_results, padj_threshold = 0.05, lfc_threshold = 0)
```

------------------------------------------------------------------------

## 8. Downstream biology

Module 9 turns the meta-signature into testable biology. Each function
takes a public prior-knowledge resource as a plain R object, so any
regulon, marker set or drug library can be used. The calls below show
how to load each resource.
[`vignette("downstream_biology")`](https://mdjubayerhossain.com/metaXpress/articles/downstream_biology.md)
runs them on bundled real data.

### 8a. Load Prior-Knowledge Resources

``` r

# CollecTRI TF-target regulons (OmniPath); OmniPath columns are read directly
collectri <- read.delim(
  "https://omnipathdb.org/interactions?datasets=collectri&genesymbols=yes&organisms=9606"
)

# MCP-counter cell-type markers: list or data.frame (gene, cell_type)
mcp <- read.delim(
  "https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt",
  check.names = FALSE, colClasses = "character"
)
markers <- data.frame(gene = mcp$ENTREZID, cell_type = mcp$`Cell population`)

# LINCS L1000 drug signatures from Enrichr, paired into up/down sets
lincs <- mx_read_signatures("LINCS_L1000_Chem_Pert_Consensus_Sigs.txt")
```

### 8b. Transcription Factor Activity

``` r

# Use id_col = "symbol" when gene_id holds Entrez/Ensembl IDs
tf <- mx_tf_activity(res, collectri, id_col = "symbol", min_size = 10)
head(tf)
mx_downstream_plot(tf, top = 30)
```

### 8c. Cell-Type Composition Shifts

``` r

comp <- mx_cell_composition(studies, markers)
comp
mx_downstream_plot(comp)
```

### 8d. Consensus Co-Expression Modules

``` r

net <- mx_coexpression_modules(studies, n_genes = 3000, consensus = "min",
                               cut_height = 0.3, meta_result = res)
net$module_stats      # pooled eigengene effect, hub genes, meta-DEG overlap
head(net$modules)     # gene-to-module assignment with kME
mx_downstream_plot(net)
```

### 8e. Drug Signature Reversal

``` r

drugs <- mx_drug_connectivity(res, lincs, id_col = "symbol")
head(drugs, 20)       # negative reversal_score = reverses the disease signature
mx_downstream_plot(drugs, top = 30)
```

------------------------------------------------------------------------

## 9. Report Generation

### Generate Reproducible Reports & Export

``` r

mx_report(res, de_results, output_dir = "results")    # HTML with methods text
mx_export(res, format = "csv", output_dir = "results")
```

------------------------------------------------------------------------

## Session Information

``` r

mx_session_info()$session_info
```
