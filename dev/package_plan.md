# metaXpress: Development Plan
> An end-to-end R package for bulk RNA-seq meta-analysis

---

## Overview

**Package name:** `metaXpress`  
**Target repository:** Bioconductor (primary), CRAN (secondary)  
**Language:** R (core), with optional Rcpp for performance-critical steps  
**License:** MIT  
**Goal:** Provide a unified, reproducible pipeline for integrating multiple bulk RNA-seq studies — from raw GEO data retrieval through meta-analysis statistics to pathway-level interpretation and reporting.

---

## Problem Statement

Existing tools address only fragments of the meta-analysis workflow:
- **DExMA** — meta-analysis statistics only, no data ingestion or normalization
- **GEDI** — harmonization only, no meta-analysis statistics
- **limma** — per-study DE only, not designed for cross-study integration
- **hCoCena** — co-expression focused, not DE-focused

**metaXpress** closes this gap with a single, opinionated, end-to-end pipeline.

---

## Target Users

- Computational biologists running cross-cohort transcriptomic studies
- Clinical researchers integrating public GEO datasets for biomarker discovery
- Bioinformaticians building reproducible meta-analysis workflows

---

## Package Architecture

```
metaXpress/
├── R/
│   ├── 01_ingest.R          # Module 1: Data ingestion & QC
│   ├── 02_harmonize.R       # Module 2: Normalization & harmonization
│   ├── 03_de.R              # Module 3: Per-study DE analysis
│   ├── 04_meta.R            # Module 4: Meta-analysis statistics
│   ├── 05_missing.R         # Module 5: Missing gene handling
│   ├── 06_pathway.R         # Module 6: Pathway meta-analysis
│   ├── 07_visualize.R       # Module 7: Visualization
│   ├── 08_report.R          # Module 8: Report generation
│   └── utils.R              # Shared utilities
├── tests/
│   └── testthat/            # Unit tests per module
├── vignettes/
│   ├── quickstart.Rmd
│   └── full_workflow.Rmd
├── data/                    # Example datasets
├── DESCRIPTION
├── NAMESPACE
└── README.md
```

---

## Module Specifications

### Module 1 — Data Ingestion & QC
**File:** `01_ingest.R`

**Purpose:** Fetch, parse, and quality-check studies from public repositories or user-supplied files.

**Key functions:**
| Function | Description |
|---|---|
| `mx_fetch_geo(accessions)` | Download count matrices + metadata from GEO via GEOquery |
| `mx_fetch_sra(srp_ids)` | Fetch raw data from SRA |
| `mx_load_local(paths, metadata)` | Load user-supplied count matrices |
| `mx_qc_study(study)` | Apply 10-point QC checklist (Heberle et al. 2025) |
| `mx_cluster_samples(study)` | Auto-cluster samples by metadata (sampleclusteR approach) |
| `mx_filter_studies(studies, qc_threshold)` | Remove studies below QC threshold |

**Inputs:** GEO accession IDs, SRA project IDs, or local file paths  
**Outputs:** Standardized `metaXpressStudy` S4 object list

**QC Criteria (10-point, from Heberle et al. 2025):**
1. Minimum sample size per group (n ≥ 3)
2. Sequencing depth (≥ 10M reads)
3. Alignment rate (≥ 70%)
4. rRNA contamination (< 10%)
5. Duplicate rate (< 50%)
6. Gene detection rate (≥ 15,000 genes)
7. Metadata completeness
8. Clearly defined case/control groups
9. No batch confounding with condition
10. Raw counts available (not pre-normalized)

---

### Module 2 — Normalization & Harmonization
**File:** `02_harmonize.R`

**Purpose:** Ensure all studies are on a comparable scale before meta-analysis.

**Key functions:**
| Function | Description |
|---|---|
| `mx_reannotate(study, org)` | Reannotate gene IDs to a common namespace (Ensembl/Symbol) |
| `mx_normalize(study, method)` | Normalize counts (TMM, DESeq2 VST, CPM) |
| `mx_correct_library_type(studies)` | Correct polyA vs rRNA-depleted bias (Bush et al. 2017) |
| `mx_remove_batch(studies, method)` | Batch effect removal (ComBat, limma::removeBatchEffect, Harmony) |
| `mx_align_genes(studies)` | Find common gene universe across all studies |

**Supported normalization methods:** TMM, DESeq2-VST, CPM, TPM, quantile  
**Supported batch methods:** ComBat, ComBat-seq, limma, Harmony

---

### Module 3 — Per-Study Differential Expression
**File:** `03_de.R`

**Purpose:** Run DE analysis independently on each study using a consistent interface.

**Key functions:**
| Function | Description |
|---|---|
| `mx_de(study, method, formula)` | Run DE analysis on a single study |
| `mx_de_all(studies, method, formula)` | Run DE on all studies in parallel |
| `mx_de_summary(de_results)` | Summarize per-study DE results |

**Supported DE methods:** DESeq2, edgeR, limma-voom  
**Inputs:** tximport-compatible (Salmon/kallisto) or raw count matrices  
**Outputs:** Standardized `metaXpressDEResult` object with log2FC, p-value, padj per gene per study

---

### Module 4 — Meta-Analysis Statistics
**File:** `04_meta.R`

**Purpose:** Integrate per-study DE results into a single meta-analysis result.

**Key functions:**
| Function | Description |
|---|---|
| `mx_meta(de_results, method)` | Run meta-analysis with chosen method |
| `mx_heterogeneity(de_results)` | Compute I² and Q-statistic per gene |
| `mx_forest(gene, de_results)` | Generate forest plot for a single gene |
| `mx_sensitivity(de_results)` | Leave-one-out sensitivity analysis |

**Supported meta-analysis methods:**

| Method | Reference | Best For |
|---|---|---|
| Fisher's combined p-value | Rau et al. 2013 | Many studies, independent p-values |
| Stouffer's Z-score | Rau et al. 2013 | Weighted by sample size |
| Inverse-normal (fused) | Prasad & Li 2021 | Small-n studies |
| Fixed effects (effect size) | Keel & Lindholm-Perry 2022 | Low heterogeneity |
| Random effects (DerSimonian-Laird) | Keel & Lindholm-Perry 2022 | High heterogeneity |
| Adaptive weighting (AWmeta) | Hu et al. 2025 | Mixed heterogeneity |

**Outputs:** Meta-analysis result table with meta-p, meta-FDR, meta-effect size, I², direction consistency score

---

### Module 5 — Missing Gene Handling
**File:** `05_missing.R`

**Purpose:** Handle genes not present in all studies (partial overlap is common in real-world meta-analyses).

**Key functions:**
| Function | Description |
|---|---|
| `mx_missing_summary(studies)` | Report gene coverage across studies |
| `mx_impute(de_results, method)` | Impute missing gene statistics |
| `mx_filter_coverage(de_results, min_studies)` | Filter genes by minimum study coverage |

**Imputation strategies (from DExMA):**
- Exclude genes below coverage threshold
- Mean imputation from available studies
- KNN-based imputation
- Study-weighted imputation

---

### Module 6 — Pathway Meta-Analysis
**File:** `06_pathway.R`

**Purpose:** Perform pathway enrichment at the meta-analysis level, not just per-study.

**Key functions:**
| Function | Description |
|---|---|
| `mx_pathway_meta(meta_results, db)` | Meta-analytic ORA/GSEA across studies |
| `mx_pathway_consensus(pathway_results)` | Find consensus pathways across studies |
| `mx_pathway_dedup(pathway_results)` | Remove redundant pathways (CPI approach) |
| `mx_pathway_heatmap(pathway_results)` | Cross-study pathway heatmap |

**Supported databases:** KEGG, Reactome, GO, MSigDB Hallmarks, WikiPathways  
**Methods:** ORA (Fisher's exact), GSEA, meta-GSEA

---

### Module 7 — Visualization
**File:** `07_visualize.R`

**Purpose:** Publication-ready figures for meta-analysis results.

**Key functions:**
| Function | Description |
|---|---|
| `mx_volcano(meta_results)` | Meta-analysis volcano plot |
| `mx_forest(gene, de_results)` | Forest plot with I² for a gene |
| `mx_heatmap(meta_results, top_n)` | Top DEG heatmap across studies |
| `mx_upset(de_results)` | UpSet plot of DEG overlap across studies |
| `mx_study_overview(studies)` | Study characteristics summary plot |
| `mx_heterogeneity_plot(meta_results)` | I² distribution across genes |

---

### Module 8 — Report Generation
**File:** `08_report.R`

**Purpose:** Generate a reproducible HTML/PDF report of the full analysis.

**Key functions:**
| Function | Description |
|---|---|
| `mx_report(pipeline_result, format)` | Render full analysis report (HTML/PDF) |
| `mx_export(meta_results, format)` | Export results (CSV, Excel, RDS) |
| `mx_session_info()` | Capture full session info for reproducibility |

---

## Data Structures

### `metaXpressStudy` (S4 class)
```r
setClass("metaXpressStudy", representation(
  counts     = "matrix",       # Raw count matrix (genes x samples)
  metadata   = "data.frame",   # Sample metadata
  accession  = "character",    # GEO/SRA accession
  organism   = "character",    # e.g., "Homo sapiens"
  qc_score   = "numeric",      # 0-10 QC score
  de_result  = "data.frame"    # Per-study DE result (populated by Module 3)
))
```

### `metaXpressResult` (S4 class)
```r
setClass("metaXpressResult", representation(
  meta_table    = "data.frame",  # Gene-level meta-analysis results
  method        = "character",   # Meta-analysis method used
  n_studies     = "integer",     # Number of studies integrated
  heterogeneity = "data.frame",  # I², Q-stat per gene
  pathway_result = "data.frame"  # Pathway meta-analysis results
))
```

---

## Typical Workflow

```r
library(metaXpress)

# 1. Fetch and QC studies
studies <- mx_fetch_geo(c("GSE12345", "GSE67890", "GSE11111"))
studies <- mx_filter_studies(studies, qc_threshold = 7)

# 2. Harmonize
studies <- mx_reannotate(studies, org = "Homo sapiens")
studies <- mx_correct_library_type(studies)
studies <- mx_remove_batch(studies, method = "ComBat-seq")

# 3. Per-study DE
de_results <- mx_de_all(studies, method = "DESeq2",
                         formula = ~ condition)

# 4. Meta-analysis
meta_result <- mx_meta(de_results, method = "random_effects")

# 5. Pathway meta-analysis
pathway_result <- mx_pathway_meta(meta_result, db = "Hallmarks")

# 6. Visualize
mx_volcano(meta_result)
mx_forest("TP53", de_results)
mx_pathway_heatmap(pathway_result)

# 7. Report
mx_report(meta_result, format = "html")
```

---

## Development Roadmap

### Phase 1 — Foundation (Months 1–2)
- [ ] Set up package skeleton (`usethis::create_package`)
- [ ] Implement `metaXpressStudy` and `metaXpressResult` S4 classes
- [ ] Module 1: GEO ingestion + QC
- [ ] Module 3: DESeq2 per-study DE wrapper
- [ ] Unit tests for Modules 1 & 3
- [ ] Quickstart vignette

### Phase 2 — Core Meta-Analysis (Months 3–4)
- [ ] Module 2: Normalization & harmonization
- [ ] Module 4: P-value combination methods (Fisher, Stouffer, inverse-normal)
- [ ] Module 5: Missing gene handling
- [ ] Module 4: Effect size models (fixed + random effects)
- [ ] Unit tests for Modules 2, 4 & 5
- [ ] Benchmark against DExMA on Alzheimer's dataset (Heberle et al. 2025)

### Phase 3 — Advanced Features (Months 5–6)
- [ ] Module 4: AWmeta adaptive weighting
- [ ] Module 6: Pathway meta-analysis
- [ ] Module 7: Full visualization suite
- [ ] Module 8: HTML/PDF report generation
- [ ] Full workflow vignette

### Phase 4 — Polish & Release (Month 7)
- [ ] Bioconductor submission checklist
- [ ] pkgdown documentation site
- [ ] Performance optimization (parallelization via BiocParallel)
- [ ] CRAN/Bioconductor submission

---

## Dependencies

### Required
| Package | Purpose |
|---|---|
| `GEOquery` | GEO data fetching |
| `DESeq2` | Per-study DE analysis |
| `edgeR` | Per-study DE analysis |
| `limma` | Per-study DE + batch correction |
| `sva` | ComBat batch correction |
| `tximport` | Transcript-level input support |
| `clusterProfiler` | Pathway enrichment |
| `msigdbr` | MSigDB gene sets |
| `ggplot2` | Visualization |
| `ComplexHeatmap` | Heatmaps |
| `BiocParallel` | Parallelization |

### Suggested
| Package | Purpose |
|---|---|
| `harmony` | Alternative batch correction |
| `AnnotationDbi` | Gene ID reannotation |
| `org.Hs.eg.db` | Human gene annotation |
| `rmarkdown` | Report generation |
| `openxlsx` | Excel export |

---

## Key References

| # | Citation | Module |
|---|---|---|
| 1 | Rau et al. (2013) BMC Bioinformatics | Module 4 |
| 2 | Prasad & Li (2021) BMC Bioinformatics | Module 4 |
| 3 | Hu et al. (2025) bioRxiv | Module 4 |
| 4 | Keel & Lindholm-Perry (2022) Front. Genetics | Module 4 |
| 5 | Villatoro-García et al. (2022) Mathematics | Module 5 |
| 6 | Stokholm et al. (2024) Current Protocols | Module 2 |
| 7 | Bush et al. (2017) BMC Bioinformatics | Module 2 |
| 8 | Ritchie et al. (2015) Nucleic Acids Research | Module 3 |
| 9 | Soneson et al. (2015) F1000Research | Module 3 |
| 10 | Zeng et al. (2018) Genes | Module 6 |
| 11 | Heberle et al. (2025) Alzheimer's & Dementia | Module 1 |
| 12 | Coke et al. (2025) bioRxiv | Module 1 |
