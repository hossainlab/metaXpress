<div align="center">

# metaXpress

**An end-to-end R/Bioconductor package for bulk RNA-seq meta-analysis**

[![R CMD Check](https://github.com/hossainlab/metaXpress/actions/workflows/check.yaml/badge.svg)](https://github.com/hossainlab/metaXpress/actions/workflows/check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-submission-blueviolet?logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCA1MCA1MCI+PHBhdGggZmlsbD0iI2ZmZiIgZD0iTTI1IDVDMTMuOSA1IDUgMTMuOSA1IDI1czguOSAyMCAyMCAyMCAyMC04LjkgMjAtMjBTMzYuMSA1IDI1IDV6bTAgMzZjLTguOCAwLTE2LTcuMi0xNi0xNlMxNi4yIDkgMjUgOXMxNiA3LjIgMTYgMTYtNy4yIDE2LTE2IDE2eiIvPjwvc3ZnPg==)](https://bioconductor.org/)
[![R version](https://img.shields.io/badge/R-%E2%89%A54.3.0-276DC3?logo=r)](https://www.r-project.org/)
[![GitHub issues](https://img.shields.io/github/issues/hossainlab/metaXpress?color=red)](https://github.com/hossainlab/metaXpress/issues)
[![GitHub last commit](https://img.shields.io/github/last-commit/hossainlab/metaXpress)](https://github.com/hossainlab/metaXpress/commits/main)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](https://github.com/hossainlab/metaXpress/pulls)

</div>

---

## Overview

Most transcriptomic meta-analysis tools address only one part of the
workflow — existing packages handle data harmonization **or** DE analysis **or**
meta-statistics, but not the full pipeline. **metaXpress** closes this gap.

It provides a single, opinionated, end-to-end pipeline:

```
GEO / SRA / Local files
        │
        ▼
┌───────────────────┐
│  Module 1: Ingest │  mx_fetch_geo · mx_load_local · mx_qc_study
└────────┬──────────┘
         │
         ▼
┌──────────────────────┐
│ Module 2: Harmonize  │  mx_reannotate · mx_normalize · mx_remove_batch
└────────┬─────────────┘
         │
         ▼
┌─────────────────────┐
│  Module 3: Per-Study │  mx_de · mx_de_all  (DESeq2 / edgeR / limma-voom)
│  Differential Expr. │
└────────┬────────────┘
         │
         ▼
┌──────────────────────┐
│ Module 5: Missing    │  mx_impute · mx_filter_coverage
│ Gene Handling        │
└────────┬─────────────┘
         │
         ▼
┌──────────────────────┐
│  Module 4: Meta-     │  mx_meta  (6 statistical methods)
│  Analysis Statistics │
└────────┬─────────────┘
         │
         ├─────────────────────────┐
         ▼                         ▼
┌──────────────────┐    ┌────────────────────┐
│  Module 6:       │    │  Module 7:         │
│  Pathway         │    │  Visualization     │
│  Meta-Analysis   │    │  (volcano · forest │
│                  │    │   heatmap · UpSet) │
└────────┬─────────┘    └────────────────────┘
         │
         ▼
┌──────────────────────┐
│  Module 8: Report    │  mx_report (HTML / PDF) · mx_export
└──────────────────────┘
```

---

## Why metaXpress?

| Feature | metaXpress | DExMA | GEDI | limma |
|---|:---:|:---:|:---:|:---:|
| GEO data ingestion | ✅ | ❌ | ❌ | ❌ |
| 10-point study QC | ✅ | ❌ | ❌ | ❌ |
| Cross-platform harmonization | ✅ | ❌ | ✅ | ❌ |
| Library-type correction | ✅ | ❌ | ✅ | ❌ |
| Per-study DE (3 engines) | ✅ | ❌ | ❌ | ✅ |
| 6 meta-analysis methods | ✅ | ✅ | ❌ | ❌ |
| Missing gene imputation | ✅ | ✅ | ❌ | ❌ |
| Pathway meta-analysis | ✅ | ❌ | ❌ | ❌ |
| Reproducible report | ✅ | ❌ | ❌ | ❌ |

---

## Installation

```r
# Bioconductor (once released)
BiocManager::install("metaXpress")

# Development version from GitHub
remotes::install_github("hossainlab/metaXpress")
```

---

## Quick Start

```r
library(metaXpress)

# ── 1. Fetch and QC studies from GEO ──────────────────────────────────────
studies <- mx_fetch_geo(c("GSE53697", "GSE95587", "GSE118553"))
studies <- mx_filter_studies(studies, qc_threshold = 7)

# ── 2. Harmonize ──────────────────────────────────────────────────────────
studies <- mx_reannotate(studies, org = "Homo sapiens", target_id = "SYMBOL")
studies <- mx_correct_library_type(studies)
studies <- mx_remove_batch(studies, method = "ComBat-seq")
studies <- mx_align_genes(studies)

# ── 3. Per-study DE ───────────────────────────────────────────────────────
studies <- mx_de_all(studies, method = "DESeq2", formula = ~ condition)
mx_de_summary(studies)

# ── 4. Handle missing genes ───────────────────────────────────────────────
de_results <- lapply(studies, function(s) s@de_result)
de_results <- mx_filter_coverage(de_results, min_studies = 2)

# ── 5. Meta-analysis ──────────────────────────────────────────────────────
meta_result <- mx_meta(de_results, method = "random_effects")
meta_result   # prints: method, n_studies, n_genes, n_sig

# ── 6. Pathway enrichment ─────────────────────────────────────────────────
meta_result <- mx_pathway_meta(meta_result, db = "Hallmarks")

# ── 7. Visualize ──────────────────────────────────────────────────────────
mx_volcano(meta_result, label_top = 15)
mx_forest("TP53", de_results, meta_result)
mx_heterogeneity_plot(meta_result)

# ── 8. Export ─────────────────────────────────────────────────────────────
mx_report(meta_result, studies, de_results, format = "html")
mx_export(meta_result, format = "excel")
```

---

## Meta-Analysis Methods

metaXpress implements six methods selectable via the `method` argument of
`mx_meta()`:

| Method | `method =` | Effect Size | Handles Heterogeneity | Min Studies | Reference |
|---|---|:---:|:---:|:---:|---|
| Fisher's combined p-value | `"fisher"` | ❌ | ❌ | 2 | Rau et al. 2013 |
| Stouffer's Z-score | `"stouffer"` | ❌ | ❌ | 2 | Rau et al. 2013 |
| Fused inverse-normal | `"inverse_normal"` | ❌ | Partial | 2 | Prasad & Li 2021 |
| Fixed effects | `"fixed_effects"` | ✅ | ❌ | 2 | Keel & Lindholm-Perry 2022 |
| **Random effects** *(default)* | `"random_effects"` | ✅ | ✅ | 3 | Keel & Lindholm-Perry 2022 |
| Adaptive weighting (AWmeta) | `"awmeta"` | ✅ | ✅ | 3 | Hu et al. 2025 |

**Choosing a method:**
- **2 studies** → use `"fisher"` or `"stouffer"`
- **3+ studies, I² < 25%** → use `"fixed_effects"`
- **3+ studies, I² ≥ 25%** → use `"random_effects"` (default)
- **Mixed/unknown heterogeneity** → use `"awmeta"`

---

## S4 Data Structures

```
mx_fetch_geo() ──► metaXpressStudy ──► mx_de() ──► mx_meta() ──► metaXpressResult
                   ├── counts                                      ├── meta_table
                   ├── metadata                                    ├── method
                   ├── accession                                   ├── n_studies
                   ├── organism                                    ├── heterogeneity
                   ├── qc_score (0–10)                            └── pathway_result
                   └── de_result
```

---

## Study QC Scoring

`mx_qc_study()` scores each study against 10 criteria
(Heberle et al. 2025), awarding 1 point per criterion:

| # | Criterion | Threshold |
|---|---|---|
| 1 | Minimum sample size | ≥ 3 replicates per group |
| 2 | Sequencing depth | Median ≥ 10M reads |
| 3 | Alignment rate | Mean ≥ 70% |
| 4 | rRNA contamination | < 10% |
| 5 | Duplicate rate | < 50% |
| 6 | Gene detection rate | ≥ 15,000 genes in ≥ 50% of samples |
| 7 | Metadata completeness | `condition` + `sample_id` present |
| 8 | Clear case/control | Exactly 2 condition levels |
| 9 | No batch confounding | Batch not perfectly correlated with condition |
| 10 | Raw counts | Integer counts (not FPKM/TPM) |

Studies scoring below `qc_threshold = 7` are removed by `mx_filter_studies()`.

---

## Function Reference

<details>
<summary><strong>Module 1 — Data Ingestion & QC</strong></summary>

| Function | Description |
|---|---|
| `mx_fetch_geo(accessions)` | Download count matrices + metadata from GEO |
| `mx_fetch_sra(srp_ids)` | Fetch from SRA *(planned)* |
| `mx_load_local(count_paths, metadata_paths)` | Load user-supplied files |
| `mx_qc_study(study)` | Apply 10-point QC checklist |
| `mx_cluster_samples(study)` | Auto-cluster samples from metadata *(planned)* |
| `mx_filter_studies(studies, qc_threshold)` | Remove studies below QC threshold |

</details>

<details>
<summary><strong>Module 2 — Normalization & Harmonization</strong></summary>

| Function | Description |
|---|---|
| `mx_reannotate(studies, org, target_id)` | Standardize gene ID namespace |
| `mx_normalize(study, method)` | TMM / VST / CPM / TPM / quantile |
| `mx_correct_library_type(studies)` | polyA vs rRNA-depleted correction |
| `mx_remove_batch(studies, method)` | ComBat-seq / ComBat / limma / harmony |
| `mx_align_genes(studies)` | Restrict to common gene universe |

</details>

<details>
<summary><strong>Module 3 — Per-Study DE</strong></summary>

| Function | Description |
|---|---|
| `mx_de(study, method, formula)` | DESeq2 / edgeR / limma-voom on one study |
| `mx_de_all(studies, method, BPPARAM)` | Parallel DE across all studies |
| `mx_de_summary(studies)` | n DEGs per study summary table |

</details>

<details>
<summary><strong>Module 4 — Meta-Analysis Statistics</strong></summary>

| Function | Description |
|---|---|
| `mx_meta(de_results, method)` | Run meta-analysis (6 methods) |
| `mx_heterogeneity(de_results)` | I², Q-stat, τ² per gene |
| `mx_sensitivity(de_results)` | Leave-one-out sensitivity analysis |

</details>

<details>
<summary><strong>Module 5 — Missing Gene Handling</strong></summary>

| Function | Description |
|---|---|
| `mx_missing_summary(de_results)` | Gene × study coverage matrix |
| `mx_impute(de_results, method)` | exclude / mean / KNN / weighted |
| `mx_filter_coverage(de_results, min_studies)` | Filter by study coverage |

</details>

<details>
<summary><strong>Module 6 — Pathway Meta-Analysis</strong></summary>

| Function | Description |
|---|---|
| `mx_pathway_meta(meta_result, db)` | ORA / GSEA with Hallmarks, KEGG, Reactome, GO |
| `mx_pathway_consensus(pathway_results)` | Pathways significant across majority of studies |
| `mx_pathway_dedup(pathway_results)` | Remove redundant pathways |
| `mx_pathway_heatmap(pathway_results)` | Cross-study pathway heatmap |

</details>

<details>
<summary><strong>Module 7 — Visualization</strong></summary>

| Function | Description |
|---|---|
| `mx_volcano(meta_result)` | Meta-analysis volcano plot |
| `mx_forest(gene, de_results)` | Forest plot with I² for a single gene |
| `mx_heatmap(meta_result, studies)` | Top DEG heatmap across studies |
| `mx_upset(de_results)` | UpSet plot of DEG overlap |
| `mx_study_overview(studies)` | QC metrics summary plot |
| `mx_heterogeneity_plot(meta_result)` | I² distribution histogram |

</details>

<details>
<summary><strong>Module 8 — Report & Export</strong></summary>

| Function | Description |
|---|---|
| `mx_report(meta_result, ..., format)` | Render HTML / PDF report |
| `mx_export(meta_result, format)` | Export CSV / Excel / RDS |
| `mx_session_info()` | Capture session info for reproducibility |

</details>

---

## Development

```r
# Load all functions without installing (standard dev workflow)
devtools::load_all()       # Ctrl+Shift+L in RStudio

# Run all tests
devtools::test()

# Run a single test file
devtools::test_file("tests/testthat/test-meta.R")

# Full R CMD check
devtools::check()

# Bioconductor-specific checks
BiocCheck::BiocCheck(".")
```

---

## Roadmap

- [x] Package scaffold & S4 classes
- [x] Module 1: GEO ingestion + 10-point QC
- [x] Module 2: Normalization & batch correction
- [x] Module 3: DESeq2 / edgeR / limma-voom wrappers
- [x] Module 4: All 6 meta-analysis methods
- [x] Module 5: Missing gene handling
- [x] Module 6: Pathway meta-analysis (ORA + GSEA)
- [x] Module 7: Visualization suite
- [x] Module 8: Report generation & export
- [ ] Example dataset (`data/metaXpress_example.rda`)
- [ ] pkgdown documentation site
- [ ] Bioconductor submission

---

## Key References

| Method | Reference |
|---|---|
| Fisher / Stouffer meta-analysis | Rau, Marot & Jaffrézic (2013) *BMC Bioinformatics* |
| Fused inverse-normal | Prasad & Li (2021) *BMC Bioinformatics* |
| Random effects (DerSimonian-Laird) | Keel & Lindholm-Perry (2022) *Front. Genetics* |
| AWmeta adaptive weighting | Hu et al. (2025) *bioRxiv* |
| Missing gene imputation | Villatoro-García et al. (2022) *Mathematics* (DExMA) |
| Library-type correction | Bush et al. (2017) *BMC Bioinformatics* |
| Study QC criteria | Heberle et al. (2025) *Alzheimer's & Dementia* |
| Sample clustering from GEO | Coke, Niranjan & Ewing (2025) *bioRxiv* |

---

## License

MIT © [DeepBio Limited](https://github.com/hossainlab)
