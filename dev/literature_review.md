# Bulk RNA-seq Meta-Analysis: Literature Review
> For developing an end-to-end R package for transcriptomic meta-analysis

---

## Table of Contents
1. [Foundational Meta-Analysis Methods](#1-foundational-meta-analysis-methods)
2. [Existing R Packages to Study](#2-existing-r-packages-to-study-and-differentiate-from)
3. [Data Harmonization & Batch Correction](#3-data-harmonization--batch-correction)
4. [Pathway-Level Meta-Analysis](#4-pathway-level-meta-analysis)
5. [Applied Examples & Validation Use Cases](#5-applied-examples--validation-use-cases)
6. [Recommended Package Architecture](#6-recommended-package-architecture)
7. [Full Reference List](#7-full-reference-list)

---

## 1. Foundational Meta-Analysis Methods

### P-value Combination Approaches

#### Rau et al. (2013) — *Differential meta-analysis of RNA-seq data from multiple studies*
- **Journal:** BMC Bioinformatics | **Citations:** 144
- **DOI:** https://doi.org/10.1186/1471-2105-15-91
- **Summary:** The seminal paper on differential meta-analysis of RNA-seq across studies. Benchmarks Fisher, Stouffer, and other p-value combination methods.
- **Key Finding:** P-value combination outperforms fixed-effect models when inter-study variability is moderate-to-high.
- **Relevance:** Essential reading — defines the statistical foundation for your meta-analysis module.

---

#### Prasad & Li (2021) — *Fused inverse-normal method for integrated differential expression analysis of RNA-seq data*
- **Journal:** BMC Bioinformatics
- **DOI:** https://doi.org/10.1186/s12859-022-04859-9
- **Summary:** Proposes the *fused inverse-normal* method, a hybrid p-value combination approach tailored for small-n RNA-seq studies.
- **Key Finding:** Effective at identifying DEGs in glioblastoma; performs well under small sample sizes common in RNA-seq experiments.
- **Relevance:** Implement as a method option in your p-value combination module.

---

#### Hu et al. (2025) — *AWmeta empowers adaptively-weighted transcriptomic meta-analysis*
- **Journal:** bioRxiv (preprint)
- **DOI:** https://doi.org/10.1101/2025.05.06.650408
- **Summary:** Novel adaptively-weighted meta-analysis method that combines p-value and effect size approaches, addressing limitations of both paradigms.
- **Relevance:** Most recent advance — worth implementing as a premium method option. Reduces experimental costs and improves reliability.

---

### Effect Size / Model-Based Approaches

#### Keel & Lindholm-Perry (2022) — *Recent developments and future directions in meta-analysis of differential gene expression in livestock RNA-Seq*
- **Journal:** Frontiers in Genetics | **Citations:** 13
- **DOI:** https://doi.org/10.3389/fgene.2022.983043
- **Summary:** Comprehensive review of meta-analysis strategies for RNA-seq including fixed effects, random effects, and p-value combination methods.
- **Key Finding:** Meta-analyses increase sample size and statistical power; addresses unique challenges in cross-study integration.
- **Relevance:** Excellent methods comparison — use to justify your design choices and method selection logic.

---

## 2. Existing R Packages to Study and Differentiate From

| Package | Paper | What It Does | Gap Your Package Can Fill |
|---|---|---|---|
| **DExMA** | Villatoro-García et al. (2022) | Gene expression meta-analysis with missing genes | Limited batch correction; no end-to-end pipeline |
| **GEDI** | Stokholm et al. (2024) | Multi-platform integration (microarray + RNA-seq) | No meta-analysis statistics |
| **hCoCena** | Oestreich et al. (2022) | Co-expression meta-analysis across datasets | Co-expression only, not DE-focused |
| **limma** | Ritchie et al. (2015) | DE analysis, some multi-study support | Not designed for cross-study meta-analysis |
| **sampleclusteR** | Coke et al. (2025) | Automated sample clustering from GEO metadata | Upstream only; no downstream analysis |

---

### DExMA (2022) — *An R Package for Performing Gene Expression Meta-Analysis with Missing Genes*
- **Journal:** Mathematics | **Citations:** 9
- **DOI:** https://doi.org/10.3390/math10183376
- **Summary:** Provides meta-analysis functions for transcriptomic data with explicit handling of missing genes across studies.
- **Relevance:** Closest existing tool to your goal — study its API design and identify gaps (no GEO integration, limited normalization, no pathway meta-analysis).

---

### GEDI (2024) — *An R Package for Integration of Transcriptomic Data from Multiple Platforms*
- **Journal:** Current Protocols
- **DOI:** https://doi.org/10.1002/cpz1.70046
- **Summary:** Automatically reannotates and removes batch effects when integrating microarray and NGS transcriptomic data.
- **Relevance:** Strong batch correction module — study its reannotation approach for cross-platform gene ID harmonization.

---

### hCoCena (2022) — *Horizontal integration and analysis of transcriptomics datasets*
- **Journal:** Bioinformatics | **Citations:** 12
- **DOI:** https://doi.org/10.1093/bioinformatics/btac589
- **Summary:** User-friendly tool for analyzing multiple transcriptomic datasets via gene co-expression analysis across conditions.
- **Relevance:** Useful reference for multi-dataset visualization design; co-expression as an optional module.

---

### limma (2015) — *limma powers differential expression analyses for RNA-sequencing and microarray studies*
- **Journal:** Nucleic Acids Research | **Citations:** 31,041
- **DOI:** https://doi.org/10.1093/nar/gkv007
- **Summary:** The gold-standard R/Bioconductor package for DE analysis. Handles complex experimental designs and small sample sizes via empirical Bayes.
- **Relevance:** Use as the per-study DE engine within your pipeline (alongside DESeq2/edgeR wrappers).

---

### sampleclusteR (2025) — *A lightweight R package for automated clustering of transcriptomics samples using metadata*
- **Journal:** bioRxiv (preprint)
- **DOI:** https://doi.org/10.1101/2025.04.10.648129
- **Summary:** Automates clustering of GEO/ArrayExpress samples based on metadata, enabling large-scale meta-analysis setup.
- **Relevance:** Directly useful for your data ingestion module — automates the tedious sample grouping step.

---

## 3. Data Harmonization & Batch Correction

#### Bush et al. (2017) — *Integration of quantitated expression estimates from polyA-selected and rRNA-depleted RNA-seq libraries*
- **Journal:** BMC Bioinformatics | **Citations:** 33
- **DOI:** https://doi.org/10.1186/s12859-017-1714-9
- **Summary:** Addresses integration of polyA-selected vs. rRNA-depleted libraries — a real-world challenge when pulling GEO datasets. Proposes ratio-based correction.
- **Key Finding:** Reference transcriptome filtering + ratio-based correction creates equivalent expression profiles across library types.
- **Relevance:** Critical for your harmonization module — GEO datasets frequently mix library preparation methods.

---

#### Panahi et al. (2025) — *Unifying RNA-seq data using meta-analysis; bioinformatics frameworks and application for plant genomics*
- **Journal:** Current Plant Biology
- **DOI:** https://doi.org/10.1016/j.cpb.2025.100523
- **Summary:** Recent review of RNA-seq meta-analysis frameworks covering cross-study variability, platform differences, and normalization strategies.
- **Relevance:** Good for framing your package's scope and normalization design, even though plant-focused — the methodological challenges are universal.

---

#### Soneson, Love & Robinson (2015) — *Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences*
- **Journal:** F1000Research | **Citations:** 3,088
- **DOI:** https://doi.org/10.12688/f1000research.7563.2
- **Summary:** Demonstrates that using transcript-level quantification (e.g., tximport) improves gene-level DE inference.
- **Relevance:** Informs your quantification input module — support tximport-style inputs (Salmon/kallisto) for best practice.

---

## 4. Pathway-Level Meta-Analysis

#### Zeng et al. (2018) — *Comparative Pathway Integrator: A Framework of Meta-Analytic Integration of Multiple Transcriptomic Studies*
- **Journal:** Genes | **Citations:** 7
- **DOI:** https://doi.org/10.1101/444604
- **Summary:** Meta-analytic pathway enrichment across multiple studies. Handles pathway redundancy when combining multiple public pathway databases.
- **Key Finding:** Discovers novel enrichment patterns and confirms previously identified functions across psychiatric disorder datasets.
- **Relevance:** Directly relevant for your pathway meta-analysis module — study its redundancy-reduction approach.

---

## 5. Applied Examples & Validation Use Cases

#### Heberle et al. (2025) — *Systematic review and meta-analysis of bulk RNAseq studies in human Alzheimer's disease brain tissue*
- **Journal:** Alzheimer's & Dementia
- **DOI:** https://doi.org/10.1002/alz.70025
- **Summary:** Real-world bulk RNA-seq meta-analysis across 24 studies. Developed 10 quality criteria for study inclusion and performed meta-analysis on 3 high-quality datasets. Identified 571 DEGs.
- **Relevance:** Use as a template for your QC module design and as a validation dataset for benchmarking your package.

---

#### Thind et al. (2021) — *Demystifying emerging bulk RNA-Seq applications*
- **Journal:** Briefings in Bioinformatics | **Citations:** 72
- **DOI:** https://doi.org/10.1093/bib/bbab259
- **Summary:** Broad overview of bulk RNA-seq applications and bioinformatics tools — covers isoform expression, alternative splicing, SNP calling.
- **Relevance:** Useful for scoping which features to include in your package beyond standard DE analysis.

---

## 6. Recommended Package Architecture

Based on the literature gaps identified above, the following modular architecture is recommended:

```
RNAmetaR (proposed package)
│
├── Module 1: Data Ingestion & QC
│   ├── GEO/SRA automated fetching
│   ├── Metadata-based sample clustering [sampleclusteR]
│   └── Study-level QC criteria (10-point checklist) [Heberle et al. 2025]
│
├── Module 2: Normalization & Harmonization
│   ├── Cross-platform gene ID reannotation [GEDI]
│   ├── Library type correction (polyA vs rRNA-depleted) [Bush et al. 2017]
│   └── Batch effect removal (ComBat, limma::removeBatchEffect)
│
├── Module 3: Per-Study DE Analysis
│   ├── DESeq2 wrapper
│   ├── edgeR wrapper
│   ├── limma-voom wrapper [Ritchie et al. 2015]
│   └── tximport support [Soneson et al. 2015]
│
├── Module 4: Meta-Analysis Statistics
│   ├── P-value combination: Fisher, Stouffer, inverse-normal [Rau et al. 2013]
│   ├── Fused inverse-normal [Prasad & Li 2021]
│   ├── Random/fixed effects models [Keel & Lindholm-Perry 2022]
│   └── Adaptive weighting (AWmeta) [Hu et al. 2025]
│
├── Module 5: Missing Gene Handling
│   └── Imputation / partial overlap strategies [DExMA]
│
├── Module 6: Pathway Meta-Analysis
│   └── Cross-study enrichment with redundancy reduction [CPI, Zeng et al. 2018]
│
└── Module 7: Visualization & Reporting
    ├── Forest plots with heterogeneity metrics (I²)
    ├── Volcano plots (per-study + meta)
    ├── Cross-study heatmaps [hCoCena]
    └── HTML/PDF report generation
```

### Key Differentiators vs. Existing Tools
- **vs. DExMA:** Full end-to-end pipeline (GEO → results), not just statistics
- **vs. GEDI:** Includes meta-analysis statistics, not just harmonization
- **vs. hCoCena:** DE-focused, not co-expression only
- **vs. limma:** Designed specifically for cross-study integration

---

## 7. Full Reference List

| # | Authors | Year | Title | Journal | DOI |
|---|---|---|---|---|---|
| 1 | Rau, Marot, Jaffrézic | 2013 | Differential meta-analysis of RNA-seq data from multiple studies | BMC Bioinformatics | [10.1186/1471-2105-15-91](https://doi.org/10.1186/1471-2105-15-91) |
| 2 | Prasad, Li | 2021 | Fused inverse-normal method for integrated differential expression analysis of RNA-seq data | BMC Bioinformatics | [10.1186/s12859-022-04859-9](https://doi.org/10.1186/s12859-022-04859-9) |
| 3 | Hu et al. | 2025 | AWmeta empowers adaptively-weighted transcriptomic meta-analysis | bioRxiv | [10.1101/2025.05.06.650408](https://doi.org/10.1101/2025.05.06.650408) |
| 4 | Keel, Lindholm-Perry | 2022 | Recent developments and future directions in meta-analysis of differential gene expression in livestock RNA-Seq | Frontiers in Genetics | [10.3389/fgene.2022.983043](https://doi.org/10.3389/fgene.2022.983043) |
| 5 | Villatoro-García et al. | 2022 | DExMA: An R Package for Performing Gene Expression Meta-Analysis with Missing Genes | Mathematics | [10.3390/math10183376](https://doi.org/10.3390/math10183376) |
| 6 | Stokholm, Rabaglino, Kadarmideen | 2024 | GEDI: An R Package for Integration of Transcriptomic Data from Multiple Platforms | Current Protocols | [10.1002/cpz1.70046](https://doi.org/10.1002/cpz1.70046) |
| 7 | Oestreich et al. | 2022 | hCoCena: horizontal integration and analysis of transcriptomics datasets | Bioinformatics | [10.1093/bioinformatics/btac589](https://doi.org/10.1093/bioinformatics/btac589) |
| 8 | Ritchie et al. | 2015 | limma powers differential expression analyses for RNA-sequencing and microarray studies | Nucleic Acids Research | [10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007) |
| 9 | Coke, Niranjan, Ewing | 2025 | sampleclusteR: A lightweight R package for automated clustering of transcriptomics samples using metadata | bioRxiv | [10.1101/2025.04.10.648129](https://doi.org/10.1101/2025.04.10.648129) |
| 10 | Bush et al. | 2017 | Integration of quantitated expression estimates from polyA-selected and rRNA-depleted RNA-seq libraries | BMC Bioinformatics | [10.1186/s12859-017-1714-9](https://doi.org/10.1186/s12859-017-1714-9) |
| 11 | Panahi et al. | 2025 | Unifying RNA-seq data using meta-analysis; bioinformatics frameworks and application for plant genomics | Current Plant Biology | [10.1016/j.cpb.2025.100523](https://doi.org/10.1016/j.cpb.2025.100523) |
| 12 | Soneson, Love, Robinson | 2015 | Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences | F1000Research | [10.12688/f1000research.7563.2](https://doi.org/10.12688/f1000research.7563.2) |
| 13 | Zeng et al. | 2018 | Comparative Pathway Integrator: A Framework of Meta-Analytic Integration of Multiple Transcriptomic Studies | Genes | [10.1101/444604](https://doi.org/10.1101/444604) |
| 14 | Heberle et al. | 2025 | Systematic review and meta-analysis of bulk RNAseq studies in human Alzheimer's disease brain tissue | Alzheimer's & Dementia | [10.1002/alz.70025](https://doi.org/10.1002/alz.70025) |
| 15 | Thind et al. | 2021 | Demystifying emerging bulk RNA-Seq applications | Briefings in Bioinformatics | [10.1093/bib/bbab259](https://doi.org/10.1093/bib/bbab259) |
