# metaXpress Package Architecture Reference

## Full Module Map

```
metaXpress/
├── R/
│   ├── AllClasses.R         # S4 class definitions
│   ├── AllGenerics.R        # S4 generic definitions
│   ├── 01_ingest.R          # Module 1: Data ingestion & QC
│   ├── 02_harmonize.R       # Module 2: Normalization & harmonization
│   ├── 03_de.R              # Module 3: Per-study DE analysis
│   ├── 04_meta.R            # Module 4: Meta-analysis statistics
│   ├── 05_missing.R         # Module 5: Missing gene handling
│   ├── 06_pathway.R         # Module 6: Pathway meta-analysis
│   ├── 07_visualize.R       # Module 7: Visualization
│   ├── 08_report.R          # Module 8: Report generation
│   └── utils.R              # Internal utilities
├── tests/testthat/
│   ├── test-ingest.R
│   ├── test-harmonize.R
│   ├── test-de.R
│   ├── test-meta.R
│   ├── test-missing.R
│   ├── test-pathway.R
│   └── test-visualize.R
├── vignettes/
│   ├── quickstart.Rmd
│   └── full_workflow.Rmd
├── data/
│   └── metaXpress_example.rda   # Small example dataset for tests/vignettes
├── inst/
│   └── extdata/                 # Example GEO-format files
├── DESCRIPTION
├── NAMESPACE
└── README.md
```

---

## S4 Class Hierarchy

### metaXpressStudy
Represents a single RNA-seq study.

| Slot | Type | Description |
|------|------|-------------|
| `counts` | matrix | Raw count matrix (genes x samples). Rownames = gene IDs, colnames = sample IDs |
| `metadata` | data.frame | Sample metadata. MUST contain `condition` column (case/control labels) |
| `accession` | character | GEO/SRA accession (e.g., "GSE12345") or "local" |
| `organism` | character | Species (e.g., "Homo sapiens") |
| `qc_score` | numeric | 0-10 QC score from mx_qc_study() |
| `de_result` | data.frame | Per-study DE result. Empty data.frame until mx_de() is called |

**Validity rules:**
- `metadata` must have a `condition` column
- `ncol(counts)` must equal `nrow(metadata)`
- `qc_score` must be between 0 and 10

### metaXpressResult
Represents the output of a meta-analysis.

| Slot | Type | Description |
|------|------|-------------|
| `meta_table` | data.frame | Gene-level results (see column spec below) |
| `method` | character | Meta-analysis method used |
| `n_studies` | integer | Number of studies integrated |
| `heterogeneity` | data.frame | I-squared and Q-stat per gene |
| `pathway_result` | data.frame | Pathway meta-analysis results (empty until mx_pathway_meta() called) |

**meta_table required columns:**
- `gene_id` — gene identifier (consistent namespace)
- `meta_log2FC` — pooled effect size (log2 fold change)
- `meta_pvalue` — combined p-value
- `meta_padj` — BH-adjusted meta p-value (USE THIS for significance)
- `i_squared` — heterogeneity statistic (0-100%)
- `q_stat` — Cochran's Q statistic
- `n_studies` — number of studies with data for this gene
- `direction_consistency` — fraction of studies with same direction as meta_log2FC

---

## Data Flow

```
GEO/SRA/Local Files
        |
        v
[Module 1: mx_fetch_geo / mx_load_local]
        |
        v
List of metaXpressStudy objects (raw counts + metadata)
        |
        v
[Module 1: mx_qc_study + mx_filter_studies]
        |
        v
QC-filtered study list
        |
        v
[Module 2: mx_reannotate + mx_normalize + mx_remove_batch]
        |
        v
Harmonized study list (common gene namespace, batch-corrected)
        |
        v
[Module 3: mx_de_all]
        |
        v
List of per-study DE results (log2FC, pvalue, padj per gene)
        |
        v
[Module 5: mx_missing_summary + mx_impute / mx_filter_coverage]
        |
        v
Aligned DE result matrix (genes x studies)
        |
        v
[Module 4: mx_meta]
        |
        v
metaXpressResult (meta_log2FC, meta_padj, I-squared per gene)
        |
        v
[Module 6: mx_pathway_meta]  [Module 7: mx_volcano / mx_forest / mx_heatmap]
        |                              |
        v                              v
Pathway enrichment results       Publication-ready figures
        |
        v
[Module 8: mx_report]
        |
        v
HTML / PDF reproducible report
```

---

## Function Signatures (Complete)

### Module 1 — R/01_ingest.R
```r
mx_fetch_geo(accessions, count_type = "raw", cache_dir = tempdir())
mx_fetch_sra(srp_ids, cache_dir = tempdir())
mx_load_local(count_paths, metadata_paths, accessions = NULL, organism = "Homo sapiens")
mx_qc_study(study)                          # Returns study with qc_score slot filled
mx_cluster_samples(study)                   # Auto-cluster samples by metadata
mx_filter_studies(studies, qc_threshold = 7)
```

### Module 2 — R/02_harmonize.R
```r
mx_reannotate(studies, org = "Homo sapiens", target_id = c("SYMBOL", "ENSEMBL", "ENTREZID"))
mx_normalize(study, method = c("TMM", "VST", "CPM", "TPM", "quantile"))
mx_correct_library_type(studies)            # polyA vs rRNA-depleted correction
mx_remove_batch(studies, method = c("ComBat-seq", "ComBat", "limma", "harmony"))
mx_align_genes(studies)                     # Returns studies restricted to common genes
```

### Module 3 — R/03_de.R
```r
mx_de(study, method = c("DESeq2", "edgeR", "limma-voom"), formula = ~ condition, ...)
mx_de_all(studies, method = "DESeq2", formula = ~ condition, BPPARAM = BiocParallel::SerialParam())
mx_de_summary(de_results)                   # Summary table: n DEGs per study
```

### Module 4 — R/04_meta.R
```r
mx_meta(de_results, method = "random_effects", min_studies = 2, alpha = 0.05, ...)
mx_heterogeneity(de_results)                # Returns data.frame with I-sq, Q, tau-sq per gene
mx_sensitivity(de_results, method = "random_effects")  # Leave-one-out analysis
mx_forest(gene, de_results, meta_result = NULL)
```

### Module 5 — R/05_missing.R
```r
mx_missing_summary(de_results)              # Returns coverage matrix (genes x studies)
mx_impute(de_results, method = c("exclude", "mean", "knn", "weighted"))
mx_filter_coverage(de_results, min_studies = 2)
```

### Module 6 — R/06_pathway.R
```r
mx_pathway_meta(meta_result, db = c("Hallmarks", "KEGG", "Reactome", "GO_BP"), ...)
mx_pathway_consensus(pathway_results)       # Pathways significant in majority of studies
mx_pathway_dedup(pathway_results)           # Remove redundant pathways
mx_pathway_heatmap(pathway_results, top_n = 30)
```

### Module 7 — R/07_visualize.R
```r
mx_volcano(meta_result, padj_threshold = 0.05, lfc_threshold = 1, label_top = 10)
mx_forest(gene, de_results, meta_result = NULL)
mx_heatmap(meta_result, studies, top_n = 50)
mx_upset(de_results, padj_threshold = 0.05, lfc_threshold = 1)
mx_study_overview(studies)                  # QC metrics summary plot
mx_heterogeneity_plot(meta_result)          # I-squared distribution
```

### Module 8 — R/08_report.R
```r
mx_report(meta_result, studies, de_results, format = c("html", "pdf"), output_dir = ".")
mx_export(meta_result, format = c("csv", "excel", "rds"), output_dir = ".")
mx_session_info()                           # Capture reproducibility info
```

---

## Dependency Map

| Package | Version | Used In | Purpose |
|---------|---------|---------|---------|
| GEOquery | >= 2.68 | Module 1 | GEO data fetching |
| DESeq2 | >= 1.40 | Module 3 | Per-study DE |
| edgeR | >= 3.42 | Module 3 | Per-study DE |
| limma | >= 3.56 | Module 3 | Per-study DE + batch correction |
| sva | >= 3.48 | Module 2 | ComBat batch correction |
| tximport | >= 1.28 | Module 1 | Salmon/kallisto input |
| clusterProfiler | >= 4.8 | Module 6 | Pathway enrichment |
| msigdbr | >= 7.5 | Module 6 | Gene set databases |
| ggplot2 | >= 3.4 | Module 7 | Visualization |
| ggrepel | >= 0.9 | Module 7 | Label collision avoidance |
| ComplexHeatmap | >= 2.16 | Module 7 | Heatmaps |
| BiocParallel | >= 1.34 | Module 3 | Parallelization |
| AnnotationDbi | >= 1.62 | Module 2 | Gene ID mapping |
| org.Hs.eg.db | >= 3.17 | Module 2 | Human gene annotation |
