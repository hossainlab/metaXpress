# metaXpress: End-to-End Bulk RNA-seq Meta-Analysis

metaXpress provides a unified, reproducible pipeline for integrating
multiple bulk RNA-seq studies into a single meta-analysis. Starting from
raw GEO accession IDs, the package handles every step: data ingestion
and QC, gene-ID harmonisation, library-type correction and batch
removal, per-study differential expression, meta-analysis statistics
(p-value combination and effect-size models), missing-gene imputation,
pathway-level enrichment, and production-quality reporting.

The pipeline is built around two S4 classes that flow through every
step:
[`metaXpressStudy`](https://hossainlab.github.io/metaXpress/reference/metaXpressStudy-class.md)
(one object per study, output of ingestion) and
[`metaXpressResult`](https://hossainlab.github.io/metaXpress/reference/metaXpressResult-class.md)
(the combined meta-analysis result). All public functions are prefixed
`mx_` and parallelised via BiocParallel.

## S4 Classes

- [`metaXpressStudy`](https://hossainlab.github.io/metaXpress/reference/metaXpressStudy-class.md):

  Holds a single RNA-seq study: raw count matrix, sample metadata,
  accession ID, organism, QC score, and per-study DE result. Created by
  [`mx_fetch_geo`](https://hossainlab.github.io/metaXpress/reference/mx_fetch_geo.md),
  [`mx_fetch_sra`](https://hossainlab.github.io/metaXpress/reference/mx_fetch_sra.md),
  or
  [`mx_load_local`](https://hossainlab.github.io/metaXpress/reference/mx_load_local.md).

- [`metaXpressResult`](https://hossainlab.github.io/metaXpress/reference/metaXpressResult-class.md):

  Holds the combined meta-analysis output: gene-level statistics (meta
  log2FC, combined p-value, FDR, I\\^2\\), heterogeneity table, and
  (optionally) pathway enrichment results. Created by
  [`mx_meta`](https://hossainlab.github.io/metaXpress/reference/mx_meta.md).

## Module 1 — Data Ingestion and QC

- [`mx_fetch_geo`](https://hossainlab.github.io/metaXpress/reference/mx_fetch_geo.md):

  Download count matrices and sample metadata from GEO for one or more
  accession IDs; applies the 10-point QC checklist automatically.

- [`mx_fetch_sra`](https://hossainlab.github.io/metaXpress/reference/mx_fetch_sra.md):

  Download raw RNA-seq data for SRA project IDs (requires sratools).

- [`mx_load_local`](https://hossainlab.github.io/metaXpress/reference/mx_load_local.md):

  Construct `metaXpressStudy` objects from local count-matrix files or
  `SummarizedExperiment` objects.

- [`mx_qc_study`](https://hossainlab.github.io/metaXpress/reference/mx_qc_study.md):

  Score a study against the 10-point QC checklist derived from Heberle
  et al. (2025).

- [`mx_cluster_samples`](https://hossainlab.github.io/metaXpress/reference/mx_cluster_samples.md):

  Automatically assign samples to case/control groups using metadata
  clustering.

- [`mx_filter_studies`](https://hossainlab.github.io/metaXpress/reference/mx_filter_studies.md):

  Retain only studies whose `qc_score` meets a minimum threshold.

## Module 2 — Normalisation and Harmonisation

- [`mx_reannotate`](https://hossainlab.github.io/metaXpress/reference/mx_reannotate.md):

  Convert gene identifiers to a unified namespace (SYMBOL, ENSEMBL, or
  ENTREZID) using AnnotationDbi.

- [`mx_normalize`](https://hossainlab.github.io/metaXpress/reference/mx_normalize.md):

  Normalise within-study counts by TMM, VST, CPM, TPM, or quantile
  methods.

- [`mx_correct_library_type`](https://hossainlab.github.io/metaXpress/reference/mx_correct_library_type.md):

  Adjust for systematic bias between polyA-selected and rRNA-depleted
  libraries.

- [`mx_remove_batch`](https://hossainlab.github.io/metaXpress/reference/mx_remove_batch.md):

  Remove cross-study batch effects via ComBat-seq, ComBat, limma, or
  Harmony.

- [`mx_align_genes`](https://hossainlab.github.io/metaXpress/reference/mx_align_genes.md):

  Restrict all studies to their common gene universe.

## Module 3 — Per-Study Differential Expression

- [`mx_de`](https://hossainlab.github.io/metaXpress/reference/mx_de.md):

  Run DE on a single `metaXpressStudy` using DESeq2, edgeR, or
  limma-voom.

- [`mx_de_all`](https://hossainlab.github.io/metaXpress/reference/mx_de_all.md):

  Apply
  [`mx_de`](https://hossainlab.github.io/metaXpress/reference/mx_de.md)
  to every study in a list, in parallel via BiocParallel.

- [`mx_de_summary`](https://hossainlab.github.io/metaXpress/reference/mx_de_summary.md):

  Tabulate significant DEG counts across studies at user-specified
  thresholds.

## Module 4 — Meta-Analysis Statistics

- [`mx_meta`](https://hossainlab.github.io/metaXpress/reference/mx_meta.md):

  Combine per-study DE results using Fisher, Stouffer, inverse-normal,
  fixed-effects, random-effects (DerSimonian-Laird), or AWmeta adaptive
  weighting. Returns a
  [`metaXpressResult`](https://hossainlab.github.io/metaXpress/reference/metaXpressResult-class.md).

- [`mx_heterogeneity`](https://hossainlab.github.io/metaXpress/reference/mx_heterogeneity.md):

  Compute per-gene Cochran's Q, I\\^2\\, \\\tau^2\\, and heterogeneity
  p-values.

- [`mx_sensitivity`](https://hossainlab.github.io/metaXpress/reference/mx_sensitivity.md):

  Leave-one-out sensitivity analysis: rerun meta-analysis excluding each
  study in turn to assess result stability.

- [`mx_study_overview`](https://hossainlab.github.io/metaXpress/reference/mx_study_overview.md):

  Summarise each study's contribution to the meta-analysis.

## Module 5 — Missing Gene Handling

- [`mx_missing_summary`](https://hossainlab.github.io/metaXpress/reference/mx_missing_summary.md):

  Report gene-by-study coverage as a binary presence/absence matrix.

- [`mx_filter_coverage`](https://hossainlab.github.io/metaXpress/reference/mx_filter_coverage.md):

  Retain only genes detected in at least `min_studies` studies.

- [`mx_impute`](https://hossainlab.github.io/metaXpress/reference/mx_impute.md):

  Impute missing per-study statistics (zero, mean, median, or Bayesian
  shrinkage).

## Module 6 — Pathway Meta-Analysis

- [`mx_pathway_meta`](https://hossainlab.github.io/metaXpress/reference/mx_pathway_meta.md):

  Run cross-study ORA or GSEA against MSigDB gene sets (Hallmarks, KEGG,
  Reactome, GO) using clusterProfiler.

- [`mx_pathway_consensus`](https://hossainlab.github.io/metaXpress/reference/mx_pathway_consensus.md):

  Identify pathways enriched across a majority of individual studies.

- [`mx_pathway_dedup`](https://hossainlab.github.io/metaXpress/reference/mx_pathway_dedup.md):

  Remove redundant pathways via kappa-coefficient overlap clustering.

- [`mx_pathway_heatmap`](https://hossainlab.github.io/metaXpress/reference/mx_pathway_heatmap.md):

  Heatmap of enrichment scores or \\-\log\_{10}\\(padj) for top
  pathways.

## Module 7 — Visualisation

- [`mx_volcano`](https://hossainlab.github.io/metaXpress/reference/mx_volcano.md):

  Volcano plot of meta log2FC vs. combined \\-\log\_{10}\\(padj) with
  optional gene labels.

- [`mx_forest`](https://hossainlab.github.io/metaXpress/reference/mx_forest.md):

  Forest plot showing per-study log2FC ± 95% CI and the pooled
  meta-estimate for a single gene.

- [`mx_heatmap`](https://hossainlab.github.io/metaXpress/reference/mx_heatmap.md):

  Heatmap of log2FC values for the top *n* meta-significant genes across
  all studies.

- [`mx_heterogeneity_plot`](https://hossainlab.github.io/metaXpress/reference/mx_heterogeneity_plot.md):

  Histogram of the I\\^2\\ distribution across all tested genes.

- [`mx_upset`](https://hossainlab.github.io/metaXpress/reference/mx_upset.md):

  UpSet plot of DEG overlap across studies.

- [`mx_study_overview`](https://hossainlab.github.io/metaXpress/reference/mx_study_overview.md):

  Multi-panel QC summary: scores, sample sizes, sequencing depths, and
  gene detection rates.

## Module 8 — Reporting and Export

- [`mx_report`](https://hossainlab.github.io/metaXpress/reference/mx_report.md):

  Render a fully reproducible HTML or PDF analysis report from a
  parameterised R Markdown template.

- [`mx_export`](https://hossainlab.github.io/metaXpress/reference/mx_export.md):

  Write the gene-level meta-analysis table to CSV, Excel, or RDS.

- [`mx_session_info`](https://hossainlab.github.io/metaXpress/reference/mx_session_info.md):

  Capture R session information, timestamp, and package version for
  reproducibility records.

## Getting started

The fastest way to explore the package is the quickstart vignette:
[`vignette("quickstart", package = "metaXpress")`](https://hossainlab.github.io/metaXpress/articles/quickstart.md).
For the complete GEO-to-report workflow see
[`vignette("full_workflow", package = "metaXpress")`](https://hossainlab.github.io/metaXpress/articles/full_workflow.md).

## References

Heberle, H. et al. (2025) A comprehensive framework for quality control
and meta-analysis of bulk RNA-seq data. *Alzheimer's & Dementia*,
**21**(1), e70025.
[doi:10.1002/alz.70025](https://doi.org/10.1002/alz.70025)

Fisher, R.A. (1932) *Statistical Methods for Research Workers*. 4th edn.
Edinburgh: Oliver & Boyd.

DerSimonian, R. & Laird, N. (1986) Meta-analysis in clinical trials.
*Controlled Clinical Trials*, **7**(3), 177–188.
[doi:10.1016/0197-2456(86)90046-2](https://doi.org/10.1016/0197-2456%2886%2990046-2)

Cochran, W.G. (1954) The combination of estimates from different
experiments. *Biometrics*, **10**(1), 101–129.
[doi:10.2307/3001666](https://doi.org/10.2307/3001666)

Higgins, J.P.T. & Thompson, S.G. (2002) Quantifying heterogeneity in a
meta-analysis. *Statistics in Medicine*, **21**(11), 1539–1558.
[doi:10.1002/sim.1186](https://doi.org/10.1002/sim.1186)

Love, M.I., Huber, W. & Anders, S. (2014) Moderated estimation of fold
change and dispersion for RNA-seq data with DESeq2. *Genome Biology*,
**15**(12), 550.
[doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)

Robinson, M.D., McCarthy, D.J. & Smyth, G.K. (2010) edgeR: a
Bioconductor package for differential expression analysis of digital
gene expression data. *Bioinformatics*, **26**(1), 139–140.
[doi:10.1093/bioinformatics/btp616](https://doi.org/10.1093/bioinformatics/btp616)

Ritchie, M.E. et al. (2015) limma powers differential expression
analyses for RNA-sequencing and microarray studies. *Nucleic Acids
Research*, **43**(7), e47.
[doi:10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007)

Johnson, W.E., Li, C. & Rabinovic, A. (2007) Adjusting batch effects in
microarray expression data using empirical Bayes methods.
*Biostatistics*, **8**(1), 118–127.
[doi:10.1093/biostatistics/kxj037](https://doi.org/10.1093/biostatistics/kxj037)

## See also

Useful links:

- <https://hossainlab.github.io/metaXpress/>

- <https://github.com/hossainlab/metaXpress>

- Report bugs at <https://github.com/hossainlab/metaXpress/issues>

## Author

Md. Jubayer Hossain <contact.jubayerhossain@gmail.com>
