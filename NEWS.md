# metaXpress 0.99.0

* Initial Bioconductor submission.
* Module 1: Data ingestion from GEO/SRA and local files with 10-point QC scoring.
* Module 2: Normalization (TMM, VST, CPM, TPM, quantile) and batch correction
  (ComBat-seq, ComBat, limma, harmony).
* Module 3: Per-study differential expression via DESeq2, edgeR, and limma-voom.
* Module 4: Meta-analysis statistics — Fisher, Stouffer, inverse-normal,
  fixed effects, random effects (DerSimonian-Laird), and AWmeta.
* Module 5: Missing gene handling with imputation and coverage filtering.
* Module 6: Pathway meta-analysis using KEGG, Reactome, GO, and MSigDB Hallmarks.
* Module 7: Publication-ready visualization (volcano, forest, heatmap, UpSet plots).
* Module 8: Reproducible HTML/PDF report generation and result export.
