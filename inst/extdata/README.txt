inst/extdata/ - real example data used by vignettes, examples and tests
=======================================================================

ad_brain_de.rds
    Per-study DESeq2 results (AD vs control) for four independent human
    post-mortem brain bulk RNA-seq cohorts, all quantified by the NCBI GEO
    RNA-seq pipeline (GRCh38.p13, NCBI RefSeq genes):
      GSE125583  fusiform gyrus         218 AD / 70 control  (~ sex + condition)
      GSE104704  lateral temporal lobe   12 AD / 10 aged control
      GSE53697   dorsolateral PFC         6 AD /  6 control (12 of 17 libraries
                                          processed by NCBI)
      GSE173955  hippocampus              8 AD / 10 control
    GSE95587 (89/117 donors shared with GSE125583) and GSE159699 (same donors
    as GSE104704) were excluded to keep studies independent.
    A list with:
      de_results  named list of data.frames (gene_id = HGNC symbol,
                  protein-coding genes passing >= 10 counts in >= smaller-
                  group-size samples; log2FC, lfcSE, stat, pvalue, padj,
                  baseMean, n_case, n_control)
      studies     per-study summary table
      provenance  build date, count source, DE engine and filter
    Built by inst/scripts/make_ad_brain_extdata.R.

case_study_downstream.rds
    Inputs for the Module 9 vignette (downstream_biology), from a
    pancreatic ductal adenocarcinoma (PDAC) tumour vs adjacent non-tumour
    meta-analysis of three GEO cohorts: GSE130688 (15/15), GSE136569 (5/5)
    and GSE171485 (6/6); NCBI GEO-generated counts, labels from GEO series
    metadata. A list with:
      meta_result  genome-wide random-effects (REML) meta-analysis
                   (29,679 genes)
      studies      the three cohorts restricted to MCP-counter markers and
                   the 1,500 most variable genes (colData$lib_size holds the
                   true library size; rowData$symbol holds gene symbols)
      markers      MCP-counter marker genes (Becht et al. 2016)
      collectri    CollecTRI TF-target network (Mueller-Dott et al. 2023)
      lincs        200 of 5,425 LINCS L1000 consensus drug signatures
                   (Enrichr), selected from the full ranking for the demo
    Built by scripts/make_downstream_extdata.R (source repository).

hallmark_gene_sets.rds
    The 50 MSigDB Hallmark gene sets (human gene symbols), a named list of
    character vectors, so that vignettes and examples run without internet.
    MSigDB Hallmark collection (Liberzon et al. 2015, Cell Systems 1:417,
    doi:10.1016/j.cels.2015.12.004), distributed under CC BY 4.0. Version
    recorded in attr(x, "source"). Built by
    inst/scripts/make_hallmark_extdata.R.
