inst/extdata/ — Example external data files
==========================================

This directory contains small example files used in vignettes and tests
that demonstrate data loading via mx_load_local().

Files:
  case_study_subset.rds
      List of three metaXpressStudy objects (GSE130688 CRC, GSE136569 PDAC,
      GSE171485 ccRCC; 52 samples) restricted to 400 genes. Used by the live
      tutorial for QC, DE and meta-analysis.

  case_study_downstream.rds
      Inputs for the Module 9 demo in the live tutorial, from the full case
      study. A list with:
        meta_result  genome-wide random-effects meta-analysis (29,560 genes)
        studies      the three cohorts restricted to MCP-counter markers and
                     the 1,500 most variable genes (colData$lib_size holds the
                     true library size; rowData$symbol holds gene symbols)
        markers      MCP-counter marker genes (Becht et al. 2016)
        collectri    CollecTRI TF-target network (Mueller-Dott et al. 2023)
        lincs        200 of 5,425 LINCS L1000 consensus drug signatures
                     (Enrichr), selected from the full ranking for the demo
      Built by scripts/make_downstream_extdata.R after
      scripts/fetch_case_study_data.R and scripts/run_case_study_applications.R.

Format notes:
  - Count matrices: genes as rows (with rownames), samples as columns.
    Accepted formats: .csv, .tsv, .txt, .gz, .rds
  - Metadata tables: samples as rows, variables as columns.
    Must contain at minimum: sample_id, condition
    Accepted formats: .csv, .tsv, .rds
