inst/extdata/ — Example external data files
==========================================

This directory contains small example files used in vignettes and tests
that demonstrate data loading via mx_load_local().

Files:
  (to be added as the package develops)

Format notes:
  - Count matrices: genes as rows (with rownames), samples as columns.
    Accepted formats: .csv, .tsv, .txt, .gz, .rds
  - Metadata tables: samples as rows, variables as columns.
    Must contain at minimum: sample_id, condition
    Accepted formats: .csv, .tsv, .rds
