# The metaXpressStudy class

Represents a single bulk RNA-seq study. `metaXpressStudy` extends
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html),
so raw counts live in the `"counts"` assay, sample metadata in
`colData()`, and optional gene annotation in `rowData()`. All standard
SummarizedExperiment accessors and subsetting (`study[genes, samples]`,
`study$condition`) work directly. Populated by
[`mx_fetch_geo`](https://mdjubayerhossain.com/metaXpress/reference/mx_fetch_geo.md),
[`mx_fetch_sra`](https://mdjubayerhossain.com/metaXpress/reference/mx_fetch_sra.md),
[`mx_load_local`](https://mdjubayerhossain.com/metaXpress/reference/mx_load_local.md),
or
[`mx_from_se`](https://mdjubayerhossain.com/metaXpress/reference/mx_from_se.md).

## Details

The `colData` must contain a `condition` column (case/control labels)
once an accession is set.

## Slots

- `accession`:

  A length-1 character giving the GEO/SRA accession ID (e.g.,
  `"GSE12345"`) or `"local"` for user-supplied data.

- `organism`:

  A length-1 character giving the species (e.g., `"Homo sapiens"`).

- `qc_score`:

  A length-1 numeric in \[0, 10\]. Filled by
  [`mx_qc_study`](https://mdjubayerhossain.com/metaXpress/reference/mx_qc_study.md);
  `NA` until then.

- `de_result`:

  A `data.frame` of per-study differential expression results. Empty
  until
  [`mx_de`](https://mdjubayerhossain.com/metaXpress/reference/mx_de.md)
  is called. Required columns when populated: `gene_id`, `log2FC`,
  `pvalue`, `padj`, `baseMean`.

## See also

[`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy.md),
[`mx_from_se`](https://mdjubayerhossain.com/metaXpress/reference/mx_from_se.md),
[`mx_fetch_geo`](https://mdjubayerhossain.com/metaXpress/reference/mx_fetch_geo.md),
[`mx_qc_study`](https://mdjubayerhossain.com/metaXpress/reference/mx_qc_study.md),
[`mx_de`](https://mdjubayerhossain.com/metaXpress/reference/mx_de.md)
