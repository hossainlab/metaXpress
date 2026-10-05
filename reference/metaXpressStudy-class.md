# The metaXpressStudy class

Represents a single bulk RNA-seq study, holding raw counts, sample
metadata, accession information, and QC results. Populated by
[`mx_fetch_geo`](https://hossainlab.github.io/metaXpress/reference/mx_fetch_geo.md),
[`mx_fetch_sra`](https://hossainlab.github.io/metaXpress/reference/mx_fetch_sra.md),
or
[`mx_load_local`](https://hossainlab.github.io/metaXpress/reference/mx_load_local.md).

## Slots

- `counts`:

  A numeric matrix of raw counts (genes x samples). Rownames must be
  gene identifiers; colnames must be sample identifiers.

- `metadata`:

  A `data.frame` of sample metadata. Must contain columns `condition`
  (case/control labels) and `sample_id`.

- `accession`:

  A length-1 character giving the GEO/SRA accession ID (e.g.,
  `"GSE12345"`) or `"local"` for user-supplied data.

- `organism`:

  A length-1 character giving the species (e.g., `"Homo sapiens"`).

- `qc_score`:

  A length-1 numeric in \[0, 10\]. Filled by
  [`mx_qc_study`](https://hossainlab.github.io/metaXpress/reference/mx_qc_study.md);
  `NA` until then.

- `de_result`:

  A `data.frame` of per-study differential expression results. Empty
  until
  [`mx_de`](https://hossainlab.github.io/metaXpress/reference/mx_de.md)
  is called. Required columns when populated: `gene_id`, `log2FC`,
  `pvalue`, `padj`, `baseMean`.

## See also

[`mx_fetch_geo`](https://hossainlab.github.io/metaXpress/reference/mx_fetch_geo.md),
[`mx_qc_study`](https://hossainlab.github.io/metaXpress/reference/mx_qc_study.md),
[`mx_de`](https://hossainlab.github.io/metaXpress/reference/mx_de.md)
