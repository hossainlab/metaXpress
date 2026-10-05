# The metaXpressResult class

Holds the output of a meta-analysis performed by
[`mx_meta`](https://hossainlab.github.io/metaXpress/reference/mx_meta.md),
including gene-level statistics, heterogeneity estimates, and
(optionally) pathway enrichment results.

## Slots

- `meta_table`:

  A `data.frame` with one row per gene. Required columns: `gene_id`,
  `meta_log2FC`, `meta_pvalue`, `meta_padj`, `i_squared`, `q_stat`,
  `n_studies`, `direction_consistency`.

- `method`:

  A length-1 character naming the meta-analysis method used. One of
  `"fisher"`, `"stouffer"`, `"inverse_normal"`, `"fixed_effects"`,
  `"random_effects"`, `"awmeta"`.

- `n_studies`:

  A length-1 integer giving the total number of studies integrated.

- `heterogeneity`:

  A `data.frame` with per-gene heterogeneity statistics: `gene_id`, `Q`,
  `df`, `I_sq`, `tau_sq`, `p_heterogeneity`.

- `pathway_result`:

  A `data.frame` of pathway meta-analysis results. Empty until
  [`mx_pathway_meta`](https://hossainlab.github.io/metaXpress/reference/mx_pathway_meta.md)
  is called.

## See also

[`mx_meta`](https://hossainlab.github.io/metaXpress/reference/mx_meta.md),
[`mx_heterogeneity`](https://hossainlab.github.io/metaXpress/reference/mx_heterogeneity.md),
[`mx_pathway_meta`](https://hossainlab.github.io/metaXpress/reference/mx_pathway_meta.md)
