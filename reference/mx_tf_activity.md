# Infer transcription factor activity from the meta-signature

Estimates the activity of each transcription factor (TF) from the
expression changes of its target genes, using a univariate linear model
(ULM): for each TF, the gene-level meta statistic is regressed on the
TF's signed target weights, and the t-value of the slope is the activity
score. This is the default method of decoupleR and is well suited to
CollecTRI or DoRothEA regulons.

## Usage

``` r
mx_tf_activity(
  meta_result,
  network,
  id_col = "gene_id",
  stat = c("z", "log2FC"),
  min_size = 5
)
```

## Arguments

- meta_result:

  A
  [`metaXpressResult`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressResult-class.md)
  object.

- network:

  A `data.frame` of TF-target interactions with columns `source` (TF),
  `target` (gene) and `mor` (mode of regulation: positive for
  activation, negative for repression). OmniPath exports
  (`source_genesymbol`, `target_genesymbol`, `is_stimulation`,
  `is_inhibition`) are converted automatically.

- id_col:

  Column of `meta_table` whose values match `network$target`. Default
  `"gene_id"`; use `"symbol"` when gene IDs are Entrez/Ensembl and a
  symbol column is present.

- stat:

  Gene-level statistic: `"z"` (signed meta z-score, default) or
  `"log2FC"`.

- min_size:

  Minimum number of measured targets per TF. Default 5.

## Value

A `data.frame` sorted by p-value with columns `source`, `score`
(positive = TF more active in cases), `pvalue`, `padj` (BH) and
`n_targets`.

## References

Badia-i-Mompel, P. et al. (2022) decoupleR: ensemble of computational
methods to infer biological activities from omics data. *Bioinformatics
Advances*, **2**(1), vbac016.
[doi:10.1093/bioadv/vbac016](https://doi.org/10.1093/bioadv/vbac016)

Müller-Dott, S. et al. (2023) Expanding the coverage of regulons from
high-confidence prior knowledge for accurate estimation of transcription
factor activities. *Nucleic Acids Research*, **51**(20), 10934–10949.
[doi:10.1093/nar/gkad841](https://doi.org/10.1093/nar/gkad841)

## See also

[`mx_drug_connectivity`](https://mdjubayerhossain.com/metaXpress/reference/mx_drug_connectivity.md),
[`mx_downstream_plot`](https://mdjubayerhossain.com/metaXpress/reference/mx_downstream_plot.md)

## Examples

``` r
mt <- data.frame(gene_id = paste0("G", 1:200),
                 meta_log2FC = c(rep(2, 20), rnorm(180, sd = 0.3)),
                 meta_pvalue = c(rep(1e-6, 20), runif(180)),
                 meta_padj = 1, i_squared = 0, n_studies = 3L)
res <- new("metaXpressResult", meta_table = mt, method = "fisher",
           n_studies = 3L)
net <- data.frame(source = rep(c("TF_A", "TF_B"), each = 20),
                  target = paste0("G", c(1:20, 101:120)), mor = 1)
mx_tf_activity(res, net)
#>   source      score       pvalue         padj n_targets
#> 1   TF_A 22.2375993 9.838827e-56 1.967765e-55        20
#> 2   TF_B -0.1869722 8.518739e-01 8.518739e-01        20
```
