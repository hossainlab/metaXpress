# Meta-analysis volcano plot

Creates a volcano plot from a `metaXpressResult` object, displaying
meta-analysis effect sizes (log2 fold-change) against statistical
significance (-log10 meta_padj). The top differentially expressed genes
are labelled. Volcano plots are a standard visualization for
differential expression results (Li 2012).

## Usage

``` r
mx_volcano(
  meta_result,
  padj_threshold = 0.05,
  lfc_threshold = 1,
  label_top = 10,
  title = NULL
)
```

## Arguments

- meta_result:

  A
  [`metaXpressResult`](https://hossainlab.github.io/metaXpress/reference/metaXpressResult-class.md)
  object.

- padj_threshold:

  Numeric. Adjusted p-value threshold for colouring significant genes.
  Default: `0.05`.

- lfc_threshold:

  Numeric. Log2 fold-change threshold. Default: `1`.

- label_top:

  Integer. Number of top significant genes to label. Default: `10`.

- title:

  Character scalar. Plot title. Default: `NULL` (auto).

## Value

A `ggplot2` object.

## References

Li, W. (2012) Volcano plots in analyzing differential expressions with
mRNA microarrays. *Journal of Bioinformatics and Computational Biology*,
**10**(6), 1231003.
[doi:10.1142/S0219720012310038](https://doi.org/10.1142/S0219720012310038)

## See also

[`mx_forest`](https://hossainlab.github.io/metaXpress/reference/mx_forest.md),
[`mx_heatmap`](https://hossainlab.github.io/metaXpress/reference/mx_heatmap.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  p <- mx_volcano(meta_result)
  print(p)
} # }
```
