# Plot downstream biology results

Draws a ranked bar chart for any result of
[`mx_tf_activity`](https://mdjubayerhossain.com/metaXpress/reference/mx_tf_activity.md),
[`mx_cell_composition`](https://mdjubayerhossain.com/metaXpress/reference/mx_cell_composition.md),
[`mx_coexpression_modules`](https://mdjubayerhossain.com/metaXpress/reference/mx_coexpression_modules.md)
or
[`mx_drug_connectivity`](https://mdjubayerhossain.com/metaXpress/reference/mx_drug_connectivity.md).
Bars show the score (TF activity, pooled Hedges' g or reversal score);
faded bars are not significant.

## Usage

``` r
mx_downstream_plot(x, top = 20, padj_threshold = 0.05, title = NULL)
```

## Arguments

- x:

  A result from one of the downstream functions.

- top:

  Number of entries to show. For TF activity and drug connectivity the
  top entries in each direction are shown (`top / 2` each). Default 20.

- padj_threshold:

  Significance threshold for bar shading. Default 0.05.

- title:

  Optional plot title.

## Value

A `ggplot2` object.

## Examples

``` r
if (FALSE) { # \dontrun{
  mx_downstream_plot(mx_tf_activity(meta_result, collectri))
} # }
```
