# Retrieve gene sets from MSigDB

Downloads a gene set collection with msigdbr and returns it as a named
list, using the identifier type that matches your gene IDs.

## Usage

``` r
mx_gene_sets(
  db = c("Hallmarks", "Reactome", "KEGG", "WikiPathways", "GO_BP"),
  species = "Homo sapiens",
  id_type = c("symbol", "entrez", "ensembl")
)
```

## Arguments

- db:

  Character scalar. One of `"Hallmarks"` (default), `"Reactome"`,
  `"KEGG"` (KEGG MEDICUS), `"WikiPathways"`, `"GO_BP"`.

- species:

  Species name understood by msigdbr. Default: `"Homo sapiens"`.

- id_type:

  Gene identifier: `"symbol"` (default), `"entrez"` or `"ensembl"`. NCBI
  GEO counts from
  [`mx_fetch_geo`](https://mdjubayerhossain.com/metaXpress/reference/mx_fetch_geo.md)
  use Entrez IDs.

## Value

A named list of character vectors (gene set name -\> gene IDs). The
MSigDB version is stored in `attr(x, "version")`.

## References

Liberzon, A. et al. (2015) The Molecular Signatures Database (MSigDB)
hallmark gene set collection. *Cell Systems*, **1**(6), 417–425.
[doi:10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004)

## See also

[`mx_pathway_meta`](https://mdjubayerhossain.com/metaXpress/reference/mx_pathway_meta.md),
[`mx_pathway_studies`](https://mdjubayerhossain.com/metaXpress/reference/mx_pathway_studies.md)

## Examples

``` r
if (requireNamespace("msigdbr", quietly = TRUE)) {
  hallmarks <- mx_gene_sets("Hallmarks")
  length(hallmarks)
  head(hallmarks$HALLMARK_INFLAMMATORY_RESPONSE)
}
#> Downloading gene sets (first use only, may take a few minutes)...
#> Warning: The `category` argument of `msigdbr()` is deprecated as of msigdbr 10.0.0.
#> ℹ Please use the `collection` argument instead.
#> ℹ The deprecated feature was likely used in the metaXpress package.
#>   Please report the issue at <https://github.com/hossainlab/metaXpress/issues>.
#> Downloading gene sets (first use only, may take a few minutes)...
#> [1] "ABCA1"  "ABI1"   "ACVR1B" "ACVR2A" "ADGRE1" "ADM"   
```
