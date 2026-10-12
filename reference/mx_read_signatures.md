# Read up/down drug signatures from a gene set library

Parses a GMT file or an Enrichr gene set library text file and pairs
sets named `"<drug> up"` and `"<drug> down"` (case-insensitive,
separated by a space, underscore or hyphen) into signatures for
[`mx_drug_connectivity`](https://mdjubayerhossain.com/metaXpress/reference/mx_drug_connectivity.md).
Works with Enrichr's `LINCS_L1000_Chem_Pert_Consensus_Sigs` library.

## Usage

``` r
mx_read_signatures(path)
```

## Arguments

- path:

  Path to a GMT or Enrichr library file (tab-separated: set name,
  description, genes...).

## Value

A named list; each element is a list with character vectors `up` and
`down`. Drugs lacking either direction are dropped.

## See also

[`mx_drug_connectivity`](https://mdjubayerhossain.com/metaXpress/reference/mx_drug_connectivity.md)

## Examples

``` r
f <- tempfile(fileext = ".gmt")
writeLines(c("drugA up\t\tG1\tG2\tG3", "drugA down\t\tG4\tG5"), f)
mx_read_signatures(f)
#> $drugA
#> $drugA$up
#> [1] "G1" "G2" "G3"
#> 
#> $drugA$down
#> [1] "G4" "G5"
#> 
#> 
```
