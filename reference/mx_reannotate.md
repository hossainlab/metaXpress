# Reannotate gene IDs to a common namespace

Converts the gene identifiers of one or more studies to a single
namespace (gene symbol, Ensembl gene ID or Entrez ID) so that studies
can be compared gene by gene.

## Usage

``` r
mx_reannotate(
  studies,
  org = "Homo sapiens",
  target_id = c("SYMBOL", "ENSEMBL", "ENTREZID"),
  from_id = NULL,
  mapping = NULL,
  duplicates = c("max_mean", "sum", "drop"),
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- studies:

  A named list of
  [`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
  objects, or a single study.

- org:

  Organism: `"Homo sapiens"` (default; org.Hs.eg.db), `"Mus musculus"`
  (org.Mm.eg.db) or `"Rattus norvegicus"` (org.Rn.eg.db). Ignored when
  `mapping` is given.

- target_id:

  Target identifier: `"SYMBOL"` (default), `"ENSEMBL"` or `"ENTREZID"`.

- from_id:

  Source identifier type, or `NULL` (default) to detect it per study.

- mapping:

  Optional `data.frame` with columns `from` and `to` giving the ID
  mapping to use instead of an `org.*.eg.db` package.

- duplicates:

  How to collapse several source IDs mapping to one target: `"max_mean"`
  (default), `"sum"` or `"drop"`.

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html).
  Default:
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html).

## Value

The studies with translated rownames (same type as the input: a list or
a single study). A per-study mapping summary is printed.

## Details

Processing steps for each study:

1.  Ensembl version suffixes are removed (`ENSG00000141510.17` becomes
    `ENSG00000141510`), as are Ensembl `_PAR_Y` pseudo-autosomal
    duplicates.

2.  The source identifier type is detected from *all* IDs by majority
    (Ensembl, Entrez or symbol), or set with `from_id`.

3.  IDs are mapped with the organism's `org.*.eg.db` package, or with a
    user-supplied `mapping` table (any organism, offline).

4.  Unmapped IDs are dropped and reported.

5.  When several source IDs map to the same target ID, they are
    collapsed with `duplicates`: `"max_mean"` keeps the most highly
    expressed source row (default, as in WGCNA's `collapseRows`);
    `"sum"` adds counts (appropriate when the source rows are disjoint
    features of one gene); `"drop"` removes all ambiguous genes.

The original identifier of every retained row is stored in
`rowData(study)$source_id`.

## References

Pages, H. et al. AnnotationDbi: manipulation of SQLite-based annotations
in Bioconductor. <https://bioconductor.org/packages/AnnotationDbi>

Miller, J.A. et al. (2011) Strategies for aggregating gene expression
data: the collapseRows R function. *BMC Bioinformatics*, **12**, 322.
[doi:10.1186/1471-2105-12-322](https://doi.org/10.1186/1471-2105-12-322)

## See also

[`mx_align_genes`](https://mdjubayerhossain.com/metaXpress/reference/mx_align_genes.md),
[`mx_gene_sets`](https://mdjubayerhossain.com/metaXpress/reference/mx_gene_sets.md)

## Examples

``` r
counts <- matrix(c(10, 20, 30, 5, 1, 2, 3, 4, 7, 7, 7, 7),
                 nrow = 3, byrow = TRUE,
                 dimnames = list(c("ENSG00000141510.17", "ENSG00000012048.23",
                                   "ENSG00000999999.1"),
                                 paste0("S", 1:4)))
meta  <- data.frame(condition = c("control", "control", "case", "case"))
study <- metaXpressStudy(counts, meta, accession = "local")
map   <- data.frame(from = c("ENSG00000141510", "ENSG00000012048"),
                    to   = c("TP53", "BRCA1"))
out <- mx_reannotate(study, mapping = map)
#>   local: ENSEMBL -> SYMBOL | 3 IDs, 1 unmapped (33.3%), 0 rows collapsed (max_mean), 2 genes
rownames(out)
#> [1] "TP53"  "BRCA1"
SummarizedExperiment::rowData(out)$source_id
#> [1] "ENSG00000141510.17" "ENSG00000012048.23"
```
