# Normalize counts within a single study

Applies one of five normalization methods to the raw count matrix of a
`metaXpressStudy` object.

## Usage

``` r
mx_normalize(study, method = c("TMM", "VST", "CPM", "TPM", "quantile"))
```

## Arguments

- study:

  A
  [`metaXpressStudy`](https://hossainlab.github.io/metaXpress/reference/metaXpressStudy-class.md)
  object with raw counts.

- method:

  Character scalar. Normalization method. One of `"TMM"` (default),
  `"VST"`, `"CPM"`, `"TPM"`, or `"quantile"`.

## Value

The input `study` with `counts` replaced by normalized values.

## Details

- **TMM**: Trimmed mean of M-values (Robinson & Oshlack 2010),
  implemented in edgeR. Recommended for between-sample normalization of
  count data.

- **VST**: Variance-stabilizing transformation (Anders & Huber 2010),
  implemented in DESeq2. Produces approximately homoscedastic values
  suitable for clustering and visualization.

- **CPM**: Counts per million via edgeR. Simple library-size scaling
  without variance stabilization.

- **TPM**: Transcripts per million (Li et al. 2010). Requires gene
  lengths in `metadata`. Preferred when comparing expression levels
  across genes.

- **quantile**: Quantile normalization (Bolstad et al. 2003),
  implemented in limma. Forces identical distributions across samples;
  operates on log1p-transformed counts.

## References

Robinson, M.D. & Oshlack, A. (2010) A scaling normalization method for
differential expression analysis of RNA-seq data. *Genome Biology*,
**11**(3), R25.
[doi:10.1186/gb-2010-11-3-r25](https://doi.org/10.1186/gb-2010-11-3-r25)

Anders, S. & Huber, W. (2010) Differential expression analysis for
sequence count data. *Genome Biology*, **11**(10), R106.
[doi:10.1186/gb-2010-11-10-r106](https://doi.org/10.1186/gb-2010-11-10-r106)

Bolstad, B.M. et al. (2003) A comparison of normalization methods for
high density oligonucleotide array data based on variance and bias.
*Bioinformatics*, **19**(2), 185–193.
[doi:10.1093/bioinformatics/19.2.185](https://doi.org/10.1093/bioinformatics/19.2.185)

## See also

[`mx_reannotate`](https://hossainlab.github.io/metaXpress/reference/mx_reannotate.md),
[`mx_remove_batch`](https://hossainlab.github.io/metaXpress/reference/mx_remove_batch.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  study <- mx_normalize(study, method = "TMM")
} # }
```
