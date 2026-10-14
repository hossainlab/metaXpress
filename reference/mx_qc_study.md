# Score a study against the 10-point QC checklist

Applies 10 data-level QC criteria. Each criterion is recorded as
`"pass"`, `"fail"` or `"not evaluable"` (the information needed is
absent, e.g. alignment rate is not deposited in GEO). `qc_score` is the
number of criteria *passed*. Criteria that could not be evaluated are
not counted as passed and the score is **not** rescaled, so a score is
only comparable between studies together with the number of evaluable
criteria. Decide on inclusion with
[`mx_filter_studies`](https://mdjubayerhossain.com/metaXpress/reference/mx_filter_studies.md),
which works on failures rather than the score.

## Usage

``` r
mx_qc_study(study, min_genes = 15000)
```

## Arguments

- study:

  A
  [`metaXpressStudy`](https://mdjubayerhossain.com/metaXpress/reference/metaXpressStudy-class.md)
  object.

- min_genes:

  Integer. Threshold for criterion 6. Default `15000`, appropriate for
  an unfiltered human or mouse gene-level matrix; lower it for matrices
  already restricted (e.g. to protein-coding genes).

## Value

The input `study` with `qc_score` filled (criteria passed, 0–10). The
per-criterion table (`criterion`, `name`, `status`, `score`, `details`)
is stored in `S4Vectors::metadata(study)$qc` and in
`attr(study, "qc_details")`.

## Details

The 10 criteria are:

1.  **Minimum sample size:** at least 3 samples in every condition group
    (Schurch et al. 2016).

2.  **Sequencing depth:** median library size \>= 10 million assigned
    reads.

3.  **Alignment rate:** mean \>= 70% (column `alignment_rate` in
    `colData`, fraction or percent).

4.  **rRNA contamination:** mean \< 10% (column `rrna_fraction`).

5.  **Duplicate rate:** mean \< 50% (column `duplicate_rate`).

6.  **Gene detection:** at least `min_genes` genes with a non-zero count
    in \>= 50% of samples.

7.  **Condition labels complete:** every sample has a non-missing,
    non-empty `condition`.

8.  **Clear contrast:** exactly two condition levels. Studies with more
    groups fail until they are subset to the two compared groups.

9.  **No batch-condition confounding:** batch and condition are not
    perfectly aliased (Leek et al. 2010). Not evaluable without a
    `batch` column.

10. **Raw counts:** non-negative integers with maximum \> 1,000, i.e.
    untransformed counts suitable for count models.

## References

Schurch, N.J. et al. (2016) How many biological replicates are needed in
an RNA-seq experiment and which differential expression tool should you
use? *RNA*, **22**(6), 839–851.
[doi:10.1261/rna.053959.115](https://doi.org/10.1261/rna.053959.115)

Leek, J.T. et al. (2010) Tackling the widespread and critical impact of
batch effects in high-throughput data. *Nature Reviews Genetics*,
**11**(10), 733–739.
[doi:10.1038/nrg2825](https://doi.org/10.1038/nrg2825)

Heberle, B.A. et al. (2024) Systematic review and meta-analysis of bulk
RNAseq studies in human Alzheimer's disease brain tissue. *bioRxiv*.
[doi:10.1101/2024.11.07.622520](https://doi.org/10.1101/2024.11.07.622520)
(complementary study-level reporting quality assessment).

## Examples

``` r
if (requireNamespace("airway", quietly = TRUE)) {
  data("airway", package = "airway")
  study <- mx_from_se(airway, condition_col = "dex", accession = "GSE52778")
  study <- mx_qc_study(study)
  study@qc_score
  S4Vectors::metadata(study)$qc
}
#>    criterion              name        status score
#> 1          1       sample size          pass     1
#> 2          2  sequencing depth          pass     1
#> 3          3    alignment rate not evaluable    NA
#> 4          4     rRNA fraction not evaluable    NA
#> 5          5    duplicate rate not evaluable    NA
#> 6          6    gene detection          pass     1
#> 7          7  condition labels          pass     1
#> 8          8    clear contrast          pass     1
#> 9          9 batch confounding not evaluable    NA
#> 10        10        raw counts          pass     1
#>                                                       details
#> 1                        Min samples per group: 4 (need >= 3)
#> 2                  Median assigned reads: 20.9M (need >= 10M)
#> 3                             'alignment_rate' not in colData
#> 4                              'rrna_fraction' not in colData
#> 5                             'duplicate_rate' not in colData
#> 6  Genes detected in >= 50% of samples: 24706 (need >= 15000)
#> 7                       Samples with a condition label: 8 / 8
#> 8                        Condition levels: 2 (need exactly 2)
#> 9                    No multi-level 'batch' column in colData
#> 10                                    Raw integer counts: yes
```
