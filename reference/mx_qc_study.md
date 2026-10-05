# Score a study against the 10-point QC checklist

Applies each of the 10 QC criteria defined by Heberle et al. (2025) and
stores the total score in the `qc_score` slot. Criteria that cannot be
evaluated (e.g., alignment rate not in metadata) are scored `NA` and do
not penalise the total.

## Usage

``` r
mx_qc_study(study)
```

## Arguments

- study:

  A
  [`metaXpressStudy`](https://hossainlab.github.io/metaXpress/reference/metaXpressStudy-class.md)
  object.

## Value

The input `study` with the `qc_score` slot filled. A `qc_details`
attribute is attached with a per-criterion breakdown.

## Details

The 10 criteria are:

1.  **Minimum sample size:** at least 3 samples per condition group,
    required for reliable dispersion estimation in negative binomial
    models (Schurch et al. 2016).

2.  **Sequencing depth:** median library size \>= 10 million reads,
    ensuring adequate power for detecting low-abundance transcripts
    (Sims et al. 2014).

3.  **Alignment rate:** mean alignment rate \>= 70%, indicating
    acceptable read quality (optional; skipped if not available).

4.  **rRNA contamination:** mean rRNA fraction \< 10%, flagging
    insufficient ribosomal depletion (optional).

5.  **Duplicate rate:** mean duplicate rate \< 50%, flagging library
    complexity issues (optional).

6.  **Gene detection:** at least 15,000 genes detected in \>= 50% of
    samples, confirming adequate transcriptome coverage.

7.  **Metadata completeness:** required columns `condition` and
    `sample_id` are present.

8.  **Clear case/control:** exactly two condition levels, ensuring a
    well-defined contrast for differential expression.

9.  **No batch-condition confounding:** batch and condition are not
    perfectly aliased, which would make batch correction impossible
    (Leek et al. 2010).

10. **Raw counts:** data are non-negative integers with maximum \>
    1,000, confirming untransformed count data suitable for count-based
    statistical models.

## References

Heberle, H. et al. (2025) A comprehensive framework for quality control
and meta-analysis of bulk RNA-seq data. *Alzheimer's & Dementia*,
**21**(1), e70025.
[doi:10.1002/alz.70025](https://doi.org/10.1002/alz.70025)

Schurch, N.J. et al. (2016) How many biological replicates are needed in
an RNA-seq experiment and which differential expression tool should you
use? *RNA*, **22**(6), 839–851.
[doi:10.1261/rna.053959.115](https://doi.org/10.1261/rna.053959.115)

Leek, J.T. et al. (2010) Tackling the widespread and critical impact of
batch effects in high-throughput data. *Nature Reviews Genetics*,
**11**(10), 733–739.
[doi:10.1038/nrg2825](https://doi.org/10.1038/nrg2825)

## Examples

``` r
if (FALSE) { # \dontrun{
  study <- mx_qc_study(study)
  attr(study, "qc_details")
} # }
```
