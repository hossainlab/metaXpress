# An introduction to bulk RNA-seq meta-analysis

This article is for researchers who are new to combining gene expression
studies. It explains:

- what bulk RNA-seq measures and why it is still widely used;
- why single studies are often not enough;
- what a meta-analysis adds, and how its main statistical ideas work;
- the decisions that make a meta-analysis credible.

Every real-data figure is computed when this article is built. It uses
four published Alzheimer’s disease (AD) brain cohorts that ship with
metaXpress. Figures that illustrate a concept are labelled as such. For
hands-on code, continue with
[`vignette("quickstart")`](https://mdjubayerhossain.com/metaXpress/articles/quickstart.md)
after reading.

## 1. What is bulk RNA-seq?

**RNA sequencing (RNA-seq)** measures which genes are active in a
sample, and how strongly, by sequencing its RNA molecules ([Wang et al.
2009](#ref-wang2009rnaseq); [Mortazavi et al.
2008](#ref-mortazavi2008mapping)).

In *bulk* RNA-seq, the RNA of many thousands to millions of cells is
extracted together. The result is one expression profile per sample (a
piece of tissue, a blood draw, a cell culture). Each gene’s value is the
average over all cells in that sample.

![\*\*Figure 1. The bulk RNA-seq workflow.\*\* Reads from each sample
are aligned or pseudo-aligned to a reference and counted per gene. The
resulting count matrix (genes × samples) is the input to differential
expression (DE) analysis.
Schematic.](introduction_files/figure-html/workflow-1.png)

**Figure 1. The bulk RNA-seq workflow.** Reads from each sample are
aligned or pseudo-aligned to a reference and counted per gene. The
resulting count matrix (genes × samples) is the input to differential
expression (DE) analysis. Schematic.

The main processing steps have well-established tools and published best
practices ([Conesa et al. 2016](#ref-conesa2016survey); [Stark et al.
2019](#ref-stark2019teenage)):

- **Read alignment or quantification:** STAR ([Dobin et al.
  2013](#ref-dobin2013star)), HISAT2 ([Kim et al.
  2019](#ref-kim2019hisat2)) or Salmon ([Patro et al.
  2017](#ref-patro2017salmon)).
- **Gene-level counting:** featureCounts ([Liao et al.
  2014](#ref-liao2014featurecounts)).
- **Statistical testing:** DESeq2 ([Love et al.
  2014](#ref-love2014deseq2)), edgeR ([Robinson et al.
  2010](#ref-robinson2010edger)) or limma-voom ([Law et al.
  2014](#ref-law2014voom); [Ritchie et al.
  2015](#ref-ritchie2015limma)). These use count-based models, and their
  output is a list of differentially expressed (DE) genes with a log2
  fold change and a p-value for each gene.

### Why bulk, when single-cell RNA-seq exists?

Single-cell RNA-seq resolves individual cell types, but bulk RNA-seq
remains the workhorse of clinical and translational research ([Haque et
al. 2017](#ref-haque2017practical); [Stark et al.
2019](#ref-stark2019teenage)):

|  | Bulk RNA-seq | Single-cell RNA-seq |
|----|----|----|
| Measures | Average expression of a sample | Expression in each cell |
| Cost per sample | Low | High |
| Genes detected per sample | ~15,000–20,000 | ~1,000–5,000 per cell |
| Typical cohort size | Tens to thousands | A few to tens of donors |
| Works on archived / degraded tissue | Often | Rarely |
| Public data available | Very large (GEO, SRA, TCGA, GTEx) | Growing |
| Best for | Disease vs control, biomarkers, cohorts, clinical samples | Cell types, cell states, development |

Bulk RNA-seq is also one of the largest public data resources in
biology:

- **Large reference projects:** The Cancer Genome Atlas ([Cancer Genome
  Atlas Research Network et al. 2013](#ref-tcga2013pancancer)) and the
  Genotype-Tissue Expression project ([GTEx Consortium
  2013](#ref-lonsdale2013gtex)).
- **Public archives:** many thousands of smaller studies in GEO ([Clough
  et al. 2024](#ref-clough2024geo)).
- **Uniformly reprocessed collections:** recount2/3 ([Collado-Torres et
  al. 2017](#ref-collado2017recount2); [Wilks et al.
  2021](#ref-wilks2021recount3)) and ARCHS4 ([Lachmann et al.
  2018](#ref-lachmann2018archs4)). These make the archived data reusable
  at scale.

This archive is what makes meta-analysis possible.

## 2. Why one study is usually not enough

Most bulk RNA-seq studies are small. A typical design has 3 to 10
samples per group. Each study also tests about 20,000 genes at once and
must correct for multiple testing, usually with the false discovery rate
(FDR) ([Benjamini and Hochberg 1995](#ref-benjamini1995fdr)). This has
two consequences ([Schurch et al. 2016](#ref-schurch2016replicates)):

1.  **Many true effects are missed.** A gene that changes by 40% may not
    reach FDR significance with 5 samples per group.
2.  **Many reported effects do not replicate.** Effects can come from:
    - the particular cohort, tissue sampling or laboratory protocol;
    - batch effects ([Leek et al. 2010](#ref-leek2010batch));
    - chance.

    Small studies with many tests are especially prone to false findings
    ([Ioannidis 2005](#ref-ioannidis2005false)).

The four AD cohorts in metaXpress show the problem clearly. They study
the same disease, with the same NCBI processing pipeline, but report
very different numbers of DE genes:

![\*\*Figure 2. The same disease, very different results.\*\* Genes at
FDR ≤ 0.05 in four independent AD brain cohorts (DESeq2), and the number
at meta-analysis FDR ≤ 0.05 when the four are combined with a
random-effects model. Real
data.](introduction_files/figure-html/per-study-1.png)

**Figure 2. The same disease, very different results.** Genes at FDR ≤
0.05 in four independent AD brain cohorts (DESeq2), and the number at
meta-analysis FDR ≤ 0.05 when the four are combined with a
random-effects model. Real data.

The counts in Figure 2 vary mainly with sample size. The largest cohort,
218 AD vs 70 controls, finds over 10,000 genes. The smallest, 6 vs 6,
finds one. None of the four on its own tells us which changes are
general features of AD.

## 3. What a meta-analysis adds

A **meta-analysis** combines the results of several independent studies
that address the same question ([Ramasamy et al.
2008](#ref-ramasamy2008key); [Tseng et al. 2012](#ref-tseng2012review)).
For gene expression, each study is first analysed on its own (the same
DE analysis as above). The per-gene results are then combined.

Combining results, rather than pooling raw counts into one large
dataset, means that differences between studies do not bias the
per-study estimates. These include platforms, depth, library protocols
and laboratory batch effects. The combination then provides:

- **Power.** Several modest, consistent effects add up.
- **Reproducibility.** A gene is reported because it changes in several
  independent cohorts, not because of one lucky study ([Sweeney et al.
  2017](#ref-sweeney2017methods)).
- **Generality.** Effects are estimated across populations, tissues and
  protocols, together with how much they vary.

The fold changes of independent cohorts are already related. The figure
below compares the two largest AD cohorts gene by gene:

![\*\*Figure 3. Independent studies agree more than their gene lists
suggest.\*\* Log2 fold changes (AD vs control) of genes significant in
either of the two cohorts. Points are coloured by whether the gene is
significant in the meta-analysis. Real
data.](introduction_files/figure-html/concordance-1.png)

**Figure 3. Independent studies agree more than their gene lists
suggest.** Log2 fold changes (AD vs control) of genes significant in
either of the two cohorts. Points are coloured by whether the gene is
significant in the meta-analysis. Real data.

[`mx_study_concordance()`](https://mdjubayerhossain.com/metaXpress/reference/mx_study_concordance.md)
computes these correlations for every pair of studies. A *negative*
correlation is a warning sign: it usually means that one study’s
case/control contrast is the wrong way round.

### How much power does combining studies give?

The next figure is a simple illustration, not real data. It shows a gene
with a true log2 fold change of 0.5 and a within-group standard
deviation of 1 (log2 scale). The curves give the chance of detecting
this gene at a significance level typical of genome-wide FDR control.

![\*\*Figure 4. Illustration: combining studies increases power.\*\*
Probability of detecting a log2 fold change of 0.5 (SD = 1) at α = 10⁻⁴,
as a function of samples per group in each study, for one study and for
an inverse-variance meta-analysis of 2, 4 or 8 equally sized studies
with no between-study heterogeneity. Theoretical (normal
approximation).](introduction_files/figure-html/power-1.png)

**Figure 4. Illustration: combining studies increases power.**
Probability of detecting a log2 fold change of 0.5 (SD = 1) at α = 10⁻⁴,
as a function of samples per group in each study, for one study and for
an inverse-variance meta-analysis of 2, 4 or 8 equally sized studies
with no between-study heterogeneity. Theoretical (normal approximation).

There is a caveat. Real studies differ from one another, and these
differences reduce the gain, which is why heterogeneity matters (Section
4).

## 4. Key statistical ideas

### Two ways to combine studies

Meta-analysis methods fall into two families, which answer different
questions:

1.  **Effect-size methods** combine each study’s log2 fold change
    $`y_i`$, weighting it by its precision (1 / squared standard error).
    They estimate *how large* the change is.
2.  **P-value combination methods** combine each study’s evidence for a
    change. Examples are Fisher’s and Stouffer’s methods, and
    order-statistic methods such as rOP ([Song and Tseng
    2014](#ref-song2014rop)) and AW-Fisher ([Li and Tseng
    2011](#ref-li2011awfisher); [Hu et al. 2026](#ref-hu2026awmeta)).
    They answer *whether* a gene changes in one, several or all studies,
    without estimating one common effect ([Rau et al.
    2014](#ref-rau2014metarnaseq); [Tseng et al.
    2012](#ref-tseng2012review)).

### Fixed effect versus random effects

A **fixed (common) effect** model assumes that every study measures the
same true effect, and that the studies differ only by sampling error.

A **random-effects** model assumes that the true effect itself varies
between studies: different tissues, populations and protocols have
genuinely different effects. It estimates the *average* effect and the
between-study variance $`\tau^2`$([DerSimonian and Laird
1986](#ref-dersimonian1986meta); [Borenstein et al.
2010](#ref-borenstein2010basic)).

For combining independent public cohorts, random effects is usually the
honest choice. metaXpress estimates $`\tau^2`$ by restricted maximum
likelihood (REML), which performs better than the classical
DerSimonian–Laird estimator ([Viechtbauer
2005](#ref-viechtbauer2005bias); [Langan et al.
2019](#ref-langan2019comparison)).

### Heterogeneity and the prediction interval

Two numbers describe how much studies disagree:

- **$`I^2`$** is the percentage of variation between studies that is due
  to real differences rather than chance ([Higgins and Thompson
  2002](#ref-higgins2002quantifying)).
- **The 95% prediction interval** gives the range in which the effect in
  a *new* study is expected to fall ([Higgins et al.
  2009](#ref-higgins2009reevaluation)).

If the prediction interval excludes zero, the change is expected to
reproduce in new cohorts. If it crosses zero, the change is significant
on average but not reliably present everywhere.

The two forest plots below show both situations in real data:

![\*\*Figure 5. Reading a forest plot.\*\* Each square is one study's
log2 fold change, with its 95% confidence interval; square size reflects
the study's weight. The red diamond is the pooled random-effects
estimate and its 95% CI, and the dotted line is the 95% prediction
interval. GFAP (astrocyte marker) is consistent across all cohorts. SST
(somatostatin interneurons) is lost strongly in two cohorts but not in
the others: high I² and a wide prediction interval. Real
data.](introduction_files/figure-html/forest-1.png)

**Figure 5. Reading a forest plot.** Each square is one study’s log2
fold change, with its 95% confidence interval; square size reflects the
study’s weight. The red diamond is the pooled random-effects estimate
and its 95% CI, and the dotted line is the 95% prediction interval. GFAP
(astrocyte marker) is consistent across all cohorts. SST (somatostatin
interneurons) is lost strongly in two cohorts but not in the others:
high I² and a wide prediction interval. Real data.

In the AD meta-analysis, 7,491 genes are significant at meta FDR ≤ 0.05.
For 2,133 of them the prediction interval excludes zero. These are the
changes most likely to transfer to a new AD cohort.

### Robustness: does one study drive the result?

A meta-analysis can be dominated by its largest study. **Leave-one-out**
analysis repeats the meta-analysis once without each study in turn.

![\*\*Figure 6. Leave-one-out analysis.\*\* Genes at meta FDR ≤ 0.05
with all four AD cohorts and with each cohort removed in turn. The
dashed line marks the robust core: genes that stay significant, in the
same direction, whichever study is removed. Real
data.](introduction_files/figure-html/loo-1.png)

**Figure 6. Leave-one-out analysis.** Genes at meta FDR ≤ 0.05 with all
four AD cohorts and with each cohort removed in turn. The dashed line
marks the robust core: genes that stay significant, in the same
direction, whichever study is removed. Real data.

Removing the large fusiform-gyrus cohort removes most significant genes,
because much of the evidence comes from that one study. The robust core
is the most defensible result to report.

## 5. What meta-analysis is used for

Gene-expression meta-analysis has produced results that single studies
could not:

- **Diagnostic signatures.** Multi-cohort analyses have identified small
  gene sets that diagnose a condition in new, independent cohorts:
  - sepsis and sterile inflammation ([Sweeney et al.
    2015](#ref-sweeney2015sepsis));
  - active tuberculosis ([Sweeney et al. 2016](#ref-sweeney2016tb)).
- **Shared disease mechanisms.** A “common rejection module” of genes is
  up-regulated in acute rejection of several transplanted organs. It
  pointed to existing drugs as candidate therapies ([Khatri et al.
  2013](#ref-khatri2013transplant)).
- **Drug repurposing.** Disease signatures from public data can be
  matched against drug-induced expression profiles to find compounds
  that reverse them ([Lamb et al. 2006](#ref-lamb2006cmap); [Sirota et
  al. 2011](#ref-sirota2011drug); [Subramanian et al.
  2017](#ref-subramanian2017l1000)). See
  [`vignette("downstream_biology")`](https://mdjubayerhossain.com/metaXpress/articles/downstream_biology.md).
- **Neurodegeneration and other complex diseases.** Meta-analysis
  distinguishes consistent changes from cohort- and region-specific ones
  ([Heberle et al. 2024](#ref-heberle2024ad)).
- **Reuse of public data.** Uniformly processed resources make it
  possible to reanalyse thousands of existing studies ([Collado-Torres
  et al. 2017](#ref-collado2017recount2); [Lachmann et al.
  2018](#ref-lachmann2018archs4); [Clough et al.
  2024](#ref-clough2024geo)).
- **Clinical translation.** Robust expression signatures are a step
  towards RNA-based diagnostics ([Byron et al.
  2016](#ref-byron2016translating)).

## 6. A checklist for a credible meta-analysis

![\*\*Figure 7. The analysis workflow and the metaXpress functions for
each step.\*\*
Schematic.](introduction_files/figure-html/checklist-1.png)

**Figure 7. The analysis workflow and the metaXpress functions for each
step.** Schematic.

The following practices come from the gene-expression meta-analysis
literature ([Ramasamy et al. 2008](#ref-ramasamy2008key); [Tseng et al.
2012](#ref-tseng2012review); [Sweeney et al.
2017](#ref-sweeney2017methods)) and from general meta-analysis
methodology.

1.  **Choose the question and the method before seeing results.** Random
    effects estimates an average effect. maxP or rOP ask whether a gene
    changes in all or most studies. Picking whichever method gives the
    most genes inflates false findings.
2.  **Make sure studies are independent.** Two GEO series from the same
    lab or brain bank may contain the same donors, which would count
    them twice. Check donor and case identifiers.
3.  **Harmonise the comparison.** Every study must compare the same
    groups in the same direction (case / control). Check with
    fold-change concordance.
4.  **Start from raw counts.** Run DE within each study, adjusting for
    known covariates such as sex or batch. Do not feed batch-corrected
    values into count models ([Nygaard et al.
    2016](#ref-nygaard2016batch)).
5.  **Report more than p-values.** Report:
    - the pooled effect with its confidence interval;
    - $`I^2`$ and $`\tau^2`$;
    - the prediction interval;
    - which genes survive leave-one-out analysis.
6.  **Be cautious with few studies.** With fewer than about five
    studies, $`\tau^2`$ is imprecise. Small-sample adjustments such as
    HKSJ ([IntHout et al. 2014](#ref-inthout2014hksj)) give
    better-calibrated confidence intervals, but can make genome-wide
    significance impossible. Tests for publication bias (e.g. Egger’s
    test ([Egger et al. 1997](#ref-egger1997bias))) need about ten
    studies ([Sterne et al. 2011](#ref-sterne2011recommendations)).
7.  **Interpret pathways at the right level.** Gene-set enrichment
    ([Subramanian et al. 2005](#ref-subramanian2005gsea); [Korotkevich
    et al. 2016](#ref-korotkevich2016fgsea); [Liberzon et al.
    2015](#ref-liberzon2015hallmark)) can be run on the meta-signature,
    or within each study and then combined ([Shen and Tseng
    2010](#ref-shen2010mape)).
8.  **Remember tissue composition.** In bulk data, a change in cell-type
    proportions looks like a change in gene expression. Check cell-type
    marker scores ([Becht et al. 2016](#ref-becht2016mcpcounter)).

## 7. Glossary

- **log2 fold change (log2FC)**: Difference in expression between groups
  on a log2 scale. +1 means twice as high in cases; −1 means half.

- **Standard error (SE)**: Uncertainty of a study’s log2FC. Smaller SE
  means more weight in the meta-analysis.

- **FDR / adjusted p-value**:

  Expected proportion of false discoveries among the genes called
  significant ([Benjamini and Hochberg 1995](#ref-benjamini1995fdr)).

- **Fixed effect**: Model assuming one true effect shared by all
  studies.

- **Random effects**: Model assuming true effects vary between studies
  around an average.

- **$`\tau^2`$ (tau-squared)**: Between-study variance of the true
  effects.

- **$`I^2`$**: Share of total variation due to between-study differences
  (0–100%).

- **Prediction interval**: Expected range of the effect in a new study.

- **Leave-one-out**: Re-running the meta-analysis without each study in
  turn.

- **Concordance**: Agreement of fold changes between two studies.

## 8. Where to go next

1.  [`vignette("quickstart")`](https://mdjubayerhossain.com/metaXpress/articles/quickstart.md):
    run your first meta-analysis in five minutes.
2.  [`vignette("case_study")`](https://mdjubayerhossain.com/metaXpress/articles/case_study.md):
    the full AD analysis, step by step, with every decision explained.
3.  [`vignette("full_workflow")`](https://mdjubayerhossain.com/metaXpress/articles/full_workflow.md):
    from GEO accession numbers to an HTML report.
4.  [`vignette("downstream_biology")`](https://mdjubayerhossain.com/metaXpress/articles/downstream_biology.md):
    transcription-factor activity, cell composition, co-expression
    modules and drug reversal.

## References

Becht, Etienne, Nicolas A. Giraldo, Laetitia Lacroix, et al. 2016.
“Estimating the Population Abundance of Tissue-Infiltrating Immune and
Stromal Cell Populations Using Gene Expression.” *Genome Biology* 17
(1): 218. <https://doi.org/10.1186/s13059-016-1070-5>.

Benjamini, Yoav, and Yosef Hochberg. 1995. “Controlling the False
Discovery Rate: A Practical and Powerful Approach to Multiple Testing.”
*Journal of the Royal Statistical Society: Series B (Methodological)* 57
(1): 289–300. <https://doi.org/10.1111/j.2517-6161.1995.tb02031.x>.

Borenstein, Michael, Larry V. Hedges, Julian P. T. Higgins, et al. 2010.
“A Basic Introduction to Fixed-Effect and Random-Effects Models for
Meta-Analysis.” *Research Synthesis Methods* 1 (2): 97–111.
<https://doi.org/10.1002/jrsm.12>.

Byron, Sara A., Kendall R. Van Keuren-Jensen, David M. Engelthaler, et
al. 2016. “Translating RNA Sequencing into Clinical Diagnostics:
Opportunities and Challenges.” *Nature Reviews Genetics* 17 (5): 257–71.
<https://doi.org/10.1038/nrg.2016.10>.

Cancer Genome Atlas Research Network, John N. Weinstein, Eric A.
Collisson, et al. 2013. “The Cancer Genome Atlas Pan-Cancer Analysis
Project.” *Nature Genetics* 45 (10): 1113–20.
<https://doi.org/10.1038/ng.2764>.

Clough, Emily, Tanya Barrett, Stephen E. Wilhite, et al. 2024. “NCBI
GEO: Archive for Gene Expression and Epigenomics Data Sets: 23-Year
Update.” *Nucleic Acids Research* 52 (D1): D138–44.
<https://doi.org/10.1093/nar/gkad965>.

Collado-Torres, Leonardo, Abhinav Nellore, Kai Kammers, et al. 2017.
“Reproducible RNA-seq Analysis Using recount2.” *Nature Biotechnology*
35 (4): 319–21. <https://doi.org/10.1038/nbt.3838>.

Conesa, Ana, Pedro Madrigal, Sonia Tarazona, et al. 2016. “A Survey of
Best Practices for RNA-seq Data Analysis.” *Genome Biology* 17 (1): 13.
<https://doi.org/10.1186/s13059-016-0881-8>.

DerSimonian, Rebecca, and Nan Laird. 1986. “Meta-Analysis in Clinical
Trials.” *Controlled Clinical Trials* 7 (3): 177–88.
<https://doi.org/10.1016/0197-2456(86)90046-2>.

Dobin, Alexander, Carrie A. Davis, Felix Schlesinger, et al. 2013.
“STAR: Ultrafast Universal RNA-seq Aligner.” *Bioinformatics* 29 (1):
15–21. <https://doi.org/10.1093/bioinformatics/bts635>.

Egger, Matthias, George Davey Smith, Martin Schneider, et al. 1997.
“Bias in Meta-Analysis Detected by a Simple, Graphical Test.” *BMJ* 315
(7109): 629–34. <https://doi.org/10.1136/bmj.315.7109.629>.

GTEx Consortium. 2013. “The Genotype-Tissue Expression (GTEx) Project.”
*Nature Genetics* 45 (6): 580–85. <https://doi.org/10.1038/ng.2653>.

Haque, Ashraful, Jessica Engel, Sarah A. Teichmann, et al. 2017. “A
Practical Guide to Single-Cell RNA-Sequencing for Biomedical Research
and Clinical Applications.” *Genome Medicine* 9 (1): 75.
<https://doi.org/10.1186/s13073-017-0467-4>.

Heberle, Bernardo Aguzzoli, Kristin L. Fox, Lucas Lobraico Libermann, et
al. 2024. *Systematic Review and Meta-Analysis of Bulk RNAseq Studies in
Human Alzheimer’s Disease Brain Tissue*. bioRxiv preprint
2024.11.07.622520. <https://doi.org/10.1101/2024.11.07.622520>.

Higgins, Julian P. T., and Simon G. Thompson. 2002. “Quantifying
Heterogeneity in a Meta-Analysis.” *Statistics in Medicine* 21 (11):
1539–58. <https://doi.org/10.1002/sim.1186>.

Higgins, Julian P. T., Simon G. Thompson, and David J. Spiegelhalter.
2009. “A Re-Evaluation of Random-Effects Meta-Analysis.” *Journal of the
Royal Statistical Society: Series A (Statistics in Society)* 172 (1):
137–59. <https://doi.org/10.1111/j.1467-985X.2008.00552.x>.

Hu, Yanshi, Zixuan Wang, Yueming Hu, et al. 2026. “AWmeta Empowers
Adaptively Weighted Transcriptomic Meta-Analysis.” *Current Issues in
Molecular Biology* 48 (5): 530. <https://doi.org/10.3390/cimb48050530>.

IntHout, Joanna, John P. A. Ioannidis, and George F. Borm. 2014. “The
Hartung-Knapp-Sidik-Jonkman Method for Random Effects Meta-Analysis Is
Straightforward and Considerably Outperforms the Standard
DerSimonian-Laird Method.” *BMC Medical Research Methodology* 14: 25.
<https://doi.org/10.1186/1471-2288-14-25>.

Ioannidis, John P. A. 2005. “Why Most Published Research Findings Are
False.” *PLOS Medicine* 2 (8): e124.
<https://doi.org/10.1371/journal.pmed.0020124>.

Khatri, Purvesh, Silke Roedder, Naoyuki Kimura, et al. 2013. “A Common
Rejection Module (CRM) for Acute Rejection Across Multiple Organs
Identifies Novel Therapeutics for Organ Transplantation.” *Journal of
Experimental Medicine* 210 (11): 2205–21.
<https://doi.org/10.1084/jem.20122709>.

Kim, Daehwan, Joseph M. Paggi, Chanhee Park, et al. 2019. “Graph-Based
Genome Alignment and Genotyping with HISAT2 and HISAT-genotype.” *Nature
Biotechnology* 37 (8): 907–15.
<https://doi.org/10.1038/s41587-019-0201-4>.

Korotkevich, Gennady, Vladimir Sukhov, Nikolay Budin, et al. 2016. *Fast
Gene Set Enrichment Analysis*. bioRxiv preprint 060012.
<https://doi.org/10.1101/060012>.

Lachmann, Alexander, Denis Torre, Alexandra B. Keenan, et al. 2018.
“Massive Mining of Publicly Available RNA-seq Data from Human and
Mouse.” *Nature Communications* 9 (1): 1366.
<https://doi.org/10.1038/s41467-018-03751-6>.

Lamb, Justin, Emily D. Crawford, David Peck, et al. 2006. “The
Connectivity Map: Using Gene-Expression Signatures to Connect Small
Molecules, Genes, and Disease.” *Science* 313 (5795): 1929–35.
<https://doi.org/10.1126/science.1132939>.

Langan, Dean, Julian P. T. Higgins, Dan Jackson, et al. 2019. “A
Comparison of Heterogeneity Variance Estimators in Simulated
Random-Effects Meta-Analyses.” *Research Synthesis Methods* 10 (1):
83–98. <https://doi.org/10.1002/jrsm.1316>.

Law, Charity W., Yunshun Chen, Wei Shi, et al. 2014. “Voom: Precision
Weights Unlock Linear Model Analysis Tools for RNA-seq Read Counts.”
*Genome Biology* 15 (2): R29.
<https://doi.org/10.1186/gb-2014-15-2-r29>.

Leek, Jeffrey T., Robert B. Scharpf, Héctor Corrada Bravo, et al. 2010.
“Tackling the Widespread and Critical Impact of Batch Effects in
High-Throughput Data.” *Nature Reviews Genetics* 11 (10): 733–39.
<https://doi.org/10.1038/nrg2825>.

Li, Jia, and George C. Tseng. 2011. “An Adaptively Weighted Statistic
for Detecting Differential Gene Expression When Combining Multiple
Transcriptomic Studies.” *The Annals of Applied Statistics* 5 (2A):
994–1019. <https://doi.org/10.1214/10-AOAS393>.

Liao, Yang, Gordon K. Smyth, and Wei Shi. 2014. “featureCounts: An
Efficient General Purpose Program for Assigning Sequence Reads to
Genomic Features.” *Bioinformatics* 30 (7): 923–30.
<https://doi.org/10.1093/bioinformatics/btt656>.

Liberzon, Arthur, Chet Birger, Helga Thorvaldsdóttir, et al. 2015. “The
Molecular Signatures Database Hallmark Gene Set Collection.” *Cell
Systems* 1 (6): 417–25. <https://doi.org/10.1016/j.cels.2015.12.004>.

Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. “Moderated
Estimation of Fold Change and Dispersion for RNA-seq Data with DESeq2.”
*Genome Biology* 15 (12): 550.
<https://doi.org/10.1186/s13059-014-0550-8>.

Mortazavi, Ali, Brian A Williams, Kenneth McCue, et al. 2008. “Mapping
and Quantifying Mammalian Transcriptomes by RNA-Seq.” *Nature Methods* 5
(7): 621–28. <https://doi.org/10.1038/nmeth.1226>.

Nygaard, Vegard, Einar Andreas Rødland, and Eivind Hovig. 2016. “Methods
That Remove Batch Effects While Retaining Group Differences May Lead to
Exaggerated Confidence in Downstream Analyses.” *Biostatistics* 17 (1):
29–39. <https://doi.org/10.1093/biostatistics/kxv027>.

Patro, Rob, Geet Duggal, Michael I Love, et al. 2017. “Salmon Provides
Fast and Bias-Aware Quantification of Transcript Expression.” *Nature
Methods* 14 (4): 417–19. <https://doi.org/10.1038/nmeth.4197>.

Ramasamy, Adaikalavan, Adrian Mondry, Chris C Holmes, et al. 2008. “Key
Issues in Conducting a Meta-Analysis of Gene Expression Microarray
Datasets.” *PLOS Medicine* 5 (9): e184.
<https://doi.org/10.1371/journal.pmed.0050184>.

Rau, Andrea, Guillemette Marot, and Florence Jaffrézic. 2014.
“Differential Meta-Analysis of RNA-seq Data from Multiple Studies.” *BMC
Bioinformatics* 15: 91. <https://doi.org/10.1186/1471-2105-15-91>.

Ritchie, Matthew E., Belinda Phipson, Di Wu, et al. 2015. “Limma Powers
Differential Expression Analyses for RNA-Sequencing and Microarray
Studies.” *Nucleic Acids Research* 43 (7): e47.
<https://doi.org/10.1093/nar/gkv007>.

Robinson, Mark D., Davis J. McCarthy, and Gordon K. Smyth. 2010. “edgeR:
A Bioconductor Package for Differential Expression Analysis of Digital
Gene Expression Data.” *Bioinformatics* 26 (1): 139–40.
<https://doi.org/10.1093/bioinformatics/btp616>.

Schurch, Nicholas J., Pietá Schofield, Marek Gierliński, et al. 2016.
“How Many Biological Replicates Are Needed in an RNA-seq Experiment and
Which Differential Expression Tool Should You Use?” *RNA* 22 (6):
839–51. <https://doi.org/10.1261/rna.053959.115>.

Shen, Kui, and George C. Tseng. 2010. “Meta-Analysis for Pathway
Enrichment Analysis When Combining Multiple Genomic Studies.”
*Bioinformatics* 26 (10): 1316–23.
<https://doi.org/10.1093/bioinformatics/btq148>.

Sirota, Marina, Joel T. Dudley, Jeewon Kim, et al. 2011. “Discovery and
Preclinical Validation of Drug Indications Using Compendia of Public
Gene Expression Data.” *Science Translational Medicine* 3 (96): 96ra77.
<https://doi.org/10.1126/scitranslmed.3001318>.

Song, Chi, and George C. Tseng. 2014. “Hypothesis Setting and Order
Statistic for Robust Genomic Meta-Analysis.” *The Annals of Applied
Statistics* 8 (2): 777–800. <https://doi.org/10.1214/13-AOAS683>.

Stark, Rory, Marta Grzelak, and James Hadfield. 2019. “RNA Sequencing:
The Teenage Years.” *Nature Reviews Genetics* 20 (11): 631–56.
<https://doi.org/10.1038/s41576-019-0150-2>.

Sterne, Jonathan A. C., Alex J. Sutton, John P. A. Ioannidis, et al.
2011. “Recommendations for Examining and Interpreting Funnel Plot
Asymmetry in Meta-Analyses of Randomised Controlled Trials.” *BMJ* 343:
d4002. <https://doi.org/10.1136/bmj.d4002>.

Subramanian, Aravind, Rajiv Narayan, Steven M. Corsello, et al. 2017. “A
Next Generation Connectivity Map: L1000 Platform and the First 1,000,000
Profiles.” *Cell* 171 (6): 1437–1452.e17.
<https://doi.org/10.1016/j.cell.2017.10.049>.

Subramanian, Aravind, Pablo Tamayo, Vamsi K. Mootha, et al. 2005. “Gene
Set Enrichment Analysis: A Knowledge-Based Approach for Interpreting
Genome-Wide Expression Profiles.” *Proceedings of the National Academy
of Sciences of the United States of America* 102 (43): 15545–50.
<https://doi.org/10.1073/pnas.0506580102>.

Sweeney, Timothy E, Lindsay Braviak, Cristina M Tato, et al. 2016.
“Genome-Wide Expression for Diagnosis of Pulmonary Tuberculosis: A
Multicohort Analysis.” *The Lancet Respiratory Medicine* 4 (3): 213–24.
<https://doi.org/10.1016/S2213-2600(16)00048-5>.

Sweeney, Timothy E., Winston A. Haynes, Francesco Vallania, et al. 2017.
“Methods to Increase Reproducibility in Differential Gene Expression via
Meta-Analysis.” *Nucleic Acids Research* 45 (1): e1.
<https://doi.org/10.1093/nar/gkw797>.

Sweeney, Timothy E., Aaditya Shidham, Hector R. Wong, et al. 2015. “A
Comprehensive Time-Course–Based Multicohort Analysis of Sepsis and
Sterile Inflammation Reveals a Robust Diagnostic Gene Set.” *Science
Translational Medicine* 7 (287): 287ra71.
<https://doi.org/10.1126/scitranslmed.aaa5993>.

Tseng, George C., Debashis Ghosh, and Eleanor Feingold. 2012.
“Comprehensive Literature Review and Statistical Considerations for
Microarray Meta-Analysis.” *Nucleic Acids Research* 40 (9): 3785–99.
<https://doi.org/10.1093/nar/gkr1265>.

Viechtbauer, Wolfgang. 2005. “Bias and Efficiency of Meta-Analytic
Variance Estimators in the Random-Effects Model.” *Journal of
Educational and Behavioral Statistics* 30 (3): 261–93.
<https://doi.org/10.3102/10769986030003261>.

Wang, Zhong, Mark Gerstein, and Michael Snyder. 2009. “RNA-Seq: A
Revolutionary Tool for Transcriptomics.” *Nature Reviews Genetics* 10
(1): 57–63. <https://doi.org/10.1038/nrg2484>.

Wilks, Christopher, Shijie C. Zheng, Feng Yong Chen, et al. 2021.
“recount3: Summaries and Queries for Large-Scale RNA-seq Expression and
Splicing.” *Genome Biology* 22 (1): 323.
<https://doi.org/10.1186/s13059-021-02533-6>.
