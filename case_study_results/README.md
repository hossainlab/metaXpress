# PDAC downstream-biology case study: outputs

Outputs of `scripts/run_case_study_applications.R` for a meta-analysis of
**pancreatic ductal adenocarcinoma (PDAC) vs adjacent non-tumour pancreas**
across three independent GEO cohorts:

| GEO | Tumour / non-tumour |
|---|---|
| GSE130688 | 15 / 15 |
| GSE136569 | 5 / 5 |
| GSE171485 | 6 / 6 |

All three cohorts are PDAC according to their GEO series records. An
earlier version of this case study labelled GSE130688 as colorectal and
GSE171485 as renal carcinoma, which was incorrect.

## Pipeline

1. **Counts.** NCBI GEO-generated raw counts (HISAT2 / featureCounts,
   GRCh38.p13).
2. **Labels.** Tumour/normal labels taken from the GEO series metadata.
3. **Per-study DE.** DESeq2 (tumour vs normal) in each cohort.
4. **Meta-analysis.** Random effects with REML tau² over 29,679 genes.
   1,787 meta-DEGs at FDR < 0.05 and |log2FC| >= 1 (1,166 up, 621 down).
5. **Module 9.** The downstream analyses below.

Rebuild with:

```sh
Rscript scripts/fetch_case_study_data.R
Rscript scripts/run_case_study_applications.R
Rscript scripts/make_downstream_extdata.R   # refresh inst/extdata demo file
```

## Files

| File | Content |
|---|---|
| `app_tf_activity.*` | CollecTRI TF activity (ULM): 552 TFs scored, 148 at FDR < 0.05 |
| `app_cell_composition.*` | MCP-counter scores, Hedges' g pooled across cohorts |
| `app_coexpression_modules.*`, `app_coexpression_module_stats.csv` | Consensus co-expression modules (3,000 genes) |
| `app_drug_connectivity.*` | LINCS L1000 reversal ranking (5,425 compounds) |
| `app_downstream_panel.*` | Combined figure |
| `app_metrics.rds` | Summary objects and session info |

## Main findings

**Transcription factors.**
- Activated: NF-kB (RELA, REL, RELB) and AP-1 (JUN) at the top, followed by
  E2F1/E2F4, MYC and TP53.
- Repressed: SNAI2, HDAC1 and HDAC4.

**Cell composition.**
- Neutrophils, monocytic cells, myeloid dendritic cells and fibroblasts are
  significantly higher in tumour (FDR < 0.05, I² <= 27%).
- This is the desmoplastic, myeloid-rich PDAC stroma. Part of the
  gene-level signature therefore reflects tissue composition.

**Co-expression modules.**

| Module | Hubs | Direction | FDR |
|---|---|---|---|
| M2 | MST1R, PLEK2, TSPAN1, B3GNT3, CDCP1 | up | 5.8e-5 |
| M7 | ONECUT2, USH1C, SLC5A1 | up | 1.0e-3 |
| M1 | STAT3, IL4R | up | 0.01 |
| M3 | CELP, PNLIPRP1, RBPJL | down (acinar loss) | 0.06 |

**Drug reversal.**
- Five of the top eight reversers inhibit the TGF-beta receptor ALK5:
  SD-208, GW-788388, LY-2109761, LY-364947 and SB-525334.
- The EGFR inhibitor erlotinib, approved for advanced PDAC, ranks fifth.

**Caveat.** Drug-connectivity p-values assume independent genes and are
anti-conservative. Use the scores for ranking only.
