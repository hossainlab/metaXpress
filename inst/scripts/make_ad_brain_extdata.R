# ==============================================================================
# Build inst/extdata/ad_brain_de.rds
# ==============================================================================
# Alzheimer's disease (AD) vs control, human postmortem brain, bulk RNA-seq.
# Four independent GEO series, all quantified by the NCBI GEO RNA-seq pipeline
# (HISAT2 -> GRCh38.p13, featureCounts on NCBI RefSeq genes), so every study
# shares one alignment/annotation/counting protocol.
#
#   GSE125583  fusiform gyrus         AD 219 / control 70  (Srinivasan 2020)
#   GSE104704  lateral temporal lobe  AD  12 / aged ctrl 10 (Nativio 2018)
#   GSE53697   dorsolateral PFC       AD   6 / control  6 (Magistri 2015;
#              NCBI processed 12 of the 17 deposited libraries)
#   GSE173955  hippocampus            AD   8 / control 10
#
# Excluded on purpose:
#   GSE95587   89/117 donors overlap GSE125583 (same brain bank and lab)
#   GSE159699  same CNDR donors as GSE104704
# Including either would count the same people twice and violate the
# independence assumption of meta-analysis.
#
# Usage (from the package root, ~15 min, ~150 MB download):
#   Rscript inst/scripts/make_ad_brain_extdata.R [cache_dir]
# ==============================================================================

suppressPackageStartupMessages({
  library(metaXpress)
  library(SummarizedExperiment)
})

args      <- commandArgs(trailingOnly = TRUE)
cache_dir <- if (length(args)) args[1] else "geo_cache/ad_brain"
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = 1800)

design <- list(
  GSE125583 = list(col = "diagnosis:ch1",
                   contrast = c("Alzheimer's disease", "control"),
                   formula = ~ sex + condition, region = "Fusiform gyrus"),
  GSE104704 = list(col = "study group:ch1",
                   contrast = c("Aged, diseased", "Old"),
                   formula = ~ condition, region = "Lateral temporal lobe"),
  GSE53697  = list(col = "disease status:ch1",
                   contrast = c("advanced Alzheimer's Disease", "control"),
                   formula = ~ condition, region = "Dorsolateral PFC"),
  GSE173955 = list(col = "title",
                   contrast = c("Alzheimer's Disease", "non-Alzheimer's Disease"),
                   formula = ~ condition, region = "Hippocampus")
)

studies <- mx_fetch_geo(names(design),
                        condition_col = vapply(design, `[[`, "", "col"),
                        cache_dir = cache_dir)

prep <- function(s, d) {
  if (s@accession == "GSE173955")
    s$condition <- sub("_Biological.*$", "", s$condition)
  if (s@accession == "GSE125583") {
    s$sex <- s$`Sex:ch1`
    s <- s[, !is.na(s$sex) & s$sex %in% c("F", "M")]
  }
  s <- s[, s$condition %in% d$contrast]

  # Protein-coding genes, gene symbols as identifiers
  rd <- rowData(s)
  s  <- s[!is.na(rd$symbol) & rd$gene_type == "protein-coding" &
            !duplicated(rd$symbol), ]
  rownames(s) <- rowData(s)$symbol

  # Expression filter: >= 10 counts in at least the smaller group size
  n_min <- min(table(s$condition))
  s[rowSums(assay(s, "counts") >= 10) >= n_min, ]
}
studies <- Map(prep, studies, design[names(studies)])

studies <- mx_de_all(studies, method = "DESeq2",
                     formula  = lapply(design, `[[`, "formula"),
                     contrast = lapply(design, `[[`, "contrast"))

de_results <- lapply(studies, function(s) {
  d <- s@de_result
  d <- d[, c("gene_id", "log2FC", "lfcSE", "stat", "pvalue", "padj",
             "baseMean", "n_case", "n_control")]
  # 6 significant digits keeps the file small without affecting inference
  for (v in c("log2FC", "lfcSE", "stat", "pvalue", "padj", "baseMean"))
    d[[v]] <- signif(d[[v]], 6)
  d
})

summary_tbl <- data.frame(
  accession = names(studies),
  region    = vapply(design, `[[`, "", "region")[names(studies)],
  n_AD      = vapply(de_results, function(d) d$n_case[1], integer(1)),
  n_control = vapply(de_results, function(d) d$n_control[1], integer(1)),
  n_genes   = vapply(de_results, nrow, integer(1)),
  n_DEG_fdr05 = vapply(de_results, function(d) sum(d$padj <= 0.05, na.rm = TRUE),
                       integer(1)),
  formula   = vapply(design, function(x) deparse(x$formula), "")[names(studies)],
  row.names = NULL
)
print(summary_tbl)

ad_brain <- list(
  de_results = de_results,
  studies    = summary_tbl,
  provenance = list(
    built      = format(Sys.Date()),
    counts     = "NCBI GEO RNA-seq pipeline, GRCh38.p13, NCBI RefSeq genes",
    de         = paste("DESeq2", packageVersion("DESeq2"),
                       "Wald test, unshrunken MLE log2FC"),
    filter     = "protein-coding; >= 10 counts in >= smaller-group-size samples",
    metaXpress = as.character(packageVersion("metaXpress"))
  )
)
saveRDS(ad_brain, "inst/extdata/ad_brain_de.rds", compress = "xz")
message("Saved inst/extdata/ad_brain_de.rds (",
        round(file.size("inst/extdata/ad_brain_de.rds") / 1024), " KB)")
