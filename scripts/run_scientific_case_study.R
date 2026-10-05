# ==============================================================================
# metaXpress Scientific Case Study: Multi-Cohort Solid Cancer Meta-Analysis
# ==============================================================================
# This script performs a complete, end-to-end scientific case study using three
# publicly available GEO bulk RNA-seq cohorts comparing cancer vs normal tissues:
#   1. GSE130688: Colorectal Adenocarcinoma vs Normal Colon (n = 30)
#   2. GSE136569: Pancreatic Ductal Adenocarcinoma vs Normal Pancreas (n = 10)
#   3. GSE171485: Clear Cell Renal Cell Carcinoma vs Normal Kidney (n = 12)
# Total: 52 patient samples across 3 distinct epithelial malignancies.
# ==============================================================================

suppressPackageStartupMessages({
  library(metaXpress)
  library(ggplot2)
})

out_dir <- "case_study_results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message("\n==================================================================")
message("🔬 STEP 1: Ingesting Public GEO Studies & 10-Point QC Scoring")
message("==================================================================")

studies <- mx_load_local(
  count_paths = c(
    "data/GSE130688_raw_counts.tsv",
    "data/GSE136569_raw_counts.tsv",
    "data/GSE171485_raw_counts.tsv"
  ),
  metadata_paths = c(
    "data/GSE130688_metadata.csv",
    "data/GSE136569_metadata.csv",
    "data/GSE171485_metadata.csv"
  ),
  organism = "Homo sapiens"
)

# Rename to clean accession identifiers
names(studies) <- c("GSE130688_CRC", "GSE136569_PDAC", "GSE171485_ccRCC")

for (s_name in names(studies)) {
  s <- studies[[s_name]]
  message(sprintf("  [%s] Samples: %d (%s), Genes: %d, QC Score: %d/10",
                  s_name, ncol(s@counts),
                  paste(table(s@metadata$condition), collapse = " vs "),
                  nrow(s@counts), s@qc_score))
}

# Generate study overview plot
overview_plots <- mx_study_overview(studies)
p_grid <- patchwork::wrap_plots(overview_plots$qc_scores,
                                overview_plots$sample_sizes,
                                overview_plots$library_sizes,
                                overview_plots$gene_detection, ncol = 2)
ggplot2::ggsave(file.path(out_dir, "case_study_study_overview.png"),
                p_grid, width = 9, height = 7, dpi = 300)
message("  Saved: case_study_study_overview.png")

# Filter studies
studies <- mx_filter_studies(studies, qc_threshold = 7)
message(sprintf("  Retained %d/%d studies passing QC threshold >= 7", length(studies), 3))

message("\n==================================================================")
message("🧬 STEP 2: Harmonization & Gene Space Alignment")
message("==================================================================")

studies <- mx_align_genes(studies)
common_genes_n <- nrow(studies[[1]]@counts)
message(sprintf("  Common genes across all 3 independent cohorts: %d", common_genes_n))

message("\n==================================================================")
message("📊 STEP 3: Per-Study Differential Expression (DESeq2 GLM)")
message("==================================================================")

studies <- mx_de_all(studies, method = "DESeq2")
de_summary <- mx_de_summary(studies)
print(de_summary)

message("\n==================================================================")
message("🔍 STEP 4: Missing Gene Analysis & Coverage Profiling")
message("==================================================================")

de_results <- lapply(studies, function(s) s@de_result)
cov_mat <- mx_missing_summary(de_results)
cov_pct <- attr(cov_mat, "coverage_pct")
message(sprintf("  100%% complete gene presence across cohorts: %d / %d (%.1f%%)",
                sum(cov_pct == 100), length(cov_pct),
                mean(cov_pct == 100) * 100))

message("\n==================================================================")
message("🧮 STEP 5: Multi-Study Meta-Analysis (Random Effects DerSimonian-Laird)")
message("==================================================================")

n_samples <- vapply(studies, function(s) ncol(s@counts), integer(1))

meta_res <- mx_meta(
  de_results = de_results,
  method = "random_effects",
  min_studies = 3,
  n_samples = n_samples
)

mt <- meta_res@meta_table

# ------------------------------------------------------------------------------
# STEP 6: Biological Annotation using NCBI GEO Entrez-to-Symbol mapping
# ------------------------------------------------------------------------------
message("\n==================================================================")
message("🏷️ STEP 6: Gene Annotation & Cancer Signature Discovery")
message("==================================================================")

if (file.exists("geo_cache/GSE53697/GSE53697_RNAseq_AD.txt.gz")) {
  geo_dict <- read.delim("geo_cache/GSE53697/GSE53697_RNAseq_AD.txt.gz")[, 1:2]
  colnames(geo_dict) <- c("entrez_id", "symbol")
  geo_dict$entrez_id <- as.character(geo_dict$entrez_id)
  
  # Merge symbols into meta table
  mt$symbol <- geo_dict$symbol[match(as.character(mt$gene_id), geo_dict$entrez_id)]
  mt$symbol[is.na(mt$symbol) | mt$symbol == ""] <- mt$gene_id[is.na(mt$symbol) | mt$symbol == ""]
} else {
  mt$symbol <- mt$gene_id
}

# Update meta_res table with symbols
meta_res@meta_table <- mt

# Identify significant Meta-DEGs
padj_cut <- 0.05
lfc_cut  <- 1.0

sig_meta <- mt[!is.na(mt$meta_padj) & mt$meta_padj <= padj_cut & abs(mt$meta_log2FC) >= lfc_cut, ]
sig_up   <- sig_meta[sig_meta$meta_log2FC > 0, ]
sig_down <- sig_meta[sig_meta$meta_log2FC < 0, ]

message(sprintf("  Total genes tested in meta-analysis: %d", nrow(mt)))
message(sprintf("  Statistically Significant Meta-DEGs (FDR < %.2f, |log2FC| >= %.1f): %d",
                padj_cut, lfc_cut, nrow(sig_meta)))
message(sprintf("    ▲ Shared Pan-Cancer Up-regulated: %d", nrow(sig_up)))
message(sprintf("    ▼ Shared Pan-Cancer Down-regulated: %d", nrow(sig_down)))

# Top 15 Up-regulated genes
sig_up_sorted <- sig_up[order(sig_up$meta_padj), ]
message("\n  Top 10 Conserved Up-regulated Pan-Cancer Genes:")
top_up_print <- head(sig_up_sorted[, c("gene_id", "symbol", "meta_log2FC", "meta_padj", "i_squared", "direction_consistency")], 10)
print(top_up_print)

# Top 15 Down-regulated genes
sig_down_sorted <- sig_down[order(sig_down$meta_padj), ]
message("\n  Top 10 Conserved Down-regulated Pan-Cancer Genes:")
top_down_print <- head(sig_down_sorted[, c("gene_id", "symbol", "meta_log2FC", "meta_padj", "i_squared", "direction_consistency")], 10)
print(top_down_print)

# ------------------------------------------------------------------------------
# STEP 7: Publication-Grade Visualizations
# ------------------------------------------------------------------------------
message("\n==================================================================")
message("📈 STEP 7: Generating Scientific Visualizations")
message("==================================================================")

# 7a. Volcano Plot with annotated Gene Symbols
# Temporarily replace gene_id with symbol for clear plot labels
meta_res_plot <- meta_res
meta_res_plot@meta_table$gene_id <- meta_res_plot@meta_table$symbol
p_volcano <- mx_volcano(meta_res_plot, padj_threshold = padj_cut,
                        lfc_threshold = lfc_cut, label_top = 12,
                        title = "Pan-Cancer vs Normal Meta-Analysis Volcano Plot (Random Effects)")
ggplot2::ggsave(file.path(out_dir, "case_study_volcano_plot.png"),
                p_volcano, width = 8, height = 6.5, dpi = 300)
message("  Saved: case_study_volcano_plot.png")

# 7b. Heterogeneity Plot (I-squared distribution)
p_het <- mx_heterogeneity_plot(meta_res) +
  ggplot2::labs(title = "Between-Study Heterogeneity (I² Index) Across 3 Cancer Cohorts",
                subtitle = "Classification according to Higgins et al. (2003)")
ggplot2::ggsave(file.path(out_dir, "case_study_heterogeneity_plot.png"),
                p_het, width = 7.5, height = 5.5, dpi = 300)
message("  Saved: case_study_heterogeneity_plot.png")

# 7c. Forest Plots for Top Pan-Cancer Markers
top_up_gene_id <- sig_up_sorted$gene_id[1]
top_up_symbol  <- sig_up_sorted$symbol[1]
p_forest_up <- mx_forest(top_up_gene_id, de_results, meta_res) +
  ggplot2::labs(title = paste0("Forest Plot: Conserved Up-regulation of ", top_up_symbol, " (", top_up_gene_id, ")"),
                subtitle = "Individual cohort effect sizes vs Pooled DerSimonian-Laird Random Effects model")
ggplot2::ggsave(file.path(out_dir, "case_study_forest_top_up.png"),
                p_forest_up, width = 7, height = 4.5, dpi = 300)
message(sprintf("  Saved: case_study_forest_top_up.png (%s)", top_up_symbol))

top_down_gene_id <- sig_down_sorted$gene_id[1]
top_down_symbol  <- sig_down_sorted$symbol[1]
p_forest_down <- mx_forest(top_down_gene_id, de_results, meta_res) +
  ggplot2::labs(title = paste0("Forest Plot: Conserved Down-regulation of ", top_down_symbol, " (", top_down_gene_id, ")"),
                subtitle = "Individual cohort effect sizes vs Pooled DerSimonian-Laird Random Effects model")
ggplot2::ggsave(file.path(out_dir, "case_study_forest_top_down.png"),
                p_forest_down, width = 7, height = 4.5, dpi = 300)
message(sprintf("  Saved: case_study_forest_top_down.png (%s)", top_down_symbol))

# ------------------------------------------------------------------------------
# STEP 8: Export Results & Summary Artifacts
# ------------------------------------------------------------------------------
message("\n==================================================================")
message("💾 STEP 8: Exporting Scientific Results")
message("==================================================================")

mx_export(meta_res, format = "csv", output_dir = out_dir, prefix = "pan_cancer_meta_results")
message("  Saved: pan_cancer_meta_results.csv")

# Save summary metrics RDS
metrics <- list(
  studies = names(studies),
  total_samples = sum(n_samples),
  common_genes = common_genes_n,
  de_summary = de_summary,
  meta_method = "random_effects",
  total_tested = nrow(mt),
  n_sig_meta = nrow(sig_meta),
  n_up = nrow(sig_up),
  n_down = nrow(sig_down),
  top_up = top_up_print,
  top_down = top_down_print,
  median_I2 = median(mt$i_squared, na.rm = TRUE)
)
saveRDS(metrics, file.path(out_dir, "case_study_metrics.rds"))
message("  Saved: case_study_metrics.rds")

message("\n==================================================================")
message("🎉 CASE STUDY EXECUTION COMPLETED SUCCESSFULLY!")
message("==================================================================\n")
