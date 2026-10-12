# ==============================================================================
# metaXpress Case Study, Part 2: Downstream Biology Applications
# ==============================================================================
# Rebuilds the pan-cancer meta-analysis (CRC, PDAC, ccRCC vs normal tissue)
# from NCBI GEO counts, then applies Module 9:
#   A. Transcription factor activity      (CollecTRI regulons)
#   B. Cell-type composition shifts       (MCP-counter markers)
#   C. Consensus co-expression modules
#   D. Drug signature reversal            (LINCS L1000 consensus signatures)
#
# Run scripts/fetch_case_study_data.R first to download the inputs.
# Outputs go to case_study_results/.
# ==============================================================================

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(ggplot2)
})

in_dir  <- "geo_cache/case_study"
out_dir <- "case_study_results"
dir.create(out_dir, showWarnings = FALSE)

save_plot <- function(p, name, width = 7, height = 6) {
  for (ext in c("png", "pdf"))
    ggsave(file.path(out_dir, paste0(name, ".", ext)), p,
           width = width, height = height, dpi = 300)
  message("  Saved: ", name)
}

# ------------------------------------------------------------------------------
# 1. Rebuild studies from NCBI counts + curated condition labels
# ------------------------------------------------------------------------------
message("\n== 1. Building studies from NCBI GEO counts ==")

labels <- readRDS("inst/extdata/case_study_subset.rds")
cohorts <- c(GSE130688 = "CRC", GSE136569 = "PDAC", GSE171485 = "ccRCC")

studies <- lapply(names(cohorts), function(acc) {
  counts <- read.delim(file.path(in_dir, paste0(acc, "_raw_counts.tsv.gz")),
                       row.names = 1, check.names = FALSE)
  meta <- as.data.frame(SummarizedExperiment::colData(labels[[acc]]))
  counts <- as.matrix(counts[, meta$sample_id])
  meta$condition <- ifelse(meta$condition == "normal", "normal", "tumor")
  mx_qc_study(metaXpressStudy(counts, meta, accession = acc))
})
names(studies) <- paste(names(cohorts), cohorts, sep = "_")

for (nm in names(studies))
  message(sprintf("  %-16s %2d samples (%s), %d genes, QC %d/10", nm,
                  ncol(studies[[nm]]),
                  paste(table(studies[[nm]]$condition), collapse = "/"),
                  nrow(studies[[nm]]), studies[[nm]]@qc_score))

# Entrez ID -> HGNC symbol (NCBI gene_info)
gene_info <- read.delim(file.path(in_dir, "Homo_sapiens.gene_info.gz"),
                        quote = "", colClasses = "character")
symbol_of <- setNames(gene_info$Symbol, gene_info$GeneID)

# ------------------------------------------------------------------------------
# 2. Per-study DE and random-effects meta-analysis
# ------------------------------------------------------------------------------
message("\n== 2. Per-study DESeq2 + random-effects meta-analysis ==")

studies <- mx_filter_studies(studies, qc_threshold = 7)
studies <- mx_align_genes(studies)
studies <- mx_de_all(studies, method = "DESeq2")
de_results <- lapply(studies, function(s) s@de_result)

meta_res <- mx_meta(de_results, method = "random_effects", min_studies = 3,
                    n_samples = vapply(studies, ncol, integer(1)))
mt <- meta_res@meta_table
mt$symbol <- unname(symbol_of[mt$gene_id])
meta_res@meta_table <- mt

# Cache pipeline objects for scripts/make_downstream_extdata.R
saveRDS(list(studies = studies, meta_res = meta_res),
        file.path(in_dir, "pipeline_objects.rds"))

sig <- !is.na(mt$meta_padj) & mt$meta_padj <= 0.05 & abs(mt$meta_log2FC) >= 1
message(sprintf("  %d genes tested, %d meta-DEGs (%d up, %d down)",
                nrow(mt), sum(sig), sum(sig & mt$meta_log2FC > 0),
                sum(sig & mt$meta_log2FC < 0)))

# ------------------------------------------------------------------------------
# A. Transcription factor activity
# ------------------------------------------------------------------------------
message("\n== A. TF activity (CollecTRI, ULM) ==")

collectri <- read.delim(file.path(in_dir, "collectri.tsv"))
tf <- mx_tf_activity(meta_res, collectri, id_col = "symbol", min_size = 10)
message(sprintf("  %d TFs scored, %d with padj <= 0.05", nrow(tf),
                sum(tf$padj <= 0.05)))
print(head(tf[order(-tf$score), ], 8), row.names = FALSE)
print(head(tf[order(tf$score), ], 8), row.names = FALSE)
write.csv(tf, file.path(out_dir, "app_tf_activity.csv"), row.names = FALSE)
save_plot(mx_downstream_plot(tf, top = 30), "app_tf_activity", height = 7)

# ------------------------------------------------------------------------------
# B. Cell-type composition
# ------------------------------------------------------------------------------
message("\n== B. Cell-type composition (MCP-counter markers) ==")

mcp <- read.delim(file.path(in_dir, "mcpcounter_genes.txt"),
                  colClasses = "character", check.names = FALSE)
markers <- data.frame(gene = mcp$ENTREZID, cell_type = mcp$`Cell population`)
comp <- mx_cell_composition(studies, markers)
print(comp, row.names = FALSE, digits = 3)
write.csv(comp, file.path(out_dir, "app_cell_composition.csv"),
          row.names = FALSE)
save_plot(mx_downstream_plot(comp), "app_cell_composition", height = 5)

# ------------------------------------------------------------------------------
# C. Consensus co-expression modules
# ------------------------------------------------------------------------------
message("\n== C. Consensus co-expression modules ==")

net <- mx_coexpression_modules(studies, n_genes = 3000, consensus = "min",
                               cut_height = 0.3, min_module_size = 20,
                               meta_result = meta_res)
net$modules$symbol <- unname(symbol_of[net$modules$gene_id])
ms <- net$module_stats
ms$hub_genes <- vapply(strsplit(ms$hub_genes, ", "), function(g) {
  sym <- symbol_of[g]
  paste(ifelse(is.na(sym), g, sym), collapse = ", ")
}, character(1))
net$module_stats <- ms
print(ms[, c("module", "n_genes", "meta_g", "padj", "i_squared",
             "pct_meta_up", "pct_meta_down", "hub_genes")],
      row.names = FALSE, digits = 3)
write.csv(net$modules, file.path(out_dir, "app_coexpression_modules.csv"),
          row.names = FALSE)
write.csv(ms, file.path(out_dir, "app_coexpression_module_stats.csv"),
          row.names = FALSE)
save_plot(mx_downstream_plot(net), "app_coexpression_modules", height = 5)

# ------------------------------------------------------------------------------
# D. Drug signature reversal
# ------------------------------------------------------------------------------
message("\n== D. Drug reversal (LINCS L1000 consensus signatures) ==")

sigs <- mx_read_signatures(file.path(in_dir, "lincs_l1000_chem_consensus.txt"))
drugs <- mx_drug_connectivity(meta_res, sigs, id_col = "symbol")
message(sprintf("  %d drug signatures scored", nrow(drugs)))
print(head(drugs, 15), row.names = FALSE, digits = 3)
write.csv(drugs, file.path(out_dir, "app_drug_connectivity.csv"),
          row.names = FALSE)
save_plot(mx_downstream_plot(drugs, top = 30), "app_drug_connectivity",
          height = 7)

# ------------------------------------------------------------------------------
# Combined four-panel figure
# ------------------------------------------------------------------------------
panel <- patchwork::wrap_plots(
  mx_downstream_plot(tf, top = 20, title = "A. Transcription factor activity"),
  mx_downstream_plot(comp, title = "B. Cell-type composition") +
    guides(alpha = "none"),
  mx_downstream_plot(net, title = "C. Consensus co-expression modules") +
    guides(alpha = "none"),
  mx_downstream_plot(drugs, top = 20, title = "D. Drug signature reversal") +
    guides(alpha = "none"),
  ncol = 2, guides = "collect"
) & theme(legend.position = "bottom")
save_plot(panel, "app_downstream_panel", width = 13, height = 11)

# ------------------------------------------------------------------------------
# Summary metrics for the report
# ------------------------------------------------------------------------------
saveRDS(list(
  n_genes_tested = nrow(mt), n_meta_degs = sum(sig),
  n_up = sum(sig & mt$meta_log2FC > 0), n_down = sum(sig & mt$meta_log2FC < 0),
  tf = tf, composition = comp, module_stats = ms, drugs = head(drugs, 50),
  session = sessionInfo()
), file.path(out_dir, "app_metrics.rds"))
message("\nDone.")
