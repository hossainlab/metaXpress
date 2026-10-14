# ==============================================================================
# Compute every number and figure reported in preprint/metaXpress.tex
# ==============================================================================
# Inputs:
#   inst/extdata/ad_brain_de.rds            (inst/scripts/make_ad_brain_extdata.R)
#   geo_cache/case_study/pipeline_objects.rds and case_study_results/
#                                           (scripts/run_case_study_applications.R)
# Outputs:
#   preprint/figures/*.pdf
#   preprint/results.txt                    key = value, cited in the manuscript
# Run from the package root:  Rscript scripts/make_preprint_results.R
# ==============================================================================

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(ggplot2)
})
set.seed(2026)
fig_dir <- "preprint/figures"
dir.create(fig_dir, showWarnings = FALSE)
res_file <- "preprint/results.txt"
cat("", file = res_file)
put <- function(key, value) {
  cat(sprintf("%s = %s\n", key, paste(format(value), collapse = "; ")),
      file = res_file, append = TRUE)
}
save_fig <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = 300)
}
quiet <- function(x) suppressMessages(suppressWarnings(x))

# ------------------------------------------------------------------------------
# 1. Statistical validation against reference implementations
# ------------------------------------------------------------------------------
ad <- readRDS("inst/extdata/ad_brain_de.rds")
de <- ad$de_results
genes <- sample(Reduce(intersect, lapply(de, `[[`, "gene_id")), 300)
maxdiff <- c()
for (tm in c("DL", "REML", "PM")) for (hk in c(FALSE, TRUE)) {
  mt <- quiet(mx_meta(de, tau_method = tm, hksj = hk))@meta_table
  d <- vapply(genes, function(g) {
    y <- vapply(de, function(x) x$log2FC[match(g, x$gene_id)], numeric(1))
    s <- vapply(de, function(x) x$lfcSE[match(g, x$gene_id)], numeric(1))
    p <- vapply(de, function(x) x$pvalue[match(g, x$gene_id)], numeric(1))
    ok <- !is.na(p)   # DESeq2 Cook's outliers (p = NA) are not measured
    y <- y[ok]; s <- s[ok]
    r <- metafor::rma(yi = y, sei = s, method = tm,
                      test = if (hk) "knha" else "z",
                      control = list(threshold = 1e-12, tol = 1e-12,
                                     maxiter = 5000, stepadj = 0.5))
    i <- match(g, mt$gene_id)
    max(abs(c(mt$meta_log2FC[i] - r$b[1], mt$meta_se[i] - r$se,
              mt$tau2[i] - r$tau2, mt$meta_pvalue[i] - r$pval)))
  }, numeric(1))
  maxdiff[paste(tm, if (hk) "HKSJ" else "Wald")] <- max(d)
}
put("validation_genes", length(genes))
put("validation_max_abs_diff_vs_metafor", signif(max(maxdiff), 2))

data("airway", package = "airway")
st <- mx_from_se(airway, condition_col = "dex", accession = "GSE52778")
st <- st[rowSums(SummarizedExperiment::assay(st)) >= 10, ]
eng <- lapply(c("DESeq2", "edgeR", "limma-voom"), function(m)
  quiet(mx_de(st, method = m, formula = ~ cell + condition,
              contrast = c("trt", "untrt")))@de_result)
dds <- DESeq2::DESeqDataSet(airway[rownames(st), ], ~ cell + dex)
dds$dex <- relevel(dds$dex, "untrt")
nat <- DESeq2::results(DESeq2::DESeq(dds, quiet = TRUE),
                       name = "dex_trt_vs_untrt")
put("airway_deseq2_max_abs_diff_log2FC",
    signif(max(abs(eng[[1]]$log2FC - nat$log2FoldChange), na.rm = TRUE), 2))
g <- Reduce(intersect, lapply(eng, `[[`, "gene_id"))
lfc <- sapply(eng, function(x) x$log2FC[match(g, x$gene_id)])
put("airway_engine_min_pearson_r", round(min(cor(lfc)), 3))

# ------------------------------------------------------------------------------
# 2. Alzheimer's disease brain case study
# ------------------------------------------------------------------------------
put("ad_studies", paste(ad$studies$accession, ad$studies$region,
                        ad$studies$n_AD, ad$studies$n_control,
                        ad$studies$n_genes, ad$studies$n_DEG_fdr05,
                        sep = "|"))
cc <- quiet(mx_study_concordance(de))
put("ad_concordance", paste(cc$study_1, cc$study_2,
                            round(cc$correlation, 2), sep = "|"))

res <- quiet(mx_meta(de, method = "random_effects", tau_method = "REML"))
mt  <- res@meta_table
sig <- mt$meta_padj <= 0.05
put("ad_genes_tested", nrow(mt))
put("ad_genes_all4", sum(mt$n_studies == 4))
put("ad_sig", sum(sig))
put("ad_sig_up", sum(sig & mt$meta_log2FC > 0))
put("ad_sig_down", sum(sig & mt$meta_log2FC < 0))
put("ad_sig_pi_excl0", sum(sig & (mt$pi_lower > 0 | mt$pi_upper < 0),
                           na.rm = TRUE))
put("ad_sig_median_I2", round(median(mt$i_squared[sig]), 1))

methods <- c("random_effects", "fixed_effects", "inverse_normal", "fisher",
             "rop", "maxp", "awmeta")
n_m <- vapply(methods, function(m) {
  args <- list(de, method = m)
  if (m == "rop") args$r <- 2
  sum(quiet(do.call(mx_meta, args))@meta_table$meta_padj <= 0.05)
}, integer(1))
put("ad_method_n_sig", paste(names(n_m), n_m, sep = "|"))

hk <- quiet(mx_meta(de, tau_method = "REML", hksj = TRUE))@meta_table
put("ad_hksj_min_p", signif(min(hk$meta_pvalue), 3))
put("ad_hksj_n_sig", sum(hk$meta_padj <= 0.05))

loo <- quiet(mx_sensitivity(de, tau_method = "REML"))
put("ad_loo_n_sig", paste(names(loo), vapply(loo, function(r)
  sum(r@meta_table$meta_padj <= 0.05), integer(1)), sep = "|"))
rob  <- attr(loo, "robustness")
core <- rob$gene_id[rob$full_padj <= 0.05 & rob$frac_loo_sig == 1]
put("ad_robust_core", length(core))

ss <- quiet(mx_small_study(de, genes = core))
put("ad_egger_tested", nrow(ss))
put("ad_egger_fdr05", sum(ss$egger_padj <= 0.05, na.rm = TRUE))

aw <- quiet(mx_meta(de, method = "awmeta"))@meta_table
aw_sig <- aw[aw$meta_padj <= 0.05, ]
put("ad_aw_sig_only_GSE125583",
    sum(aw_sig$aw_studies == "GSE125583"))
put("ad_aw_sig_all4", sum(aw_sig$aw_n_selected == 4))

known <- c("GFAP", "C4B", "CD44", "VIP", "VGF", "BDNF", "CRH", "NPTX2",
           "SST", "PVALB")
k <- mt[match(known, mt$gene_id), ]
for (i in seq_len(nrow(k)))
  put(paste0("ad_gene_", k$gene_id[i]),
      sprintf("lfc=%.2f ci=[%.2f,%.2f] pi=[%.2f,%.2f] padj=%.2g I2=%.0f",
              k$meta_log2FC[i], k$ci_lower[i], k$ci_upper[i], k$pi_lower[i],
              k$pi_upper[i], k$meta_padj[i], k$i_squared[i]))

res <- quiet(mx_pathway_meta(res, db = "Hallmarks"))
pr  <- res@pathway_result
put("ad_gsea_top", paste(head(pr$pathway_id, 8), round(head(pr$NES, 8), 2),
                         signif(head(pr$padj, 8), 2), sep = "|"))
ps <- quiet(mx_pathway_studies(de, db = "Hallmarks"))
put("ad_pathway_studies_top",
    paste(head(ps$pathway_id, 6), head(ps$n_up, 6), head(ps$n_down, 6),
          signif(head(ps$meta_padj, 6), 2), sep = "|"))

# Figures
p_vol <- mx_volcano(res, padj_threshold = 0.05, lfc_threshold = 0.5,
                    label_top = 12, title = "")
save_fig(p_vol, "fig2_ad_volcano", 6.5, 5)
f1 <- mx_forest("GFAP", de, res)
f2 <- mx_forest("SST", de, res)
save_fig(patchwork::wrap_plots(f1, f2, ncol = 1), "fig3_ad_forest", 7, 6.5)
save_fig(mx_heterogeneity_plot(res), "figS1_ad_heterogeneity", 6.5, 4.5)
save_fig(mx_pathway_heatmap(attr(ps, "per_study"), top_n = 15,
                            value = "NES"), "fig4_ad_pathways", 7, 5.5)
loo_df <- data.frame(
  analysis = c("All 4 studies", paste("without", names(loo))),
  n_sig = c(sum(sig), vapply(loo, function(r)
    sum(r@meta_table$meta_padj <= 0.05), integer(1))))
loo_df$analysis <- factor(loo_df$analysis, levels = rev(loo_df$analysis))
p_loo <- ggplot(loo_df, aes(n_sig, analysis)) + geom_col(fill = "#377EB8") +
  geom_vline(xintercept = length(core), linetype = "dashed",
             colour = "#E41A1C") +
  labs(x = "Genes with meta FDR <= 0.05", y = NULL) + theme_bw()
save_fig(p_loo, "figS2_ad_leave_one_out", 6, 3)

# ------------------------------------------------------------------------------
# 3. PDAC downstream-biology application
# ------------------------------------------------------------------------------
obj <- readRDS("geo_cache/case_study/pipeline_objects.rds")
pdac_de <- lapply(obj$studies, function(s) s@de_result)
put("pdac_studies", paste(names(obj$studies),
                          vapply(obj$studies, function(s)
                            sum(s$condition == "tumor"), integer(1)),
                          vapply(obj$studies, function(s)
                            sum(s$condition == "normal"), integer(1)),
                          vapply(pdac_de, function(d)
                            sum(d$padj <= 0.05, na.rm = TRUE), integer(1)),
                          sep = "|"))
pm <- obj$meta_res@meta_table
psig <- pm$meta_padj <= 0.05 & abs(pm$meta_log2FC) >= 1
put("pdac_genes_tested", nrow(pm))
put("pdac_meta_deg_lfc1", sum(psig))
put("pdac_meta_deg_up", sum(psig & pm$meta_log2FC > 0))
put("pdac_meta_deg_down", sum(psig & pm$meta_log2FC < 0))
pcc <- quiet(mx_study_concordance(pdac_de))
put("pdac_concordance", paste(pcc$study_1, pcc$study_2,
                              round(pcc$correlation, 2), sep = "|"))
tf <- read.csv("case_study_results/app_tf_activity.csv")
put("pdac_tf_scored", nrow(tf))
put("pdac_tf_sig", sum(tf$padj <= 0.05))
comp <- read.csv("case_study_results/app_cell_composition.csv")
put("pdac_comp", paste(comp$cell_type, round(comp$meta_g, 2),
                       signif(comp$padj, 2), round(comp$i_squared), sep = "|"))
mods <- read.csv("case_study_results/app_coexpression_module_stats.csv")
put("pdac_modules", paste(mods$module, mods$n_genes, round(mods$meta_g, 2),
                          signif(mods$padj, 2), mods$hub_genes, sep = "|"))
drugs <- read.csv("case_study_results/app_drug_connectivity.csv")
put("pdac_drugs_n", nrow(drugs))
put("pdac_drugs_top10", paste(head(drugs$drug, 10), collapse = ", "))
file.copy("case_study_results/app_downstream_panel.pdf",
          file.path(fig_dir, "fig5_pdac_downstream.pdf"), overwrite = TRUE)
file.copy("scripts/make_architecture_diagram.pdf",
          file.path(fig_dir, "fig1_architecture.pdf"), overwrite = TRUE)

put("metaXpress_version", as.character(packageVersion("metaXpress")))
put("R_version", R.version.string)
message("Done: ", res_file)
