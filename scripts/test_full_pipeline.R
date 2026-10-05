# =============================================================================
# metaXpress Full Pipeline Test
# =============================================================================
# Tests the complete workflow using attached datasets.


library(metaXpress)


# ── 1. Load studies with raw integer counts (ENTREZID)
studies <- mx_load_local(
  count_paths = c("data/GSE130688_raw_counts.tsv", "data/GSE136569_raw_counts.tsv"),
  metadata_paths = c("data/GSE130688_metadata.csv", "data/GSE136569_metadata.csv"),
  organism = "Homo sapiens"
)



# ── 2. QC summary ────────────────────────────────────────────────────────────
mx_qc_study(studies$GSE136569_raw_counts)
attr(studies$GSE136569_raw_counts, "qc_details")



for (nm in names(studies)) {
  s <- studies[[nm]]
  message(sprintf("  %s: QC score = %d/10, samples = %d, genes = %d",
                  nm, s@qc_score, ncol(s@counts), nrow(s@counts)))
}

# Filter to high-quality studies (score >= 5)
studies <- mx_filter_studies(studies, qc_threshold = 5)

message(sprintf("  Retained %d studies after QC filtering.", length(studies)))

# ── 3. Gene alignment ────────────────────────────────────────────────────────
message("\n=== Step 3: Aligning genes ===")

studies <- mx_align_genes(studies)

message(sprintf("  Common genes: %d", nrow(studies[[1]]@counts)))

# ── 4. Per-study Differential Expression ──────────────────────────────────────
message("\n=== Step 4: Per-study DE (DESeq2) ===")

studies <- mx_de_all(studies, method = "DESeq2")

summary_tbl <- mx_de_summary(studies)
print(summary_tbl)

# ── 5. Missing gene handling ──────────────────────────────────────────────────
message("\n=== Step 5: Missing gene summary ===")

de_results <- lapply(studies, function(s) s@de_result)

cov_mat <- mx_missing_summary(de_results)
cov_pct <- attr(cov_mat, "coverage_pct")
message(sprintf("  Genes in all studies: %d / %d (%.1f%%)",
                sum(cov_pct == 100), length(cov_pct),
                sum(cov_pct == 100) / length(cov_pct) * 100))

# ── 6. Meta-analysis ──────────────────────────────────────────────────────────
message("\n=== Step 6: Meta-analysis (random effects) ===")

n_samples <- vapply(studies, function(s) ncol(s@counts), integer(1))

meta_result <- mx_meta(
  de_results,
  method      = "random_effects",
  min_studies = 2,
  n_samples   = n_samples
)

mt <- meta_result@meta_table
sig <- mt[!is.na(mt$meta_padj) & mt$meta_padj <= 0.05 &
            abs(mt$meta_log2FC) >= 1, ]

message(sprintf("  Total genes tested: %d", nrow(mt)))
message(sprintf("  Significant (padj<=0.05, |LFC|>=1): %d", nrow(sig)))
message(sprintf("    Up-regulated: %d", sum(sig$meta_log2FC > 0)))
message(sprintf("    Down-regulated: %d", sum(sig$meta_log2FC < 0)))

# Top 10 DEGs
top10 <- head(mt[order(mt$meta_padj), ], 10)
message("\n  Top 10 genes by meta_padj:")
print(top10[, c("gene_id", "meta_log2FC", "meta_padj", "i_squared",
                "n_studies", "direction_consistency")])

# ── 7. Visualization ──────────────────────────────────────────────────────────
message("\n=== Step 7: Visualization ===")

# Volcano plot
p_volcano <- mx_volcano(meta_result, label_top = 15)
message("  Volcano plot created.")

# Heterogeneity plot
p_het <- mx_heterogeneity_plot(meta_result)
message("  I-squared distribution plot created.")

# Forest plot for top gene
if (nrow(sig) > 0) {
  top_gene <- sig$gene_id[which.min(sig$meta_padj)]
  p_forest <- mx_forest(top_gene, de_results, meta_result)
  message(sprintf("  Forest plot created for %s.", top_gene))
}

# ── 8. Export ──────────────────────────────────────────────────────────────────
message("\n=== Step 8: Export ===")

mx_export(meta_result, format = "csv", output_dir = tempdir(),
          prefix = "metaXpress_test")
message(sprintf("  Results exported to %s", tempdir()))

message("\n=== Pipeline completed successfully! ===\n")
