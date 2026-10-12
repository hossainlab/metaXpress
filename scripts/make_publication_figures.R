# ==============================================================================
# metaXpress: Publication Figures Generation (Journal-Grade Minimalist Design)
# ==============================================================================
# Design Rules:
# 1. NO titles, subtitles, or captions inside the plot graphics (handled by LaTeX).
# 2. Vector PDF and 600 DPI PNG outputs for infinite resolution and zero pixelation.
# 3. Clean, breathable typography, thin line weights, generous whitespace.
# 4. Colorblind-safe, muted scientific palettes (Nature/Cell style).
# 5. ZERO label collisions, text overlaps, or margin clippings.
# ==============================================================================

suppressPackageStartupMessages({
  library(metaXpress)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
})

out_dirs <- c("case_study_results", "preprint/figures")
for (d in out_dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

save_both <- function(filename_base, plot_obj, width, height) {
  for (d in out_dirs) {
    png_path <- file.path(d, paste0(filename_base, ".png"))
    pdf_path <- file.path(d, paste0(filename_base, ".pdf"))
    
    # Save crisp vector PDF with cairo
    ggplot2::ggsave(pdf_path, plot_obj, width = width, height = height, device = cairo_pdf)
    # Save high-res raster PNG
    ggplot2::ggsave(png_path, plot_obj, width = width, height = height, dpi = 600)
    message("Saved: ", pdf_path, " (", width, "x", height, " in)")
  }
}

# ------------------------------------------------------------------------------
# Minimalist Academic Theme (No inside-plot titles/subtitles)
# ------------------------------------------------------------------------------
theme_journal <- function(base_size = 9, base_family = "sans") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 9, color = "#1E293B", face = "plain"),
      axis.text = ggplot2::element_text(size = 8, color = "#334155"),
      axis.line = ggplot2::element_line(color = "#475569", linewidth = 0.4),
      axis.ticks = ggplot2::element_line(color = "#475569", linewidth = 0.4),
      legend.title = ggplot2::element_text(size = 8.5, face = "bold", color = "#1E293B"),
      legend.text = ggplot2::element_text(size = 8, color = "#334155"),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.key.size = grid::unit(0.35, "cm"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 8.5, face = "bold", color = "#1E293B"),
      plot.margin = ggplot2::margin(6, 8, 6, 8)
    )
}

# ==============================================================================
# Figure 2: Study Overview & QC Benchmarking (2x2 Balanced Multi-Panel)
# ==============================================================================
message("Generating Figure 2: Study Overview & QC Benchmark...")

df_overview <- data.frame(
  cohort = factor(c("GSE130688\n(CRC)", "GSE136569\n(PDAC)", "GSE171485\n(ccRCC)"),
                  levels = c("GSE130688\n(CRC)", "GSE136569\n(PDAC)", "GSE171485\n(ccRCC)")),
  qc_score = c(10, 10, 10),
  n_normal = c(15, 5, 6),
  n_tumor = c(15, 5, 6),
  median_lib_m = c(42.1, 50.2, 36.8),
  aligned_genes_k = c(39.376, 39.376, 39.376)
)

# Panel A: QC Scores
p2a <- ggplot2::ggplot(df_overview, ggplot2::aes(x = cohort, y = qc_score)) +
  ggplot2::geom_col(fill = "#2563EB", width = 0.45) +
  ggplot2::geom_hline(yintercept = 7, linetype = "dashed", color = "#DC2626", linewidth = 0.45) +
  ggplot2::annotate("text", x = 3.3, y = 7.5, label = "Pass cutoff (7/10)",
                    color = "#DC2626", size = 2.8, hjust = 1, fontface = "italic") +
  ggplot2::scale_y_continuous(limits = c(0, 11), breaks = seq(0, 10, 2)) +
  ggplot2::labs(x = NULL, y = "QC score (0-10)") +
  theme_journal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, lineheight = 0.9))

# Panel B: Sample Distribution (Normal vs Tumor)
df_samp <- data.frame(
  cohort = rep(df_overview$cohort, 2),
  group = factor(rep(c("Normal", "Tumor"), each = 3), levels = c("Normal", "Tumor")),
  count = c(df_overview$n_normal, df_overview$n_tumor)
)

p2b <- ggplot2::ggplot(df_samp, ggplot2::aes(x = cohort, y = count, fill = group)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.6), width = 0.5) +
  ggplot2::scale_fill_manual(values = c("Normal" = "#64748B", "Tumor" = "#DC2626"), name = NULL) +
  ggplot2::scale_y_continuous(limits = c(0, 18), breaks = seq(0, 15, 5)) +
  ggplot2::labs(x = NULL, y = "Sample count (n)") +
  theme_journal() +
  ggplot2::theme(
    legend.position = c(0.82, 0.85),
    legend.background = ggplot2::element_rect(fill = "white", color = "#E2E8F0", linewidth = 0.3),
    axis.text.x = ggplot2::element_text(size = 8, lineheight = 0.9)
  )

# Panel C: Sequencing Depth
p2c <- ggplot2::ggplot(df_overview, ggplot2::aes(x = cohort, y = median_lib_m)) +
  ggplot2::geom_col(fill = "#D97706", width = 0.45) +
  ggplot2::geom_hline(yintercept = 20, linetype = "dashed", color = "#475569", linewidth = 0.45) +
  ggplot2::annotate("text", x = 3.3, y = 22.5, label = "ENCODE target (20M)",
                    color = "#475569", size = 2.8, hjust = 1, fontface = "italic") +
  ggplot2::scale_y_continuous(limits = c(0, 58), breaks = seq(0, 50, 10)) +
  ggplot2::labs(x = NULL, y = "Median depth (M reads)") +
  theme_journal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, lineheight = 0.9))

# Panel D: Quantified Genes
p2d <- ggplot2::ggplot(df_overview, ggplot2::aes(x = cohort, y = aligned_genes_k)) +
  ggplot2::geom_col(fill = "#059669", width = 0.45) +
  ggplot2::scale_y_continuous(limits = c(0, 45), breaks = seq(0, 40, 10)) +
  ggplot2::labs(x = NULL, y = "Quantified genes (\u00d710\u00b3)") +
  theme_journal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, lineheight = 0.9))

# Assemble 2x2 grid with ample room for labels
fig2 <- (p2a + p2b) / (p2c + p2d) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 11, face = "bold"))

save_both("case_study_study_overview", fig2, width = 6.4, height = 3.9)

# ==============================================================================
# Figure 3: Volcano Plot (Nature/Cell Style, Explicit Non-Crossing Leader Lines)
# ==============================================================================
message("Generating Figure 3: Publication Volcano Plot...")

meta_df <- read.csv("case_study_results/pan_cancer_meta_results.csv", stringsAsFactors = FALSE)
meta_df$symbol[is.na(meta_df$symbol) | meta_df$symbol == ""] <- meta_df$gene_id[is.na(meta_df$symbol) | meta_df$symbol == ""]

meta_df$neg_log10_padj <- -log10(meta_df$meta_padj)
meta_df$neg_log10_padj[meta_df$neg_log10_padj > 9.5] <- 9.5

padj_cut <- 0.05
lfc_cut  <- 1.0

meta_df$status <- "Non-significant"
meta_df$status[!is.na(meta_df$meta_padj) & meta_df$meta_padj <= padj_cut & meta_df$meta_log2FC >= lfc_cut] <- "Up-regulated (n = 377)"
meta_df$status[!is.na(meta_df$meta_padj) & meta_df$meta_padj <= padj_cut & meta_df$meta_log2FC <= -lfc_cut] <- "Down-regulated (n = 181)"

meta_df$status <- factor(meta_df$status, levels = c("Up-regulated (n = 377)", "Down-regulated (n = 181)", "Non-significant"))

# Define explicit label targets to ensure 100% clean non-crossing leader lines
# Up-regulated genes
df_up <- data.frame(
  symbol = c("MUC4", "IGF2BP3", "MST1R", "XDH", "CDCP1"),
  x_target = c(5.2, 5.2, 5.2, 5.2, 5.2),
  y_target = c(7.2, 6.4, 5.6, 4.8, 4.0),
  stringsAsFactors = FALSE
)
df_up <- merge(df_up, meta_df[, c("symbol", "meta_log2FC", "neg_log10_padj")], by = "symbol")

# Down-regulated genes
df_down <- data.frame(
  symbol = c("FAM107A", "BDKRB1", "SYPL2", "TTPA"),
  x_target = c(-5.2, -5.2, -5.2, -5.2),
  y_target = c(5.8, 5.0, 4.2, 3.4),
  stringsAsFactors = FALSE
)
df_down <- merge(df_down, meta_df[, c("symbol", "meta_log2FC", "neg_log10_padj")], by = "symbol")

fig3 <- ggplot2::ggplot(meta_df, ggplot2::aes(x = meta_log2FC, y = neg_log10_padj)) +
  # Background non-significant points
  ggplot2::geom_point(data = meta_df[meta_df$status == "Non-significant", ],
                      color = "#CBD5E1", alpha = 0.35, size = 0.7) +
  # Foreground significant points
  ggplot2::geom_point(data = meta_df[meta_df$status == "Up-regulated (n = 377)", ],
                      color = "#DC2626", alpha = 0.75, size = 1.3) +
  ggplot2::geom_point(data = meta_df[meta_df$status == "Down-regulated (n = 181)", ],
                      color = "#2563EB", alpha = 0.75, size = 1.3) +
  # Threshold lines (subtle)
  ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dotted",
                      color = "#64748B", linewidth = 0.4) +
  ggplot2::geom_hline(yintercept = -log10(padj_cut), linetype = "dotted",
                      color = "#64748B", linewidth = 0.4) +
  # Non-crossing leader segments for Up genes
  ggplot2::geom_segment(
    data = df_up,
    ggplot2::aes(x = meta_log2FC, y = neg_log10_padj, xend = x_target - 0.2, yend = y_target),
    color = "#475569", linewidth = 0.35
  ) +
  ggplot2::geom_text(
    data = df_up,
    ggplot2::aes(x = x_target, y = y_target, label = symbol),
    color = "#991B1B", fontface = "italic", size = 3.2, hjust = 0
  ) +
  # Non-crossing leader segments for Down genes
  ggplot2::geom_segment(
    data = df_down,
    ggplot2::aes(x = meta_log2FC, y = neg_log10_padj, xend = x_target + 0.2, yend = y_target),
    color = "#475569", linewidth = 0.35
  ) +
  ggplot2::geom_text(
    data = df_down,
    ggplot2::aes(x = x_target, y = y_target, label = symbol),
    color = "#1E40AF", fontface = "italic", size = 3.2, hjust = 1
  ) +
  ggplot2::scale_x_continuous(limits = c(-8.0, 8.0), breaks = seq(-8, 8, 2)) +
  ggplot2::scale_y_continuous(limits = c(0, 10.0), breaks = seq(0, 10, 2)) +
  ggplot2::labs(
    x = expression("Combined effect size ("*log[2]*" fold change)"),
    y = expression("-log"[10]*" (Random effects adjusted "*italic(P)*"-value)")
  ) +
  theme_journal() +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 9.5),
    axis.text = ggplot2::element_text(size = 8.5)
  )

save_both("case_study_volcano_plot", fig3, width = 5.8, height = 4.4)

# ==============================================================================
# Figure 4: Multi-Cohort Forest Plots (Lancet/BMJ Style, NO INSIDE TITLES)
# ==============================================================================
message("Generating Figure 4: Publication Forest Plots...")

build_journal_forest <- function(df_cohorts, df_pooled, accent_col = "#DC2626") {
  
  df_cohorts$study_clean <- factor(df_cohorts$study, levels = rev(df_cohorts$study))
  
  # Bounds
  x_plot_min <- -4.0
  x_plot_max <-  4.0
  
  # Distinct table columns with generous spacing to eliminate any overlap
  col_effect_x <- 4.6   # Starts effect size column
  col_ci_x     <- 6.3   # Starts 95% CI column
  col_stat_x   <- 10.2  # Starts Weight / P-value column
  total_x_max  <- 13.0
  
  p <- ggplot2::ggplot() +
    # Null reference line
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "#94A3B8", linewidth = 0.4) +
    # Cohort error bars
    ggplot2::geom_errorbar(
      data = df_cohorts,
      ggplot2::aes(y = study_clean, xmin = ci_lo, xmax = ci_hi),
      orientation = "y", width = 0.18, color = "#1E293B", linewidth = 0.6
    ) +
    # Cohort square points sized by weight
    ggplot2::geom_point(
      data = df_cohorts,
      ggplot2::aes(y = study_clean, x = log2FC, size = weight_pct),
      shape = 15, color = "#1E293B"
    ) +
    # Pooled summary diamond
    ggplot2::geom_polygon(
      data = df_pooled$diamond,
      ggplot2::aes(x = x, y = y),
      fill = accent_col, color = "#0F172A", linewidth = 0.5
    ) +
    # Separate columns for effect size, CI, and weight/P
    # 1. log2FC
    ggplot2::geom_text(
      data = df_cohorts,
      ggplot2::aes(y = study_clean, x = col_effect_x, label = sprintf("%+.2f", log2FC)),
      hjust = 0, size = 2.8, color = "#1E293B", family = "mono"
    ) +
    ggplot2::geom_text(
      data = df_pooled$summary_row,
      ggplot2::aes(y = label, x = col_effect_x, label = sprintf("%+.2f", log2FC)),
      hjust = 0, size = 2.8, fontface = "bold", color = accent_col, family = "mono"
    ) +
    # 2. 95% CI
    ggplot2::geom_text(
      data = df_cohorts,
      ggplot2::aes(y = study_clean, x = col_ci_x, label = sprintf("[%+.2f, %+.2f]", ci_lo, ci_hi)),
      hjust = 0, size = 2.8, color = "#1E293B", family = "mono"
    ) +
    ggplot2::geom_text(
      data = df_pooled$summary_row,
      ggplot2::aes(y = label, x = col_ci_x, label = sprintf("[%+.2f, %+.2f]", ci_lo, ci_hi)),
      hjust = 0, size = 2.8, fontface = "bold", color = accent_col, family = "mono"
    ) +
    # 3. Weight / P-value
    ggplot2::geom_text(
      data = df_cohorts,
      ggplot2::aes(y = study_clean, x = col_stat_x, label = sprintf("%5.1f%%", weight_pct)),
      hjust = 0, size = 2.8, color = "#475569", family = "mono"
    ) +
    ggplot2::geom_text(
      data = df_pooled$summary_row,
      ggplot2::aes(y = label, x = col_stat_x, label = sprintf("p = %.1e", padj)),
      hjust = 0, size = 2.8, fontface = "bold", color = accent_col, family = "mono"
    ) +
    # Column headers at top of panel
    ggplot2::annotate("text", x = col_effect_x, y = 4.8,
                      label = expression(bold(log[2]*"FC")),
                      hjust = 0, size = 2.7, color = "#334155") +
    ggplot2::annotate("text", x = col_ci_x, y = 4.8,
                      label = expression(bold("95% CI")),
                      hjust = 0, size = 2.7, color = "#334155") +
    ggplot2::annotate("text", x = col_stat_x, y = 4.8,
                      label = expression(bold("Weight / P")),
                      hjust = 0, size = 2.7, color = "#334155") +
    ggplot2::scale_size_continuous(range = c(2.5, 4.8), guide = "none") +
    ggplot2::scale_x_continuous(limits = c(x_plot_min, total_x_max), breaks = seq(-4, 4, 2)) +
    ggplot2::scale_y_discrete(limits = c("GSE171485 (ccRCC)", "GSE136569 (PDAC)", "GSE130688 (CRC)", "Pooled (RE)")) +
    ggplot2::coord_cartesian(ylim = c(0.6, 5.2), clip = "off") +
    ggplot2::labs(x = expression("Effect size ("*log[2]*" fold change, 95% CI)"), y = NULL) +
    theme_journal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8.5, face = "plain"),
      plot.margin = ggplot2::margin(10, 15, 6, 8)
    )
  
  return(p)
}

# MST1R Data
mst1r_cohorts <- data.frame(
  study = c("GSE130688 (CRC)", "GSE136569 (PDAC)", "GSE171485 (ccRCC)"),
  log2FC = c(2.296, 2.190, 1.302),
  se = c(0.456, 0.546, 0.975),
  weight_pct = c(52.2, 36.4, 11.4),
  stringsAsFactors = FALSE
)
mst1r_cohorts$ci_lo <- mst1r_cohorts$log2FC - 1.96 * mst1r_cohorts$se
mst1r_cohorts$ci_hi <- mst1r_cohorts$log2FC + 1.96 * mst1r_cohorts$se

mst1r_pooled_row <- data.frame(
  label = "Pooled (RE)",
  log2FC = 2.144,
  ci_lo = 2.144 - 1.96 * 0.329,
  ci_hi = 2.144 + 1.96 * 0.329,
  padj = 5.50e-7,
  weight_pct = 100
)

mst1r_diamond <- data.frame(
  x = c(mst1r_pooled_row$ci_lo, mst1r_pooled_row$log2FC, mst1r_pooled_row$ci_hi, mst1r_pooled_row$log2FC),
  y = c(4, 4 + 0.18, 4, 4 - 0.18)
)

mst1r_pooled <- list(summary_row = mst1r_pooled_row, diamond = mst1r_diamond)

p_forest_mst1r <- build_journal_forest(
  df_cohorts = mst1r_cohorts, df_pooled = mst1r_pooled, accent_col = "#DC2626"
)

# FAM107A Data
fam_cohorts <- data.frame(
  study = c("GSE130688 (CRC)", "GSE136569 (PDAC)", "GSE171485 (ccRCC)"),
  log2FC = c(-2.066, -2.337, -1.280),
  se = c(0.563, 0.588, 0.647),
  weight_pct = c(37.4, 34.3, 28.3),
  stringsAsFactors = FALSE
)
fam_cohorts$ci_lo <- fam_cohorts$log2FC - 1.96 * fam_cohorts$se
fam_cohorts$ci_hi <- fam_cohorts$log2FC + 1.96 * fam_cohorts$se

fam_pooled_row <- data.frame(
  label = "Pooled (RE)",
  log2FC = -1.937,
  ci_lo = -1.937 - 1.96 * 0.344,
  ci_hi = -1.937 + 1.96 * 0.344,
  padj = 2.39e-5,
  weight_pct = 100
)

fam_diamond <- data.frame(
  x = c(fam_pooled_row$ci_lo, fam_pooled_row$log2FC, fam_pooled_row$ci_hi, fam_pooled_row$log2FC),
  y = c(4, 4 + 0.18, 4, 4 - 0.18)
)

fam_pooled <- list(summary_row = fam_pooled_row, diamond = fam_diamond)

p_forest_fam <- build_journal_forest(
  df_cohorts = fam_cohorts, df_pooled = fam_pooled, accent_col = "#2563EB"
)

fig4 <- (p_forest_mst1r / p_forest_fam) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 11, face = "bold"))

save_both("case_study_forest_combined", fig4, width = 7.0, height = 4.6)

# ==============================================================================
# Figure 5: Between-Study Heterogeneity (I^2) Profile (Stacked Horizontal)
# ==============================================================================
message("Generating Figure 5: Heterogeneity Distribution...")

meta_valid_i2 <- meta_df[!is.na(meta_df$i_squared), ]
meta_sig_i2   <- meta_valid_i2[meta_valid_i2$status != "Non-significant", ]

# Panel A: Transcriptome-wide I2 histogram
p5a <- ggplot2::ggplot(meta_valid_i2, ggplot2::aes(x = i_squared)) +
  ggplot2::geom_histogram(binwidth = 5, fill = "#3B82F6", color = "white", linewidth = 0.2, boundary = 0) +
  ggplot2::geom_vline(xintercept = c(25, 50, 75), linetype = "dashed",
                      color = c("#0284C7", "#D97706", "#DC2626"), linewidth = 0.45) +
  ggplot2::scale_x_continuous(breaks = seq(0, 100, 25), limits = c(-5, 105)) +
  ggplot2::scale_y_continuous(labels = scales::comma, limits = c(0, 19500), expand = c(0, 0)) +
  ggplot2::labs(
    x = expression("Heterogeneity index ("*italic(I)^2*", %)"),
    y = "Number of genes"
  ) +
  theme_journal()

# Panel B: Compositional breakdown comparison (Horizontal Stacked Bar - Zero Collision)
het_stack_df <- data.frame(
  category = factor(c(
    rep("All evaluated genes\n(n = 29,560)", 4),
    rep("Pan-Cancer Meta-DEGs\n(n = 558)", 4)
  ), levels = c("Pan-Cancer Meta-DEGs\n(n = 558)", "All evaluated genes\n(n = 29,560)")),
  bracket = factor(rep(c("Low (\u226425%)", "Moderate (25-50%)", "Substantial (50-75%)", "High (>75%)"), 2),
                   levels = c("High (>75%)", "Substantial (50-75%)", "Moderate (25-50%)", "Low (\u226425%)")),
  pct = c(
    mean(meta_valid_i2$i_squared <= 25) * 100,
    mean(meta_valid_i2$i_squared > 25 & meta_valid_i2$i_squared <= 50) * 100,
    mean(meta_valid_i2$i_squared > 50 & meta_valid_i2$i_squared <= 75) * 100,
    mean(meta_valid_i2$i_squared > 75) * 100,
    mean(meta_sig_i2$i_squared <= 25) * 100,
    mean(meta_sig_i2$i_squared > 25 & meta_sig_i2$i_squared <= 50) * 100,
    mean(meta_sig_i2$i_squared > 50 & meta_sig_i2$i_squared <= 75) * 100,
    mean(meta_sig_i2$i_squared > 75) * 100
  )
)

# Reverse bracket factor order for proper left-to-right stacking
het_stack_df$bracket <- factor(het_stack_df$bracket,
                               levels = rev(c("Low (\u226425%)", "Moderate (25-50%)", "Substantial (50-75%)", "High (>75%)")))

p5b <- ggplot2::ggplot(het_stack_df, ggplot2::aes(y = category, x = pct, fill = bracket)) +
  ggplot2::geom_col(position = "fill", width = 0.55, color = "white", linewidth = 0.4) +
  ggplot2::scale_fill_manual(
    values = c("Low (\u226425%)" = "#2563EB",
               "Moderate (25-50%)" = "#F59E0B",
               "Substantial (50-75%)" = "#EA580C",
               "High (>75%)" = "#DC2626"),
    breaks = c("Low (\u226425%)", "Moderate (25-50%)", "Substantial (50-75%)", "High (>75%)"),
    labels = c("Low", "Moderate", "Substantial", "High"),
    name = "Heterogeneity"
  ) +
  ggplot2::scale_x_continuous(labels = scales::percent_format(scale = 100), expand = c(0, 0)) +
  # Print dominant labels inside bars cleanly
  ggplot2::annotate("text", y = 1, x = 0.48, label = "95.5% Low", color = "white",
                    fontface = "bold", size = 3.2) +
  ggplot2::annotate("text", y = 2, x = 0.33, label = "66.3% Low", color = "white",
                    fontface = "bold", size = 3.2) +
  ggplot2::labs(x = "Proportion of genes (%)", y = NULL) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE)) +
  theme_journal() +
  ggplot2::theme(
    legend.position = "top",
    legend.justification = "right",
    legend.title = ggplot2::element_text(size = 8, face = "bold"),
    legend.text = ggplot2::element_text(size = 7.5),
    legend.key.size = grid::unit(0.25, "cm"),
    legend.margin = ggplot2::margin(0, 0, 2, 0),
    axis.text.y = ggplot2::element_text(size = 8.5, lineheight = 0.95)
  )

fig5 <- (p5a | p5b) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 11, face = "bold"))

save_both("case_study_heterogeneity_plot", fig5, width = 7.2, height = 3.2)

message("\nAll journal-grade publication figures successfully updated in vector PDF and PNG!")
