test_that("mx_de returns study with de_result slot populated", {
  s  <- make_test_study(n_genes = 100, n_samples = 6)
  s2 <- mx_de(s, method = "DESeq2")
  expect_s4_class(s2, "metaXpressStudy")
  expect_true(nrow(s2@de_result) > 0)
})

test_that("mx_de result has required columns", {
  s   <- make_test_study(n_genes = 100, n_samples = 6)
  s2  <- mx_de(s, method = "DESeq2")
  required <- c("gene_id", "log2FC", "pvalue", "padj", "baseMean")
  expect_true(all(required %in% colnames(s2@de_result)))
})

test_that("mx_de padj values are in [0, 1]", {
  s  <- make_test_study(n_genes = 100, n_samples = 6)
  s2 <- mx_de(s, method = "DESeq2")
  padj_vals <- s2@de_result$padj[!is.na(s2@de_result$padj)]
  expect_true(all(padj_vals >= 0 & padj_vals <= 1))
})

test_that("mx_de errors without condition column", {
  s <- make_test_study()
  s$condition <- NULL
  expect_error(mx_de(s, method = "DESeq2"), "condition")
})

test_that("mx_de_all runs on all studies and returns same length list", {
  studies <- list(
    A = make_test_study(n_genes = 80, n_samples = 6, accession = "A", seed = 1),
    B = make_test_study(n_genes = 80, n_samples = 6, accession = "B", seed = 2)
  )
  result <- mx_de_all(studies, method = "DESeq2")
  expect_equal(length(result), 2)
  expect_true(all(vapply(result, function(s) nrow(s@de_result) > 0,
                          logical(1))))
})

test_that("mx_de_summary returns a data.frame with correct columns", {
  studies <- list(
    A = mx_de(make_test_study(n_genes = 100, n_samples = 6, accession = "A",
                               seed = 1), method = "DESeq2"),
    B = mx_de(make_test_study(n_genes = 100, n_samples = 6, accession = "B",
                               seed = 2), method = "DESeq2")
  )
  summ <- mx_de_summary(studies)
  expect_true(is.data.frame(summ))
  expect_true(all(c("study", "n_total", "n_up", "n_down") %in%
                    colnames(summ)))
  expect_equal(nrow(summ), 2)
})

test_that("mx_de honours explicit contrast and records group sizes", {
  s  <- make_test_study(n_genes = 100, n_samples = 6)
  fwd <- mx_de(s, method = "DESeq2", contrast = c("case", "control"))
  rev <- mx_de(s, method = "DESeq2", contrast = c("control", "case"))
  expect_equal(fwd@de_result$log2FC, -rev@de_result$log2FC, tolerance = 1e-3)
  expect_equal(attr(fwd@de_result, "contrast"),
               c(case = "case", control = "control"))
  expect_equal(unique(fwd@de_result$n_case), 3L)
  expect_error(mx_de(s, contrast = c("case", "nope")), "not found")
})

test_that("all engines return log2-scale SE consistent with stat", {
  s <- make_test_study(n_genes = 150, n_samples = 8)
  for (m in c("DESeq2", "edgeR", "limma-voom")) {
    de <- mx_de(s, method = m)@de_result
    expect_true(all(c("lfcSE", "stat", "n_case", "n_control") %in% names(de)),
                info = m)
    ok <- !is.na(de$lfcSE) & !is.na(de$stat)
    if (m != "limma-voom")
      expect_equal(de$log2FC[ok] / de$lfcSE[ok], de$stat[ok],
                   tolerance = 1e-4, info = m)
  }
})

test_that("mx_de adjusts for covariates in formula", {
  s <- make_test_study(n_genes = 100, n_samples = 8)
  s$batch <- rep(c("b1", "b2"), 4)
  for (m in c("DESeq2", "edgeR", "limma-voom"))
    expect_s4_class(mx_de(s, method = m, formula = ~ batch + condition),
                    "metaXpressStudy")
  expect_error(mx_de(s, formula = ~ nope + condition), "not found")
  expect_error(mx_de(s, formula = ~ batch), "condition")
})

test_that("reference detection does not match substrings", {
  expect_equal(.identify_reference_level(c("abnormal", "control")), "control")
  expect_equal(.identify_reference_level(c("AD", "Normal")), "Normal")
})

test_that("mx_de matches native DESeq2 on airway (GSE52778)", {
  skip_if_not_installed("airway")
  data("airway", package = "airway", envir = environment())
  st <- mx_from_se(airway, condition_col = "dex", accession = "GSE52778")
  st <- st[rowSums(SummarizedExperiment::assay(st)) >= 10, ]
  de <- suppressMessages(mx_de(st, formula = ~ cell + condition,
                               contrast = c("trt", "untrt")))@de_result
  dds <- DESeq2::DESeqDataSet(airway[rownames(st), ], ~ cell + dex)
  dds$dex <- stats::relevel(dds$dex, "untrt")
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  ref <- DESeq2::results(dds, name = "dex_trt_vs_untrt")
  expect_equal(de$log2FC, unname(ref$log2FoldChange))
  expect_equal(de$padj, unname(ref$padj))
})
