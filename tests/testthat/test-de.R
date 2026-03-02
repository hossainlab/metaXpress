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
  s@metadata$condition <- NULL
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
