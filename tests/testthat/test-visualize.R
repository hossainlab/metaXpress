test_that("mx_volcano returns a ggplot object", {
  de     <- make_test_de_results(n_genes = 100, n_studies = 3)
  result <- mx_meta(de, method = "random_effects")
  p      <- mx_volcano(result)
  expect_s3_class(p, "ggplot")
})

test_that("mx_volcano errors on wrong input type", {
  expect_error(mx_volcano(list()), "metaXpressResult")
})

test_that("mx_forest returns a ggplot for a gene present in all studies", {
  de <- make_test_de_results(n_genes = 50, n_studies = 3)
  names(de) <- c("StudyA", "StudyB", "StudyC")
  p  <- mx_forest("GENE1", de)
  expect_s3_class(p, "ggplot")
})

test_that("mx_forest errors on gene not found in any study", {
  de <- make_test_de_results(n_genes = 50, n_studies = 3)
  expect_error(mx_forest("NOTEXIST", de), "not found")
})

test_that("mx_heterogeneity_plot returns a ggplot", {
  de     <- make_test_de_results(n_genes = 50, n_studies = 3)
  result <- mx_meta(de, method = "random_effects")
  p      <- mx_heterogeneity_plot(result)
  expect_s3_class(p, "ggplot")
})

test_that("mx_heterogeneity_plot errors on wrong input", {
  expect_error(mx_heterogeneity_plot("not_a_result"), "metaXpressResult")
})

test_that("mx_upset returns a plot object when genes are significant", {
  set.seed(1)
  de <- make_test_de_results(n_genes = 100, n_studies = 3)
  # Force some significant genes
  for (i in seq_along(de)) {
    de[[i]]$padj[1:10]   <- 0.001
    de[[i]]$log2FC[1:10] <- 2
  }
  names(de) <- c("StudyA", "StudyB", "StudyC")
  # Either ComplexHeatmap UpSet or ggplot fallback
  p <- mx_upset(de, padj_threshold = 0.05, lfc_threshold = 1)
  expect_true(!is.null(p))
})

test_that("mx_upset errors when no significant genes", {
  de <- make_test_de_results(n_genes = 50, n_studies = 2)
  for (i in seq_along(de)) de[[i]]$padj <- 1
  expect_error(mx_upset(de), "No significant genes")
})

test_that("mx_study_overview returns a named list of plots", {
  s1 <- make_test_study(n_genes = 100, accession = "GSE001")
  s2 <- make_test_study(n_genes = 100, accession = "GSE002", seed = 99)
  plots <- mx_study_overview(list(S1 = s1, S2 = s2))
  expect_type(plots, "list")
  expect_true("qc_scores" %in% names(plots))
  expect_true("sample_sizes" %in% names(plots))
  expect_s3_class(plots$sample_sizes, "ggplot")
})

test_that("mx_study_overview errors on empty list", {
  expect_error(mx_study_overview(list()), "non-empty")
})
