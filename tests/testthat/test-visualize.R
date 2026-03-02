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
