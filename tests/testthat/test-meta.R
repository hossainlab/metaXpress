test_that("mx_meta returns a metaXpressResult object", {
  de <- make_test_de_results()
  result <- mx_meta(de, method = "random_effects")
  expect_s4_class(result, "metaXpressResult")
})

test_that("meta_table has all required columns", {
  de <- make_test_de_results()
  result <- mx_meta(de, method = "fisher")
  required <- c("gene_id", "meta_log2FC", "meta_pvalue", "meta_padj",
                 "i_squared", "n_studies", "direction_consistency")
  expect_true(all(required %in% colnames(result@meta_table)))
})

test_that("meta_padj values are in [0, 1]", {
  de <- make_test_de_results()
  for (method in c("fisher", "stouffer", "random_effects", "fixed_effects")) {
    result <- mx_meta(de, method = method)
    padj <- result@meta_table$meta_padj
    expect_true(all(padj >= 0 & padj <= 1, na.rm = TRUE),
                info = paste("method:", method))
  }
})

test_that("mx_meta n_studies slot matches input length", {
  de     <- make_test_de_results(n_studies = 3)
  result <- mx_meta(de, method = "random_effects")
  expect_equal(result@n_studies, 3L)
})

test_that("mx_meta warns with fewer than 3 studies for random_effects", {
  de <- make_test_de_results(n_studies = 2)
  expect_warning(mx_meta(de, method = "random_effects"),
                 "fewer than 3 studies")
})

test_that("mx_meta errors on invalid method", {
  de <- make_test_de_results()
  expect_error(mx_meta(de, method = "not_a_method"))
})

test_that("mx_meta errors on missing required columns", {
  de <- make_test_de_results()
  de[[1]]$padj <- NULL
  expect_error(mx_meta(de, method = "fisher"), "missing columns")
})

test_that("mx_meta min_studies filters genes correctly", {
  # 3 studies, but gene GENE1 only in study 1
  de <- make_test_de_results(n_genes = 50, n_studies = 3)
  de[[2]] <- de[[2]][de[[2]]$gene_id != "GENE1", ]
  de[[3]] <- de[[3]][de[[3]]$gene_id != "GENE1", ]

  result_strict <- mx_meta(de, method = "fisher", min_studies = 3)
  result_loose  <- mx_meta(de, method = "fisher", min_studies = 2)

  expect_false("GENE1" %in% result_strict@meta_table$gene_id)
  expect_true("GENE1"  %in% result_loose@meta_table$gene_id)
})

test_that("mx_meta direction_consistency is in [0, 1]", {
  de <- make_test_de_results()
  result <- mx_meta(de, method = "stouffer")
  dc <- result@meta_table$direction_consistency
  expect_true(all(dc >= 0 & dc <= 1, na.rm = TRUE))
})

test_that("mx_heterogeneity returns data.frame with correct columns", {
  de  <- make_test_de_results()
  het <- mx_heterogeneity(de)
  expect_true(is.data.frame(het))
  expect_true(all(c("gene_id", "Q", "I_sq", "tau_sq", "p_heterogeneity") %in%
                    colnames(het)))
})

test_that("mx_sensitivity returns a list of length n-1", {
  de  <- make_test_de_results(n_genes = 50, n_studies = 4)
  loo <- mx_sensitivity(de, method = "fisher")
  expect_equal(length(loo), 4)
  expect_true(all(vapply(loo, is, logical(1), "metaXpressResult")))
})

test_that("mx_sensitivity errors with fewer than 3 studies", {
  de <- make_test_de_results(n_studies = 2)
  expect_error(mx_sensitivity(de), "at least 3 studies")
})
