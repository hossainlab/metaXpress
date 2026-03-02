test_that("mx_missing_summary returns a binary matrix", {
  de  <- make_test_de_results(n_genes = 50, n_studies = 3)
  mat <- mx_missing_summary(de)
  expect_true(is.matrix(mat))
  expect_true(all(mat %in% c(0L, 1L)))
  expect_equal(ncol(mat), 3)
})

test_that("mx_missing_summary has coverage_pct attribute", {
  de  <- make_test_de_results(n_genes = 50, n_studies = 3)
  mat <- mx_missing_summary(de)
  pct <- attr(mat, "coverage_pct")
  expect_true(!is.null(pct))
  expect_true(all(pct >= 0 & pct <= 100))
})

test_that("mx_filter_coverage removes genes below threshold", {
  de <- make_test_de_results(n_genes = 50, n_studies = 3)
  # Remove GENE1 from studies 2 and 3 so it's only in 1 study
  de[[2]] <- de[[2]][de[[2]]$gene_id != "GENE1", ]
  de[[3]] <- de[[3]][de[[3]]$gene_id != "GENE1", ]

  filtered <- mx_filter_coverage(de, min_studies = 2)
  gene_sets <- lapply(filtered, function(d) d$gene_id)
  all_kept  <- Reduce(union, gene_sets)
  expect_false("GENE1" %in% all_kept)
})

test_that("mx_filter_coverage retains genes present in enough studies", {
  de <- make_test_de_results(n_genes = 50, n_studies = 3)
  filtered <- mx_filter_coverage(de, min_studies = 3)
  # All returned genes should be in all 3 studies
  cov_mat <- mx_missing_summary(de)
  kept    <- Reduce(union, lapply(filtered, function(d) d$gene_id))
  expect_true(all(rowSums(cov_mat[kept, , drop = FALSE]) >= 3))
})

test_that("mx_impute exclude method returns only common genes", {
  de <- make_test_de_results(n_genes = 50, n_studies = 3)
  de[[2]] <- de[[2]][de[[2]]$gene_id != "GENE1", ]

  result <- mx_impute(de, method = "exclude")
  all_genes <- lapply(result, function(d) d$gene_id)
  expect_false("GENE1" %in% all_genes[[1]])
  expect_true(length(unique(lengths(all_genes))) == 1)  # all same length
})

test_that("mx_impute mean method adds rows for missing genes", {
  de <- make_test_de_results(n_genes = 50, n_studies = 3)
  de[[2]] <- de[[2]][de[[2]]$gene_id != "GENE1", ]

  result <- mx_impute(de, method = "mean")
  # All studies should now have same gene count
  n_genes <- vapply(result, nrow, integer(1))
  expect_equal(length(unique(n_genes)), 1)
})
