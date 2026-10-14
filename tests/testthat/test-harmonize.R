test_that("mx_align_genes finds common genes across studies", {
  s1 <- make_test_study(n_genes = 100, accession = "A")
  s2 <- make_test_study(n_genes = 80, accession = "B")
  # Give s2 a 20-gene overlap with s1
  rownames(s2)[seq_len(60)] <- paste0("UNIQUE_", seq_len(60))

  aligned <- mx_align_genes(list(A = s1, B = s2))
  n_common <- nrow(aligned$A)
  expect_equal(n_common, nrow(aligned$B))
  expect_true(n_common == 20)
})

test_that("mx_align_genes errors when no genes are shared", {
  s1 <- make_test_study(n_genes = 50, accession = "A")
  s2 <- make_test_study(n_genes = 50, accession = "B")
  rownames(s2) <- paste0("DIFF_", seq_len(50))

  expect_error(
    mx_align_genes(list(A = s1, B = s2)),
    "No common genes"
  )
})

test_that("mx_align_genes errors on empty list", {
  expect_error(mx_align_genes(list()), "non-empty list")
})

test_that("mx_normalize returns a study with same dimensions", {
  s <- make_test_study(n_genes = 50, n_samples = 6)
  s_norm <- mx_normalize(s, method = "CPM")
  expect_s4_class(s_norm, "metaXpressStudy")
  expect_equal(dim(s_norm), dim(s))
})

test_that("mx_normalize CPM produces non-negative values", {
  s      <- make_test_study(n_genes = 50, n_samples = 6)
  s_norm <- mx_normalize(s, method = "CPM")
  expect_true(all(SummarizedExperiment::assay(s_norm, "cpm") >= 0))
  # raw counts are never overwritten
  expect_identical(SummarizedExperiment::assay(s_norm, "counts"),
                   SummarizedExperiment::assay(s, "counts"))
})

test_that("mx_normalize TMM returns a numeric matrix", {
  s      <- make_test_study(n_genes = 50, n_samples = 6)
  s_norm <- mx_normalize(s, method = "TMM")
  expect_true(is.matrix(SummarizedExperiment::assay(s_norm, "logCPM_TMM")))
  expect_true(is.numeric(SummarizedExperiment::assay(s_norm, "logCPM_TMM")))
})

test_that("mx_remove_batch keeps raw counts and preserves condition", {
  skip_if_not_installed("limma")
  studies <- list(
    A = make_test_study(n_genes = 60, n_samples = 6, accession = "A", seed = 1),
    B = make_test_study(n_genes = 60, n_samples = 6, accession = "B", seed = 2))
  out <- mx_remove_batch(studies, method = "limma")
  expect_true("log2cpm_batch_corrected" %in%
                SummarizedExperiment::assayNames(out$A))
  expect_identical(SummarizedExperiment::assay(out$A, "counts"),
                   SummarizedExperiment::assay(studies$A, "counts"))
  expect_error(mx_remove_batch(studies, method = "harmony"), "removed")
})
