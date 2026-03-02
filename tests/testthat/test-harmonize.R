test_that("mx_align_genes finds common genes across studies", {
  s1 <- make_test_study(n_genes = 100, accession = "A")
  s2 <- make_test_study(n_genes = 80, accession = "B")
  # Give s2 a 20-gene overlap with s1
  rownames(s2@counts)[seq_len(60)] <- paste0("UNIQUE_", seq_len(60))

  aligned <- mx_align_genes(list(A = s1, B = s2))
  n_common <- nrow(aligned$A@counts)
  expect_equal(n_common, nrow(aligned$B@counts))
  expect_true(n_common == 20)
})

test_that("mx_align_genes errors when no genes are shared", {
  s1 <- make_test_study(n_genes = 50, accession = "A")
  s2 <- make_test_study(n_genes = 50, accession = "B")
  rownames(s2@counts) <- paste0("DIFF_", seq_len(50))

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
  expect_equal(dim(s_norm@counts), dim(s@counts))
})

test_that("mx_normalize CPM produces non-negative values", {
  s      <- make_test_study(n_genes = 50, n_samples = 6)
  s_norm <- mx_normalize(s, method = "CPM")
  expect_true(all(s_norm@counts >= 0))
})

test_that("mx_normalize TMM returns a numeric matrix", {
  s      <- make_test_study(n_genes = 50, n_samples = 6)
  s_norm <- mx_normalize(s, method = "TMM")
  expect_true(is.matrix(s_norm@counts))
  expect_true(is.numeric(s_norm@counts))
})
