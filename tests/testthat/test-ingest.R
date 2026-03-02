test_that("make_test_study produces a valid metaXpressStudy", {
  s <- make_test_study()
  expect_s4_class(s, "metaXpressStudy")
  expect_equal(s@accession, "TEST001")
  expect_true(nrow(s@counts) == 200)
  expect_true(ncol(s@counts) == 6)
})

test_that("mx_qc_study fills qc_score slot", {
  s   <- make_test_study()
  s_qc <- mx_qc_study(s)
  expect_s4_class(s_qc, "metaXpressStudy")
  expect_false(is.na(s_qc@qc_score))
  expect_true(s_qc@qc_score >= 0 && s_qc@qc_score <= 10)
})

test_that("mx_qc_study attaches qc_details attribute", {
  s    <- make_test_study()
  s_qc <- mx_qc_study(s)
  details <- attr(s_qc, "qc_details")
  expect_true(is.data.frame(details))
  expect_equal(nrow(details), 10)
  expect_true(all(c("criterion", "score", "details") %in% colnames(details)))
})

test_that("mx_qc_study scores criterion 1 correctly", {
  # 3 samples per group — should pass
  s   <- make_test_study(n_samples = 6)
  s_qc <- mx_qc_study(s)
  expect_equal(attr(s_qc, "qc_details")$score[1], 1L)

  # 2 samples per group — should fail
  s2 <- make_test_study(n_samples = 4)
  s2@metadata$condition <- c("control", "control", "case", "case")
  s2_qc <- mx_qc_study(s2)
  expect_equal(attr(s2_qc, "qc_details")$score[1], 0L)
})

test_that("mx_qc_study scores criterion 8 correctly (exactly 2 condition levels)", {
  s <- make_test_study()
  s_qc <- mx_qc_study(s)
  expect_equal(attr(s_qc, "qc_details")$score[8], 1L)

  # 3 levels — should fail
  s3 <- make_test_study(n_samples = 6)
  s3@metadata$condition <- c("A", "B", "C", "A", "B", "C")
  s3_qc <- mx_qc_study(s3)
  expect_equal(attr(s3_qc, "qc_details")$score[8], 0L)
})

test_that("mx_filter_studies removes low-QC studies", {
  s1 <- mx_qc_study(make_test_study(accession = "A"))
  s2 <- mx_qc_study(make_test_study(accession = "B"))
  s1@qc_score <- 8  # Force a passing score (fixture only scores 6 due to gene count)
  s2@qc_score <- 3  # Force a low score

  filtered <- mx_filter_studies(list(A = s1, B = s2), qc_threshold = 7)
  expect_equal(length(filtered), 1)
  expect_equal(names(filtered), "A")
})

test_that("mx_filter_studies warns on NA qc_score", {
  s1 <- make_test_study(accession = "A")
  s2 <- make_test_study(accession = "B")
  # Neither has been QC'd — both have NA qc_score
  expect_warning(
    mx_filter_studies(list(A = s1, B = s2)),
    "NA QC scores"
  )
})

test_that("mx_load_local errors on mismatched path lengths", {
  expect_error(
    mx_load_local(count_paths = "a.csv", metadata_paths = c("b.csv", "c.csv")),
    "same length"
  )
})
