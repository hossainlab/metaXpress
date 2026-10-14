test_that("make_test_study produces a valid metaXpressStudy", {
  s <- make_test_study()
  expect_s4_class(s, "metaXpressStudy")
  expect_equal(s@accession, "TEST001")
  expect_true(nrow(s) == 200)
  expect_true(ncol(s) == 6)
})

test_that("mx_qc_study fills qc_score slot with criteria passed", {
  s    <- make_test_study()
  s_qc <- suppressWarnings(mx_qc_study(s))
  expect_s4_class(s_qc, "metaXpressStudy")
  qc <- S4Vectors::metadata(s_qc)$qc
  expect_equal(s_qc@qc_score, sum(qc$status == "pass"))
  expect_true(s_qc@qc_score >= 0 && s_qc@qc_score <= 10)
})

test_that("mx_qc_study reports pass/fail/not evaluable without rescaling", {
  s    <- make_test_study()
  s_qc <- suppressWarnings(mx_qc_study(s))
  qc   <- S4Vectors::metadata(s_qc)$qc
  expect_equal(nrow(qc), 10)
  expect_true(all(c("criterion", "name", "status", "score", "details") %in%
                    colnames(qc)))
  # alignment/rRNA/duplicates and batch are absent in the fixture
  expect_equal(qc$status[c(3, 4, 5, 9)], rep("not evaluable", 4))
  expect_identical(attr(s_qc, "qc_details"), qc)
  # fixture has 200 genes: gene detection fails, no rescaling hides it
  expect_equal(qc$status[6], "fail")
  expect_equal(suppressWarnings(mx_qc_study(s, min_genes = 100))@qc_score,
               s_qc@qc_score + 1)
})

test_that("mx_qc_study criterion 1 and 7 behave", {
  s2 <- make_test_study(n_samples = 4)
  s2$condition <- c("control", "control", "case", "case")
  qc <- S4Vectors::metadata(suppressWarnings(mx_qc_study(s2)))$qc
  expect_equal(qc$status[1], "fail")

  s3 <- make_test_study(n_samples = 6)
  s3$condition[1] <- NA
  qc <- S4Vectors::metadata(suppressWarnings(mx_qc_study(s3)))$qc
  expect_equal(qc$status[7], "fail")
})

test_that("mx_qc_study criterion 8 needs exactly 2 levels", {
  s3 <- make_test_study(n_samples = 6)
  s3$condition <- c("A", "B", "C", "A", "B", "C")
  qc <- S4Vectors::metadata(suppressWarnings(mx_qc_study(s3)))$qc
  expect_equal(qc$status[8], "fail")
})

test_that("mx_qc_study detects batch aliased with condition", {
  s <- make_test_study(n_samples = 6)
  s$batch <- s$condition
  qc <- S4Vectors::metadata(suppressWarnings(mx_qc_study(s)))$qc
  expect_equal(qc$status[9], "fail")
  s$batch <- rep(c("b1", "b2"), 3)
  qc <- S4Vectors::metadata(suppressWarnings(mx_qc_study(s)))$qc
  expect_equal(qc$status[9], "pass")
})

test_that("mx_filter_studies removes studies failing required criteria", {
  ok  <- suppressWarnings(mx_qc_study(make_test_study(accession = "A")))
  bad <- make_test_study(n_samples = 4, accession = "B")
  bad$condition <- c("control", "control", "case", "case")
  bad <- suppressWarnings(mx_qc_study(bad))
  expect_message(out <- mx_filter_studies(list(A = ok, B = bad)),
                 "required criterion failed: 1")
  expect_equal(names(out), "A")
  # gene-detection failure is not required by default
  expect_length(mx_filter_studies(list(A = ok), max_failed = 1), 1)
  expect_length(suppressMessages(
    mx_filter_studies(list(A = ok), max_failed = 0)), 0)
})

test_that("mx_filter_studies warns when QC was not run", {
  s1 <- make_test_study(accession = "A")
  expect_warning(suppressMessages(mx_filter_studies(list(A = s1))),
                 "no QC results")
})

test_that("mx_cluster_samples adds clusters without touching condition", {
  skip_if_not_installed("stringdist")
  s <- make_test_study(n_samples = 6)
  s$title <- c(paste("healthy donor", 1:3), paste("tumour patient", 1:3))
  out <- suppressMessages(mx_cluster_samples(s, columns = "title"))
  expect_identical(out$condition, s$condition)
  expect_equal(length(unique(out$cluster[1:3])), 1)
  expect_false(out$cluster[1] == out$cluster[4])
})

test_that("mx_load_local errors on mismatched path lengths", {
  expect_error(
    mx_load_local(count_paths = "a.csv", metadata_paths = c("b.csv", "c.csv")),
    "same length"
  )
})
