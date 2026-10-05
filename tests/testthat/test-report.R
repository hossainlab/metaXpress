test_that("mx_session_info returns a named list with key entries", {
  info <- mx_session_info()
  expect_type(info, "list")
  expect_true("timestamp" %in% names(info))
  expect_true("session_info" %in% names(info))
  expect_true("metaXpress_version" %in% names(info))
  expect_s3_class(info$timestamp, "POSIXct")
})

test_that("mx_export saves RDS file and returns path invisibly", {
  de     <- make_test_de_results(n_genes = 50, n_studies = 3)
  res    <- mx_meta(de, method = "fisher")
  tmpdir <- tempdir()
  out    <- mx_export(res, format = "rds", output_dir = tmpdir,
                      prefix = "test_export")
  expect_true(file.exists(out))
  loaded <- readRDS(out)
  expect_s4_class(loaded, "metaXpressResult")
  unlink(out)
})

test_that("mx_export saves CSV and returns path", {
  de     <- make_test_de_results(n_genes = 30, n_studies = 2)
  res    <- mx_meta(de, method = "stouffer")
  tmpdir <- tempdir()
  out    <- mx_export(res, format = "csv", output_dir = tmpdir,
                      prefix = "test_csv")
  expect_true(file.exists(out))
  csv_df <- read.csv(out)
  expect_true("gene_id" %in% colnames(csv_df))
  unlink(out)
})

test_that("mx_export errors on invalid format", {
  de  <- make_test_de_results(n_genes = 20, n_studies = 2)
  res <- mx_meta(de, method = "fisher")
  expect_error(mx_export(res, format = "tsv"), "arg")
})

test_that("mx_export errors on non-metaXpressResult input", {
  expect_error(mx_export(list(), format = "rds"), "metaXpressResult")
})

test_that("mx_report errors when meta_result is not metaXpressResult", {
  expect_error(mx_report("not_result", list(), list()), "metaXpressResult")
})
