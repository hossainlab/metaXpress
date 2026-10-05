test_that("mx_pathway_consensus returns correct data.frame structure", {
  # Build two mock pathway result data.frames
  pr1 <- data.frame(
    pathway_id   = c("PATH_A", "PATH_B", "PATH_C"),
    pathway_name = c("Path A", "Path B", "Path C"),
    padj         = c(0.01, 0.04, 0.2),
    stringsAsFactors = FALSE
  )
  pr2 <- data.frame(
    pathway_id   = c("PATH_A", "PATH_C", "PATH_D"),
    pathway_name = c("Path A", "Path C", "Path D"),
    padj         = c(0.02, 0.3, 0.01),
    stringsAsFactors = FALSE
  )

  res <- mx_pathway_consensus(list(pr1, pr2), min_fraction = 0.5,
                               padj_threshold = 0.05)
  expect_s3_class(res, "data.frame")
  expect_true("pathway_id" %in% colnames(res))
  expect_true("fraction_studies" %in% colnames(res))
  # PATH_A is significant in both studies => should appear
  expect_true("PATH_A" %in% res$pathway_id)
  # PATH_C is not significant in either (padj > 0.05) => should not appear
  expect_false("PATH_C" %in% res$pathway_id)
})

test_that("mx_pathway_consensus errors on non-list input", {
  expect_error(mx_pathway_consensus("not_a_list"), "list")
})

test_that("mx_pathway_consensus returns empty data.frame when no consensus", {
  pr1 <- data.frame(pathway_id = "PATH_A", padj = 0.6,
                    stringsAsFactors = FALSE)
  pr2 <- data.frame(pathway_id = "PATH_B", padj = 0.6,
                    stringsAsFactors = FALSE)
  res <- mx_pathway_consensus(list(pr1, pr2), min_fraction = 0.5,
                               padj_threshold = 0.05)
  expect_equal(nrow(res), 0)
})

test_that("mx_pathway_dedup returns data.frame", {
  pr <- data.frame(
    pathway_id   = c("PATH_A", "PATH_B"),
    pathway_name = c("Path A", "Path B"),
    padj         = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )
  # Without msigdbr available this warns and returns input unchanged
  res <- suppressWarnings(mx_pathway_dedup(pr))
  expect_s3_class(res, "data.frame")
})

test_that("mx_pathway_dedup errors on non-data.frame input", {
  expect_error(mx_pathway_dedup(list()), "data.frame")
})

test_that("mx_pathway_dedup handles empty input", {
  empty <- data.frame(pathway_id = character(0), stringsAsFactors = FALSE)
  res <- mx_pathway_dedup(empty)
  expect_equal(nrow(res), 0)
})

test_that("mx_pathway_heatmap returns a ggplot for a single data.frame", {
  pr <- data.frame(
    pathway_id   = paste0("PATH_", LETTERS[1:5]),
    pathway_name = paste0("Path ", LETTERS[1:5]),
    padj         = c(0.001, 0.01, 0.02, 0.04, 0.05),
    NES          = c(1.5, -1.2, 0.8, -0.5, 1.1),
    stringsAsFactors = FALSE
  )
  p <- mx_pathway_heatmap(pr, top_n = 5)
  expect_s3_class(p, "ggplot")
})

test_that("mx_pathway_heatmap returns a ggplot for a list of data.frames", {
  pr1 <- data.frame(
    pathway_id = paste0("PATH_", 1:3),
    padj       = c(0.01, 0.03, 0.04),
    stringsAsFactors = FALSE
  )
  pr2 <- data.frame(
    pathway_id = paste0("PATH_", 2:4),
    padj       = c(0.02, 0.01, 0.05),
    stringsAsFactors = FALSE
  )
  p <- mx_pathway_heatmap(list(S1 = pr1, S2 = pr2), top_n = 3)
  expect_s3_class(p, "ggplot")
})
