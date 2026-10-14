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
  # Without gene sets this warns and returns input unchanged
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

make_pw_de <- function(seed = 1, shift = 0.8, k = 3) {
  set.seed(seed)
  genes <- paste0("G", 1:1000)
  de <- lapply(seq_len(k), function(s) {
    lfc <- rnorm(1000, c(rep(shift, 40), rep(0, 960)), 0.3)
    p <- 2 * pnorm(-abs(lfc / 0.3))
    data.frame(gene_id = genes, log2FC = lfc, lfcSE = 0.3, pvalue = p,
               padj = p.adjust(p, "BH"), n_case = 5, n_control = 5)
  })
  names(de) <- paste0("S", seq_len(k))
  list(de = de, sets = list(UP = genes[1:40], NULLSET = genes[501:560],
                            OVERLAP = genes[1:38]))
}

test_that("mx_pathway_meta GSEA detects an enriched set", {
  skip_if_not_installed("fgsea")
  x   <- make_pw_de()
  res <- suppressMessages(mx_meta(x$de))
  res <- mx_pathway_meta(res, gene_sets = x$sets, min_size = 10)
  pr  <- res@pathway_result
  expect_true(all(c("pathway_id", "NES", "pvalue", "padj") %in% names(pr)))
  expect_lt(pr$padj[pr$pathway_id == "UP"], 0.01)
  expect_gt(pr$NES[pr$pathway_id == "UP"], 0)
  expect_gt(pr$padj[pr$pathway_id == "NULLSET"], 0.05)
})

test_that("mx_pathway_meta ORA matches phyper and splits direction", {
  x   <- make_pw_de()
  res <- suppressMessages(mx_meta(x$de))
  pr  <- mx_pathway_meta(res, method = "ORA", gene_sets = x$sets,
                         min_size = 10)@pathway_result
  up  <- pr[pr$pathway_id == "UP" & pr$direction == "up", ]
  mt  <- res@meta_table
  N   <- sum(!is.na(mt$meta_padj))
  n   <- sum(mt$meta_padj <= 0.05 & mt$meta_log2FC > 0)
  expect_equal(up$pvalue, phyper(up$n_overlap - 1, 40, N - 40, n,
                                 lower.tail = FALSE))
  expect_lt(up$padj, 1e-10)
})

test_that("mx_pathway_studies combines per-study GSEA directionally", {
  skip_if_not_installed("fgsea")
  x  <- make_pw_de()
  ps <- suppressMessages(mx_pathway_studies(x$de, gene_sets = x$sets,
                                            min_size = 10))
  expect_true(all(c("meta_z", "meta_padj", "NES_S1") %in% names(ps)))
  expect_equal(ps$pathway_id[1], "UP")
  expect_equal(ps$n_up[ps$pathway_id == "UP"], 3L)
  expect_length(attr(ps, "per_study"), 3)

  # opposite directions in two studies cancel
  y <- x$de
  y$S2$log2FC <- -y$S2$log2FC
  y$S2$stat   <- NULL
  y <- y[c("S1", "S2")]
  ps2 <- suppressMessages(mx_pathway_studies(y, gene_sets = x$sets,
                                             min_size = 10))
  expect_gt(ps2$meta_pvalue[ps2$pathway_id == "UP"], 0.05)
})

test_that("mx_pathway_dedup removes overlapping sets", {
  x  <- make_pw_de()
  df <- data.frame(pathway_id = c("UP", "OVERLAP", "NULLSET"),
                   padj = c(0.001, 0.002, 0.5))
  out <- mx_pathway_dedup(df, gene_sets = x$sets)
  expect_equal(out$pathway_id, c("UP", "NULLSET"))
})
