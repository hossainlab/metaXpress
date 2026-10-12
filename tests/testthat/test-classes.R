test_that("metaXpressStudy extends SummarizedExperiment", {
  s <- make_test_study(n_genes = 20, n_samples = 6)
  expect_s4_class(s, "SummarizedExperiment")
  expect_equal(SummarizedExperiment::assayNames(s), "counts")
  expect_equal(dim(s), c(20L, 6L))
  expect_equal(s$condition, rep(c("control", "case"), each = 3))
})

test_that("subsetting keeps class and extra slots", {
  s <- make_test_study(n_genes = 20, n_samples = 6, accession = "GSE1")
  sub <- s[1:5, 2:4]
  expect_s4_class(sub, "metaXpressStudy")
  expect_equal(dim(sub), c(5L, 3L))
  expect_equal(sub@accession, "GSE1")
  expect_equal(sub$sample_id, c("S2", "S3", "S4"))
})

test_that("validity requires a condition column once accession is set", {
  counts <- matrix(1L, 2, 2, dimnames = list(c("G1", "G2"), c("S1", "S2")))
  expect_error(
    metaXpressStudy(counts, data.frame(sample_id = c("S1", "S2")),
                    accession = "GSE1"),
    "condition"
  )
})

test_that("validity rejects out-of-range qc_score", {
  counts <- matrix(1L, 2, 2, dimnames = list(c("G1", "G2"), c("S1", "S2")))
  meta   <- data.frame(condition = c("control", "case"))
  expect_error(metaXpressStudy(counts, meta, accession = "GSE1", qc_score = 11),
               "qc_score")
})

test_that("constructor accepts empty counts", {
  s <- metaXpressStudy(
    metadata  = data.frame(sample_id = character(0), condition = character(0)),
    accession = "GSE1"
  )
  expect_equal(dim(s), c(0L, 0L))
})

test_that("mx_from_se converts a SummarizedExperiment", {
  counts <- matrix(rpois(24, 100), nrow = 4,
                   dimnames = list(paste0("G", 1:4), paste0("S", 1:6)))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = counts),
    colData = data.frame(dex = rep(c("untrt", "trt"), each = 3)),
    rowData = data.frame(gene_length = 1000 * (1:4))
  )
  s <- mx_from_se(se, condition_col = "dex", accession = "SRP1")
  expect_s4_class(s, "metaXpressStudy")
  expect_equal(s$condition, rep(c("untrt", "trt"), each = 3))
  expect_equal(s@accession, "SRP1")
  expect_equal(SummarizedExperiment::rowData(s)$gene_length, 1000 * (1:4))

  expect_error(mx_from_se(se), "condition_col")
})

test_that("as() round-trips through SummarizedExperiment", {
  s  <- make_test_study(n_genes = 10, n_samples = 4, accession = "GSE9")
  s@qc_score <- 8
  se <- as(s, "SummarizedExperiment")
  expect_identical(class(se)[1], "SummarizedExperiment")
  expect_equal(S4Vectors::metadata(se)$accession, "GSE9")

  back <- as(se, "metaXpressStudy")
  expect_s4_class(back, "metaXpressStudy")
  expect_equal(back@accession, "GSE9")
  expect_equal(SummarizedExperiment::assay(back, "counts"),
               SummarizedExperiment::assay(s, "counts"))
})

test_that("a DESeqDataSet can enter the pipeline", {
  skip_if_not_installed("DESeq2")
  s   <- make_test_study(n_genes = 50, n_samples = 6)
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = SummarizedExperiment::assay(s, "counts"),
    colData   = S4Vectors::DataFrame(condition = factor(s$condition)),
    design    = ~ condition
  )
  s2 <- as(dds, "metaXpressStudy")
  expect_s4_class(s2, "metaXpressStudy")
  expect_equal(dim(s2), dim(s))
})
