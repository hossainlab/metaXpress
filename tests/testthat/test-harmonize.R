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

make_id_study <- function(ids, counts) {
  m <- matrix(counts, nrow = length(ids), byrow = TRUE,
              dimnames = list(ids, paste0("S", 1:4)))
  metaXpressStudy(m, data.frame(condition = rep(c("control", "case"),
                                                each = 2)),
                  accession = "X")
}

test_that("ID helpers strip versions and detect type by majority", {
  expect_equal(.clean_gene_ids(c("ENSG00000141510.17", "ENSG00000002586.20_PAR_Y",
                                 "TP53", "ENSMUSG00000059552.13")),
               c("ENSG00000141510", "ENSG00000002586", "TP53",
                 "ENSMUSG00000059552"))
  expect_equal(.detect_id_type(c("7157", "672", "TP53")), "ENTREZID")
  expect_equal(.detect_id_type(c("ENSG00000141510.1", "ENSG00000012048", "X")),
               "ENSEMBL")
  expect_equal(.detect_id_type(c("TP53", "BRCA1", "1")), "SYMBOL")
})

test_that("mx_reannotate maps with a user table and collapses duplicates", {
  s <- make_id_study(c("ENSG1.1", "ENSG2.4", "ENSG3.1", "ENSG4.2"),
                     c(10, 10, 10, 10,  1, 1, 1, 1,  5, 5, 5, 5,  2, 2, 2, 2))
  map <- data.frame(from = c("ENSG1", "ENSG2", "ENSG3"),
                    to = c("GENEA", "GENEA", "GENEB"))
  out <- suppressMessages(mx_reannotate(s, mapping = map, from_id = "ENSEMBL"))
  expect_equal(rownames(out), c("GENEA", "GENEB"))
  # max_mean keeps the more expressed of the two GENEA rows
  expect_equal(unname(SummarizedExperiment::assay(out)["GENEA", 1]), 10)
  expect_equal(SummarizedExperiment::rowData(out)$source_id,
               c("ENSG1.1", "ENSG3.1"))

  summed <- suppressMessages(mx_reannotate(s, mapping = map,
                                           duplicates = "sum"))
  expect_equal(unname(SummarizedExperiment::assay(summed)["GENEA", 1]), 11)
  dropped <- suppressMessages(mx_reannotate(s, mapping = map,
                                            duplicates = "drop"))
  expect_equal(rownames(dropped), "GENEB")

  lst <- suppressMessages(mx_reannotate(list(a = s), mapping = map))
  expect_named(lst, "a")
})

test_that("mx_reannotate maps real Ensembl IDs with org.Hs.eg.db", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("airway")
  data("airway", package = "airway", envir = environment())
  st <- mx_from_se(airway[1:2000, ], condition_col = "dex",
                   accession = "GSE52778")
  rownames(st) <- paste0(rownames(st), ".5")         # versioned IDs
  out <- suppressMessages(mx_reannotate(st, org = "Homo sapiens"))
  expect_false(anyDuplicated(rownames(out)) > 0)
  expect_true("TSPAN6" %in% rownames(out))           # ENSG00000000003
  expect_gt(nrow(out), 1000)
  expect_true(all(grepl("^ENSG", SummarizedExperiment::rowData(out)$source_id)))
})

test_that("mx_reannotate errors helpfully for unsupported organisms", {
  s <- make_id_study(paste0("G", 1:2), rep(1, 8))
  expect_error(mx_reannotate(s, org = "Danio rerio"), "mapping")
})
