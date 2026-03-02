# ============================================================================
# Shared test fixtures — sourced automatically by testthat before all tests
# ============================================================================

#' Build a minimal metaXpressStudy for testing (no network calls)
make_test_study <- function(n_genes = 200, n_samples = 6,
                             accession = "TEST001",
                             seed = 42) {
  set.seed(seed)
  gene_ids  <- paste0("GENE", seq_len(n_genes))
  sample_ids <- paste0("S", seq_len(n_samples))

  counts <- matrix(
    as.integer(abs(rnorm(n_genes * n_samples, mean = 500, sd = 200))),
    nrow = n_genes, ncol = n_samples,
    dimnames = list(gene_ids, sample_ids)
  )
  # Inflate counts for realistic library sizes (>10M)
  counts <- counts * 20000L

  metadata <- data.frame(
    sample_id = sample_ids,
    condition = rep(c("control", "case"), each = n_samples / 2),
    stringsAsFactors = FALSE
  )

  methods::new("metaXpressStudy",
               counts    = counts,
               metadata  = metadata,
               accession = accession,
               organism  = "Homo sapiens",
               qc_score  = NA_real_,
               de_result = data.frame())
}

#' Build a minimal list of DE data.frames for testing
make_test_de_results <- function(n_genes = 100, n_studies = 3, seed = 42) {
  set.seed(seed)
  lapply(seq_len(n_studies), function(i) {
    pv <- runif(n_genes)
    data.frame(
      gene_id  = paste0("GENE", seq_len(n_genes)),
      log2FC   = rnorm(n_genes, mean = 0, sd = 1.5),
      pvalue   = pv,
      padj     = p.adjust(pv, method = "BH"),
      baseMean = abs(rnorm(n_genes, mean = 200, sd = 50)),
      lfcSE    = abs(rnorm(n_genes, mean = 0.3, sd = 0.1)),
      method   = "DESeq2",
      stringsAsFactors = FALSE
    )
  })
}
