make_meta_result <- function(mt) {
  methods::new("metaXpressResult", meta_table = mt, method = "fisher",
               n_studies = 3L)
}

# 20 up genes (G1-G20), 20 down genes (G21-G40), 260 null genes
make_signature_result <- function(seed = 1) {
  set.seed(seed)
  n <- 300
  lfc <- c(rep(2, 20), rep(-2, 20), rnorm(n - 40, sd = 0.2))
  p   <- c(rep(1e-8, 40), runif(n - 40))
  make_meta_result(data.frame(
    gene_id = paste0("G", seq_len(n)), symbol = paste0("SYM", seq_len(n)),
    meta_log2FC = lfc, meta_pvalue = p, meta_padj = p.adjust(p, "BH"),
    i_squared = 0, n_studies = 3L, stringsAsFactors = FALSE
  ))
}

# Studies where marker/module genes G1-G25 rise in cases
make_shift_studies <- function(n_studies = 3, n_genes = 300, shift = 3) {
  studies <- lapply(seq_len(n_studies), function(i) {
    set.seed(100 + i)
    n <- 8
    cond <- rep(c("control", "case"), each = n / 2)
    latent <- c(rnorm(n / 2, 0), rnorm(n / 2, shift))
    # Deep background keeps the shifted genes from distorting library sizes
    base <- matrix(rnbinom(n_genes * n, mu = 5000, size = 10), n_genes, n)
    base[1:25, ] <- round(200 * 2^(outer(rep(1, 25), latent) +
                                     matrix(rnorm(25 * n, sd = 0.2), 25, n)))
    dimnames(base) <- list(paste0("G", seq_len(n_genes)),
                           paste0("S", i, "_", seq_len(n)))
    metaXpressStudy(base, data.frame(sample_id = colnames(base),
                                     condition = cond),
                    accession = paste0("GSE", i))
  })
  names(studies) <- paste0("study", seq_len(n_studies))
  studies
}

test_that("mx_tf_activity recovers an activated and a repressed TF", {
  res <- make_signature_result()
  net <- rbind(
    data.frame(source = "TF_UP",   target = paste0("G", 1:20),  mor = 1),
    data.frame(source = "TF_DOWN", target = paste0("G", 1:20),  mor = -1),
    data.frame(source = "TF_NULL", target = paste0("G", 101:120), mor = 1)
  )
  tf <- mx_tf_activity(res, net)
  expect_s3_class(tf, "data.frame")
  expect_equal(attr(tf, "mx_type"), "tf_activity")
  expect_gt(tf$score[tf$source == "TF_UP"], 5)
  expect_lt(tf$score[tf$source == "TF_DOWN"], -5)
  expect_gt(tf$pvalue[tf$source == "TF_NULL"], 0.01)
})

test_that("mx_tf_activity reads OmniPath columns and honours id_col", {
  res <- make_signature_result()
  omni <- data.frame(source_genesymbol = "TF_A",
                     target_genesymbol = paste0("SYM", 1:20),
                     is_stimulation = "False", is_inhibition = "True")
  tf <- mx_tf_activity(res, omni, id_col = "symbol")
  expect_lt(tf$score, 0)
  expect_error(mx_tf_activity(res, omni), "id_col")
})

test_that("mx_drug_connectivity ranks the reversing drug first", {
  res  <- make_signature_result()
  sigs <- list(
    mimic    = list(up = paste0("G", 1:20),  down = paste0("G", 21:40)),
    reverser = list(up = paste0("G", 21:40), down = paste0("G", 1:20)),
    random   = list(up = paste0("G", 101:120), down = paste0("G", 201:220))
  )
  dc <- mx_drug_connectivity(res, sigs)
  expect_equal(dc$drug[1], "reverser")
  expect_lt(dc$cmap_score[dc$drug == "reverser"], 0)
  expect_gt(dc$reversal_score[dc$drug == "mimic"], 0)
  expect_error(mx_drug_connectivity(res, sigs, min_size = 50), "measured genes")
})

test_that("mx_read_signatures pairs up/down sets", {
  f <- tempfile(fileext = ".txt")
  writeLines(c("Drug A Up\t\tG1\tG2", "Drug A Down\t\tG3,1.0\tG4",
               "Drug_B_up\t\tG5", "Lonely Up\t\tG6"), f)
  sigs <- mx_read_signatures(f)
  expect_named(sigs, "Drug A")
  expect_equal(sigs[["Drug A"]]$down, c("G3", "G4"))
})

test_that("mx_cell_composition detects a shifted cell type", {
  studies <- make_shift_studies()
  markers <- list(shifted = paste0("G", 1:25), stable = paste0("G", 200:230))
  comp <- mx_cell_composition(studies, markers)
  expect_equal(attr(comp, "mx_type"), "cell_composition")
  expect_true(all(c("g_study1", "g_study2", "g_study3") %in% colnames(comp)))
  expect_gt(comp$meta_g[comp$cell_type == "shifted"], 1)
  expect_lt(comp$pvalue[comp$cell_type == "shifted"], 0.001)
  expect_gt(comp$pvalue[comp$cell_type == "stable"], 0.01)
  expect_length(attr(comp, "scores"), 3)

  df_markers <- data.frame(gene = paste0("G", 1:25), cell_type = "shifted")
  expect_equal(nrow(mx_cell_composition(studies, df_markers)), 1)
})

test_that("mx_coexpression_modules finds the co-regulated module", {
  studies <- make_shift_studies()
  res <- make_signature_result()
  net <- mx_coexpression_modules(studies, n_genes = 300, min_module_size = 10,
                                 meta_result = res)
  expect_named(net, c("modules", "module_stats", "eigengenes"))
  m <- unique(net$modules$module[net$modules$gene_id %in% paste0("G", 1:25)])
  expect_length(m, 1)
  top <- net$module_stats[net$module_stats$module == m, ]
  expect_gt(top$meta_g, 1)
  expect_gt(top$pct_meta_up, 50)
  expect_match(top$hub_genes, "^G")
  expect_error(mx_coexpression_modules(studies[1]), "at least two")
})

test_that("mx_downstream_plot returns ggplots for every result type", {
  res <- make_signature_result()
  studies <- make_shift_studies()
  tf <- mx_tf_activity(res, data.frame(source = "TF", target = paste0("G", 1:20)))
  expect_s3_class(mx_downstream_plot(tf), "ggplot")
  comp <- mx_cell_composition(studies, list(a = paste0("G", 1:25)))
  expect_s3_class(mx_downstream_plot(comp), "ggplot")
  net <- mx_coexpression_modules(studies, n_genes = 300, min_module_size = 10)
  expect_s3_class(mx_downstream_plot(net), "ggplot")
  expect_error(mx_downstream_plot(data.frame(a = 1)), "downstream")
})
