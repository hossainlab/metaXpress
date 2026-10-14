test_that("mx_meta returns a metaXpressResult object", {
  de <- make_test_de_results()
  result <- mx_meta(de, method = "random_effects")
  expect_s4_class(result, "metaXpressResult")
})

test_that("meta_table has all required columns", {
  de <- make_test_de_results()
  result <- mx_meta(de, method = "fisher")
  required <- c("gene_id", "meta_log2FC", "meta_pvalue", "meta_padj",
                 "i_squared", "n_studies", "direction_consistency")
  expect_true(all(required %in% colnames(result@meta_table)))
})

test_that("meta_padj values are in [0, 1]", {
  de <- make_test_de_results()
  for (method in c("fisher", "stouffer", "random_effects", "fixed_effects")) {
    result <- mx_meta(de, method = method)
    padj <- result@meta_table$meta_padj
    expect_true(all(padj >= 0 & padj <= 1, na.rm = TRUE),
                info = paste("method:", method))
  }
})

test_that("mx_meta n_studies slot matches input length", {
  de     <- make_test_de_results(n_studies = 3)
  result <- mx_meta(de, method = "random_effects")
  expect_equal(result@n_studies, 3L)
})

test_that("mx_meta warns with fewer than 3 studies for random_effects", {
  de <- make_test_de_results(n_studies = 2)
  expect_warning(mx_meta(de, method = "random_effects"),
                 "fewer than 3 studies")
})

test_that("mx_meta errors on invalid method", {
  de <- make_test_de_results()
  expect_error(mx_meta(de, method = "not_a_method"))
})

test_that("mx_meta errors on missing required columns", {
  de <- make_test_de_results()
  de[[1]]$padj <- NULL
  expect_error(mx_meta(de, method = "fisher"), "missing required columns")
})

test_that("mx_meta min_studies filters genes correctly", {
  # 3 studies; remove GENE1 only from study 3 → it appears in 2/3 studies
  de <- make_test_de_results(n_genes = 50, n_studies = 3)
  de[[3]] <- de[[3]][de[[3]]$gene_id != "GENE1", ]

  result_strict <- mx_meta(de, method = "fisher", min_studies = 3)
  result_loose  <- mx_meta(de, method = "fisher", min_studies = 2)

  expect_false("GENE1" %in% result_strict@meta_table$gene_id)
  expect_true("GENE1"  %in% result_loose@meta_table$gene_id)
})

test_that("mx_meta direction_consistency is in [0, 1]", {
  de <- make_test_de_results()
  result <- mx_meta(de, method = "stouffer")
  dc <- result@meta_table$direction_consistency
  expect_true(all(dc >= 0 & dc <= 1, na.rm = TRUE))
})

test_that("mx_heterogeneity returns data.frame with correct columns", {
  de  <- make_test_de_results()
  het <- mx_heterogeneity(de)
  expect_true(is.data.frame(het))
  expect_true(all(c("gene_id", "Q", "I_sq", "tau_sq", "p_heterogeneity") %in%
                    colnames(het)))
})

test_that("mx_sensitivity returns a list of length n-1", {
  de  <- make_test_de_results(n_genes = 50, n_studies = 4)
  loo <- mx_sensitivity(de, method = "fisher")
  expect_equal(length(loo), 4)
  expect_true(all(vapply(loo, is, logical(1), "metaXpressResult")))
})

test_that("mx_sensitivity errors with fewer than 3 studies", {
  de <- make_test_de_results(n_studies = 2)
  expect_error(mx_sensitivity(de), "at least 3 studies")
})

# --- Statistical correctness against reference implementations ------------

make_het_de <- function(G = 60, K = 5, seed = 7) {
  set.seed(seed)
  de <- lapply(seq_len(K), function(s) {
    lfc <- rnorm(G, rep(c(1, 0), c(10, G - 10)) +
                   rnorm(G, 0, 0.4) * (s %% 2), 0.3)
    se  <- runif(G, 0.1, 0.6)
    p   <- 2 * pnorm(-abs(lfc / se))
    d <- data.frame(gene_id = paste0("G", seq_len(G)), log2FC = lfc,
                    lfcSE = se, pvalue = p, padj = p.adjust(p, "BH"),
                    n_case = 4 + s, n_control = 4)
    d[-sample(G, 5), ]
  })
  names(de) <- LETTERS[seq_len(K)]
  de
}

test_that("random/fixed effects match metafor::rma for all estimators", {
  skip_if_not_installed("metafor")
  de <- make_het_de()
  for (tm in c("DL", "REML", "PM")) for (hk in c(FALSE, TRUE)) {
    mt <- suppressMessages(mx_meta(de, tau_method = tm, hksj = hk))@meta_table
    for (i in c(1, 5, 20, 40)) {
      g  <- mt$gene_id[i]
      y  <- vapply(de, function(d) d$log2FC[match(g, d$gene_id)], numeric(1))
      s  <- vapply(de, function(d) d$lfcSE[match(g, d$gene_id)], numeric(1))
      ok <- !is.na(y)
      ref <- metafor::rma(yi = y[ok], sei = s[ok], method = tm,
                          test = if (hk) "knha" else "z",
                          control = list(threshold = 1e-12, tol = 1e-12, maxiter = 5000,
                                         stepadj = 0.5))
      info <- paste(tm, hk, g)
      expect_equal(mt$meta_log2FC[i], unname(ref$b[1, 1]), tolerance = 1e-7,
                   info = info)
      expect_equal(mt$meta_se[i], ref$se, tolerance = 1e-7, info = info)
      expect_equal(mt$meta_pvalue[i], ref$pval, tolerance = 1e-6, info = info)
      expect_equal(mt$tau2[i], ref$tau2, tolerance = 1e-7, info = info)
      expect_equal(mt$q_stat[i], ref$QE, tolerance = 1e-8, info = info)
      expect_equal(mt$ci_lower[i], ref$ci.lb, tolerance = 1e-7, info = info)
    }
  }
  fe <- suppressMessages(mx_meta(de, method = "fixed_effects"))@meta_table
  g  <- fe$gene_id[3]
  y  <- vapply(de, function(d) d$log2FC[match(g, d$gene_id)], numeric(1))
  s  <- vapply(de, function(d) d$lfcSE[match(g, d$gene_id)], numeric(1))
  ok <- !is.na(y)
  ref <- metafor::rma(yi = y[ok], sei = s[ok], method = "FE")
  expect_equal(fe$meta_log2FC[3], unname(ref$b[1, 1]), tolerance = 1e-8)
  expect_equal(fe$meta_pvalue[3], ref$pval, tolerance = 1e-8)
})

test_that("prediction interval uses t with k - 2 df (Higgins 2009)", {
  de <- make_het_de()
  mt <- suppressMessages(mx_meta(de, tau_method = "DL"))@meta_table
  i  <- which(mt$n_studies >= 3)[1]
  half <- qt(0.975, mt$n_studies[i] - 2) * sqrt(mt$tau2[i] + mt$meta_se[i]^2)
  expect_equal(mt$pi_upper[i], mt$meta_log2FC[i] + half)
})

test_that("p-value combination methods match closed forms", {
  de <- make_het_de(K = 3)
  g  <- "G1"
  p  <- vapply(de, function(d) d$pvalue[match(g, d$gene_id)], numeric(1))
  l  <- vapply(de, function(d) d$log2FC[match(g, d$gene_id)], numeric(1))
  ok <- !is.na(p); p <- p[ok]; l <- l[ok]; k <- length(p)
  get <- function(...) {
    mt <- suppressMessages(mx_meta(de, ...))@meta_table
    mt$meta_pvalue[mt$gene_id == g]
  }
  expect_equal(get(method = "fisher", directional = FALSE),
               pchisq(-2 * sum(log(p)), 2 * k, lower.tail = FALSE))
  z <- qnorm(ifelse(l > 0, p / 2, 1 - p / 2), lower.tail = FALSE)
  expect_equal(get(method = "stouffer"), 2 * pnorm(-abs(sum(z) / sqrt(k))))
  n <- vapply(de, function(d) d$n_case[1] + d$n_control[1], numeric(1))[ok]
  w <- sqrt(n)
  expect_equal(get(method = "inverse_normal"),
               2 * pnorm(-abs(sum(w * z) / sqrt(sum(w^2)))))
  expect_equal(get(method = "minp"), 1 - (1 - min(p))^k)
  expect_equal(get(method = "maxp"), max(p)^k)
  expect_equal(get(method = "rop", r = 2), unname(pbeta(sort(p)[2], 2, k - 1)))
})

test_that("directional combination does not accumulate opposite effects", {
  opp <- list(
    a = data.frame(gene_id = "X", log2FC = 2, lfcSE = 0.5,
                   pvalue = 1e-4, padj = 1e-3),
    b = data.frame(gene_id = "X", log2FC = -2, lfcSE = 0.5,
                   pvalue = 1e-4, padj = 1e-3))
  stouffer <- suppressMessages(mx_meta(opp, "stouffer"))@meta_table
  expect_equal(stouffer$meta_pvalue, 1)
  dir <- suppressMessages(mx_meta(opp, "fisher"))@meta_table
  undir <- suppressMessages(mx_meta(opp, "fisher",
                                    directional = FALSE))@meta_table
  expect_gt(dir$meta_pvalue, undir$meta_pvalue)
  expect_equal(dir$direction_consistency, 0.5)
})

test_that("awmeta is rejected with an informative error", {
  expect_error(mx_meta(make_test_de_results(), method = "awmeta"),
               "removed")
})

test_that("effect-size methods require lfcSE", {
  de <- make_test_de_results()
  de[[1]]$lfcSE <- NULL
  expect_error(suppressMessages(mx_meta(de, "random_effects")), "lfcSE")
  expect_s4_class(suppressMessages(mx_meta(de, "stouffer")),
                  "metaXpressResult")
})

test_that("mx_sensitivity attaches a robustness table", {
  de  <- make_het_de(K = 4)
  loo <- suppressMessages(mx_sensitivity(de))
  rob <- attr(loo, "robustness")
  expect_true(all(c("gene_id", "n_loo_sig", "frac_loo_sig") %in% names(rob)))
  expect_true(all(rob$frac_loo_sig >= 0 & rob$frac_loo_sig <= 1))
})

test_that("mx_study_concordance flags a flipped contrast", {
  set.seed(1)
  base <- rnorm(300)
  de <- lapply(c(1, 1, -1), function(dir) {
    lfc <- dir * base + rnorm(300, 0, 0.3)
    p <- 2 * pnorm(-abs(lfc / 0.3))
    data.frame(gene_id = paste0("G", 1:300), log2FC = lfc, pvalue = p,
               padj = p.adjust(p, "BH"))
  })
  names(de) <- c("A", "B", "C")
  expect_warning(cc <- mx_study_concordance(de), "contrast")
  expect_true(cc$correlation[cc$study_1 == "A" & cc$study_2 == "B"] > 0.8)
  expect_match(cc$flag[cc$study_2 == "C"][1], "NEGATIVE")
})
