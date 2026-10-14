# ==============================================================================
# Build inst/extdata/case_study_downstream.rds for the downstream_biology vignette
# ==============================================================================
# Compact, offline inputs so the vignette can run Module 9 on the real
# three-cohort PDAC meta-analysis. Requires scripts/fetch_case_study_data.R and
# scripts/run_case_study_applications.R to have been run first.
#
# Contents:
#   meta_result  genome-wide random-effects meta-analysis (29,560 genes, with
#                gene symbols)
#   studies      the three cohorts restricted to MCP-counter markers plus the
#                1,500 most variable genes; colData$lib_size keeps true depth,
#                rowData$symbol holds HGNC symbols
#   markers      MCP-counter marker genes (Becht et al. 2016), Entrez IDs
#   collectri    CollecTRI TF-target network (source, target, mor)
#   lincs        200 LINCS L1000 consensus signatures: the 100 strongest
#                reversers, the 50 strongest mimickers and 50 random others
#                from the full 5,425-compound ranking (demo subset only)
# ==============================================================================

suppressPackageStartupMessages(devtools::load_all(quiet = TRUE))
in_dir <- "geo_cache/case_study"
set.seed(2026)

obj <- readRDS(file.path(in_dir, "pipeline_objects.rds"))

# Genome-wide meta result, trimmed to the columns the downstream functions use
meta_result <- obj$meta_res
keep_cols <- c("gene_id", "symbol", "meta_log2FC", "meta_se", "meta_pvalue",
               "meta_padj", "i_squared", "n_studies", "direction_consistency")
meta_result@meta_table <- meta_result@meta_table[, keep_cols]
num <- vapply(meta_result@meta_table, is.double, logical(1))
meta_result@meta_table[num] <- lapply(meta_result@meta_table[num], signif, 4)
meta_result@heterogeneity <- data.frame()

mcp <- read.delim(file.path(in_dir, "mcpcounter_genes.txt"),
                  colClasses = "character", check.names = FALSE)
markers <- data.frame(gene = mcp$ENTREZID, cell_type = mcp$`Cell population`)

# Study subsets: markers + most variable genes, keeping true library sizes
lcpm <- lapply(obj$studies, function(s) {
  x <- SummarizedExperiment::assay(s, "counts")
  log2(t(t(x) / colSums(x)) * 1e6 + 1)
})
v <- rowMeans(vapply(lcpm, function(x) apply(x, 1, var), numeric(nrow(lcpm[[1]]))))
top <- names(sort(v, decreasing = TRUE))[1:1500]
genes <- union(intersect(markers$gene, rownames(obj$studies[[1]])), top)
gene_info <- read.delim(file.path(in_dir, "Homo_sapiens.gene_info.gz"),
                        quote = "", colClasses = "character")
symbol_of <- setNames(gene_info$Symbol, gene_info$GeneID)
studies <- lapply(obj$studies, function(s) {
  s$lib_size <- colSums(SummarizedExperiment::assay(s, "counts"))
  s <- s[genes, ]
  SummarizedExperiment::rowData(s)$symbol <- unname(symbol_of[genes])
  s@de_result <- data.frame()
  s
})

collectri <- metaXpress:::.standardize_network(
  read.delim(file.path(in_dir, "collectri.tsv"))
)

sigs  <- mx_read_signatures(file.path(in_dir, "lincs_l1000_chem_consensus.txt"))
drugs <- read.csv("case_study_results/app_drug_connectivity.csv")
pick  <- c(head(drugs$drug, 100), tail(drugs$drug, 50))
pick  <- c(pick, sample(setdiff(drugs$drug, pick), 50))
lincs <- sigs[pick]

out <- list(meta_result = meta_result, studies = studies, markers = markers,
            collectri = collectri, lincs = lincs)
saveRDS(out, "inst/extdata/case_study_downstream.rds", compress = "xz")
message("Saved inst/extdata/case_study_downstream.rds (",
        round(file.size("inst/extdata/case_study_downstream.rds") / 1024), " KB)")
