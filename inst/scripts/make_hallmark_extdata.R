# ==============================================================================
# Build inst/extdata/hallmark_gene_sets.rds
# ==============================================================================
# The 50 MSigDB Hallmark gene sets (Liberzon et al. 2015, Cell Systems 1:417)
# with human gene symbols, shipped so that vignettes and examples run offline.
# MSigDB Hallmark gene sets are distributed under the Creative Commons
# Attribution 4.0 license (https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp).
#
# Usage (from the package root):  Rscript inst/scripts/make_hallmark_extdata.R
# ==============================================================================
suppressPackageStartupMessages(devtools::load_all(quiet = TRUE))
hallmarks <- mx_gene_sets("Hallmarks", species = "Homo sapiens",
                          id_type = "symbol")
attr(hallmarks, "source") <- paste(
  "MSigDB Hallmark collection, version", attr(hallmarks, "version"),
  "(msigdbr", as.character(packageVersion("msigdbr")), "); CC BY 4.0;",
  "Liberzon et al. (2015) Cell Systems 1:417, doi:10.1016/j.cels.2015.12.004")
saveRDS(hallmarks, "inst/extdata/hallmark_gene_sets.rds", compress = "xz")
message(length(hallmarks), " gene sets; ",
        round(file.size("inst/extdata/hallmark_gene_sets.rds") / 1024), " KB")
