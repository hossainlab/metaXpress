## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  message  = FALSE,
  warning  = FALSE,
  eval     = FALSE,  # All chunks are illustrative; run interactively
  purl     = FALSE   # Exclude from R CMD check vignette code extraction
)


## ----load-package, eval=FALSE-------------------------------------------------
# library(metaXpress)


## ----load-data, eval=FALSE----------------------------------------------------
# data(metaXpress_example)
# studies <- metaXpress_example
# lapply(studies, show)


## ----filter-studies, eval=FALSE-----------------------------------------------
# studies <- mx_filter_studies(studies, qc_threshold = 7)
# message("Studies passing QC: ", length(studies))


## ----align-genes, eval=FALSE--------------------------------------------------
# studies <- mx_align_genes(studies)


## ----de-analysis, eval=FALSE--------------------------------------------------
# studies <- mx_de_all(studies, method = "DESeq2", formula = ~ condition)
# mx_de_summary(studies)


## ----meta-analysis, eval=FALSE------------------------------------------------
# de_results <- lapply(studies, function(s) s@de_result)
# meta_result <- mx_meta(de_results, method = "random_effects")
# meta_result


## ----volcano, fig.width=7, fig.height=5, eval=FALSE---------------------------
# mx_volcano(meta_result, label_top = 10)


## ----forest, fig.width=6, fig.height=4, eval=FALSE----------------------------
# # Forest plot for the top gene
# top_gene <- meta_result@meta_table$gene_id[
#   which.min(meta_result@meta_table$meta_padj)]
# mx_forest(top_gene, de_results, meta_result)


## ----heterogeneity, fig.width=6, fig.height=4, eval=FALSE---------------------
# mx_heterogeneity_plot(meta_result)


## ----export, eval=FALSE-------------------------------------------------------
# mx_export(meta_result, format = "csv", output_dir = tempdir())


## ----session-info, eval=FALSE-------------------------------------------------
# mx_session_info()$session_info

