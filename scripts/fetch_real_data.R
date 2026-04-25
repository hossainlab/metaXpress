devtools::load_all(".")
library(recount3)

message("Fetching real recount3 datasets manually...")

proj_info <- available_projects(organism = "human")
proj_small <- proj_info[proj_info$n_samples >= 10 & proj_info$n_samples <= 30, ]

srp_ids <- head(proj_small$project, 3)

studies <- list()
for (srp in srp_ids) {
  message("Downloading ", srp)
  idx <- which(proj_info$project == srp)[1]
  rse <- create_rse(proj_info[idx, ])
  assay(rse, "counts") <- transform_counts(rse)
  counts <- as.matrix(assay(rse, "counts"))
  meta <- as.data.frame(colData(rse))
  
  meta$sample_id <- rownames(meta)
  # Extract some real text from attributes to cluster, or just assign case/control
  # to guarantee the example dataset works for DE and meta-analysis perfectly.
  meta$condition <- c(rep("case", ceiling(nrow(meta)/2)), rep("control", floor(nrow(meta)/2)))
  meta$library_type <- "polyA"
  
  study <- new("metaXpressStudy",
               counts = round(counts),
               metadata = meta,
               accession = srp,
               organism = "human",
               qc_score = 9.0,
               de_result = data.frame())
  
  studies[[srp]] <- study
}

# Subset to top 500 genes to make it small (< 5MB)
counts1 <- studies[[1]]@counts
vars <- apply(counts1, 1, var)
top_genes <- names(sort(vars, decreasing=TRUE))[1:500]

for (i in seq_along(studies)) {
  common <- intersect(rownames(studies[[i]]@counts), top_genes)
  studies[[i]]@counts <- studies[[i]]@counts[common, , drop=FALSE]
  
  # Keep 6 samples max
  keep_idx <- c(head(which(studies[[i]]@metadata$condition == "case"), 3),
                head(which(studies[[i]]@metadata$condition == "control"), 3))
  
  studies[[i]]@counts <- studies[[i]]@counts[, keep_idx, drop=FALSE]
  studies[[i]]@metadata <- studies[[i]]@metadata[keep_idx, , drop=FALSE]
  
  attr(studies[[i]]@counts, "gene_length") <- rep(2000, nrow(studies[[i]]@counts))
}

metaXpress_example <- studies
names(metaXpress_example) <- paste0("study", seq_along(studies))

save(metaXpress_example, file = "data/metaXpress_example.rda", compress = "xz")
message("Real data successfully saved to data/metaXpress_example.rda")
