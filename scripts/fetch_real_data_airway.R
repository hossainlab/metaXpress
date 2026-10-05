devtools::load_all(".")
library(airway)

message("Generating example data from real airway RNA-seq dataset...")
data(airway, package="airway")

counts_all <- as.matrix(assay(airway, "counts"))
meta_all <- as.data.frame(colData(airway))

# Subsetting to the top 500 variable genes for file size
vars <- apply(counts_all, 1, var)
top_genes <- names(sort(vars, decreasing=TRUE))[1:500]
counts_all <- counts_all[top_genes, ]

# Create Study 1: Cell lines N61311, N052611
idx1 <- meta_all$cell %in% c("N61311", "N052611")
meta1 <- meta_all[idx1, ]
counts1 <- counts_all[, idx1]

df1 <- data.frame(
  sample_id = rownames(meta1),
  condition = ifelse(meta1$dex == "trt", "case", "control"),
  library_type = "polyA",
  stringsAsFactors = FALSE
)
attr(counts1, "gene_length") <- rep(2000, 500)

study1 <- new("metaXpressStudy",
              counts = counts1,
              metadata = df1,
              accession = "GSE52778_Part1",
              organism = "human",
              qc_score = 9.5,
              de_result = data.frame())

# Create Study 2: Cell lines N080611, N061011
idx2 <- meta_all$cell %in% c("N080611", "N061011")
meta2 <- meta_all[idx2, ]
counts2 <- counts_all[, idx2]

df2 <- data.frame(
  sample_id = rownames(meta2),
  condition = ifelse(meta2$dex == "trt", "case", "control"),
  library_type = "polyA",
  stringsAsFactors = FALSE
)
attr(counts2, "gene_length") <- rep(2000, 500)

study2 <- new("metaXpressStudy",
              counts = counts2,
              metadata = df2,
              accession = "GSE52778_Part2",
              organism = "human",
              qc_score = 9.0,
              de_result = data.frame())

# Create Study 3: Slightly modified Study 2 to act as a 3rd cohort
# To make meta-analysis work across 3 studies
counts3 <- counts2
# Add some slight noise
noise <- matrix(rnbinom(length(counts3), mu = 5, size = 2), nrow=nrow(counts3))
counts3 <- counts3 + noise
df3 <- df2

attr(counts3, "gene_length") <- rep(2000, 500)

study3 <- new("metaXpressStudy",
              counts = counts3,
              metadata = df3,
              accession = "GSE52778_Part3",
              organism = "human",
              qc_score = 8.5,
              de_result = data.frame())

metaXpress_example <- list(study1 = study1, study2 = study2, study3 = study3)

save(metaXpress_example, file = "data/metaXpress_example.rda", compress = "xz")
message("Real airway data successfully saved to data/metaXpress_example.rda")
