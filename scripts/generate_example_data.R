library(metaXpress)

set.seed(42)

# Generate genes
n_genes <- 500
genes <- paste0("GENE", sprintf("%04d", 1:n_genes))
sig_genes <- genes[1:30] # Signal in first 30 genes

studies <- list()

# study1: 500 genes, 6 samples (3 case, 3 control)
n_control1 <- 3
n_case1 <- 3
n_samples1 <- 6
base_expr1 <- rnbinom(500, mu = 500, size = 2)
counts1 <- matrix(rnbinom(500 * 6, mu = base_expr1, size = 2), nrow = 500, ncol = 6)
rownames(counts1) <- genes
colnames(counts1) <- paste0("Study1_S", 1:6)
counts1[1:30, 4:6] <- counts1[1:30, 4:6] * matrix(runif(30 * 3, 2, 5), nrow=30)
counts1 <- round(counts1)

meta1 <- data.frame(
  sample_id = colnames(counts1),
  condition = c(rep("control", 3), rep("case", 3)),
  library_type = sample(c("polyA", "rRNA_depleted"), 6, replace=TRUE),
  stringsAsFactors = FALSE
)
attr(counts1, "gene_length") <- sample(500:5000, 500, replace=TRUE)

studies$study1 <- new("metaXpressStudy", counts = counts1, metadata = meta1,
                      accession = "SIM001", organism = "human", qc_score = 8.5,
                      de_result = data.frame())

# study2: 500 genes, 6 samples (3 case, 3 control)
base_expr2 <- rnbinom(500, mu = 400, size = 2)
counts2 <- matrix(rnbinom(500 * 6, mu = base_expr2, size = 2), nrow = 500, ncol = 6)
rownames(counts2) <- genes
colnames(counts2) <- paste0("Study2_S", 1:6)
counts2[1:30, 4:6] <- counts2[1:30, 4:6] * matrix(runif(30 * 3, 2, 5), nrow=30)
counts2 <- round(counts2)

meta2 <- data.frame(
  sample_id = colnames(counts2),
  condition = c(rep("control", 3), rep("case", 3)),
  library_type = sample(c("polyA", "rRNA_depleted"), 6, replace=TRUE),
  stringsAsFactors = FALSE
)
attr(counts2, "gene_length") <- sample(500:5000, 500, replace=TRUE)

studies$study2 <- new("metaXpressStudy", counts = counts2, metadata = meta2,
                      accession = "SIM002", organism = "human", qc_score = 9.0,
                      de_result = data.frame())

# study3: 480 genes, 8 samples (4 case, 4 control)
# remove 20 random non-significant genes to simulate partial overlap
genes3 <- genes[-(31:50)]
base_expr3 <- rnbinom(480, mu = 600, size = 2)
counts3 <- matrix(rnbinom(480 * 8, mu = base_expr3, size = 2), nrow = 480, ncol = 8)
rownames(counts3) <- genes3
colnames(counts3) <- paste0("Study3_S", 1:8)
counts3[1:30, 5:8] <- counts3[1:30, 5:8] * matrix(runif(30 * 4, 2, 5), nrow=30)
counts3 <- round(counts3)

meta3 <- data.frame(
  sample_id = colnames(counts3),
  condition = c(rep("control", 4), rep("case", 4)),
  library_type = sample(c("polyA", "rRNA_depleted"), 8, replace=TRUE),
  stringsAsFactors = FALSE
)
attr(counts3, "gene_length") <- sample(500:5000, 480, replace=TRUE)

studies$study3 <- new("metaXpressStudy", counts = counts3, metadata = meta3,
                      accession = "SIM003", organism = "human", qc_score = 7.8,
                      de_result = data.frame())

metaXpress_example <- studies
dir.create("data", showWarnings = FALSE)
save(metaXpress_example, file = "data/metaXpress_example.rda", compress = "xz")
cat("Successfully created data/metaXpress_example.rda\n")
