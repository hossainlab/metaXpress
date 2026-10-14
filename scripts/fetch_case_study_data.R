# ==============================================================================
# Rebuild the PDAC Module 9 demo input data from public sources
# ==============================================================================
# Downloads NCBI-processed raw counts (GRCh38.p13, Entrez IDs) for the three
# GEO PDAC cohorts, the matching NCBI gene annotation, and the external resources
# used by the downstream biology analyses:
#   - CollecTRI TF-target regulons (OmniPath; Mueller-Dott et al. 2023)
#   - MCP-counter cell-type marker genes (Becht et al. 2016)
#   - LINCS L1000 chemical perturbation consensus signatures (Enrichr)
# Sample condition labels are read from the GEO series metadata by
# scripts/run_case_study_applications.R.
# Output: geo_cache/case_study/ (git-ignored)
# ==============================================================================

out_dir <- "geo_cache/case_study"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = 1200)

fetch <- function(url, dest) {
  dest <- file.path(out_dir, dest)
  if (!file.exists(dest)) {
    message("Downloading ", basename(dest))
    utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  }
  dest
}

geo <- "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&format=file"
for (acc in c("GSE130688", "GSE136569", "GSE171485")) {
  fetch(sprintf("%s&acc=%s&file=%s_raw_counts_GRCh38.p13_NCBI.tsv.gz", geo, acc, acc),
        sprintf("%s_raw_counts.tsv.gz", acc))
}
fetch(paste0("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/",
             "Homo_sapiens.gene_info.gz"),
      "Homo_sapiens.gene_info.gz")

fetch(paste0("https://omnipathdb.org/interactions?datasets=collectri",
             "&genesymbols=yes&organisms=9606"),
      "collectri.tsv")

fetch("https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt",
      "mcpcounter_genes.txt")

fetch(paste0("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text",
             "&libraryName=LINCS_L1000_Chem_Pert_Consensus_Sigs"),
      "lincs_l1000_chem_consensus.txt")

message("Done: ", out_dir)
