#===========================================================
# 06_CAF_enrichment.R
# GO and Reactome enrichment of CAF candidate genes
#===========================================================

source("scripts/00_project_setup.R")

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

candidate_file <- file.path(
  table_dir,
  "CAF_candidate_genes_P001.csv"
)

enrichment_dir <- file.path(
  "results",
  "05_enrichment"
)

figure_dir <- "figures"

###Read the candidate genes
caf_candidates <- read.csv(
  candidate_file,
  row.names = 1
)

cat("Candidate genes:", nrow(caf_candidates), "\n")

###Extract gene symbols
gene_symbols <- rownames(caf_candidates)

length(gene_symbols)

###Convert to Entrez IDs
gene_map <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

###GO Biological Process
ego_bp <- enrichGO(
  gene          = gene_map$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  readable      = TRUE
)

###Inspect the results
head(as.data.frame(ego_bp))

nrow(as.data.frame(ego_bp))

# Read differential expression results
#---------------------------------------

de_results <- read.csv(
  "results/expression/GSE23066_DE_results.csv",
  row.names = 1,
  check.names = FALSE
)

cat("Genes in DE table:", nrow(de_results), "\n")
