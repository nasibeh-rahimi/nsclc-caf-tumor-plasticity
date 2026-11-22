############################################################
# 02_probe_annotation.R
# Affymetrix probe ID -> gene symbol
# collapse multiple probes per gene
############################################################

source("scripts/00_project_setup.R")

library(AnnotationDbi)
library(hgu133plus2.db)
library(dplyr)


# ============================
# 1. Load RMA matrix
# ============================

expr_rma <- readRDS(
  "data/processed/expr_GSE33532_RMA.rds"
)

cat("RMA matrix:")
print(dim(expr_rma))


# ============================
# 2. Probe annotation
# ============================

probe_ids <- rownames(expr_rma)


gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = probe_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)


# ============================
# 3. Add gene names
# ============================

expr_df <- as.data.frame(expr_rma)

expr_df$gene <- gene_symbols


# remove probes without genes

expr_df <- expr_df[
  !is.na(expr_df$gene) &
    expr_df$gene != "",
]


cat(
  "After removing unmapped probes:",
  nrow(expr_df),
  "\n"
)



# ============================
# 4. Collapse duplicate genes
# keep highest variance probe
# ============================


expression_matrix <- expr_df[
  ,
  !colnames(expr_df) %in% "gene"
]


probe_variance <- apply(
  expression_matrix,
  1,
  var,
  na.rm = TRUE
)


expr_df$variance <- probe_variance



expr_df <- expr_df[
  order(
    expr_df$gene,
    -expr_df$variance
  ),
]


expr_df <- expr_df[
  !duplicated(expr_df$gene),
]



# ============================
# 5. Create final gene matrix
# ============================

expr_gene <- expr_df[
  ,
  !colnames(expr_df) %in% c("gene","variance")
]


rownames(expr_gene) <- expr_df$gene


expr_gene <- as.matrix(expr_gene)



cat(
  "FINAL gene matrix:",
  dim(expr_gene),
  "\n"
)



# ============================
# 6. Save
# ============================


saveRDS(
  expr_gene,
  "data/processed/expr_GSE33532_gene_symbol.rds"
)


write.csv(
  expr_gene,
  "results/expression/expr_GSE33532_gene_symbol.csv"
)


cat("Annotation completed\n")
