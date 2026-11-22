############################################################
# 04_GSE23066_annotation.R
# Probe annotation and duplicate collapsing
############################################################

source("scripts/00_project_setup.R")

# ============================
# 1. Load RMA matrix
# ============================

expr_rma <- readRDS(
  "data/processed/expr_GSE23066_RMA.rds"
)

cat("RMA matrix dimensions:\n")
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
# 3. Build annotation table
# ============================

expr_df <- as.data.frame(expr_rma)
expr_df$gene <- gene_symbols

# Remove probes without gene symbols
expr_df <- expr_df[
  !is.na(expr_df$gene) &
    expr_df$gene != "",
]

cat("After removing unmapped probes:\n")
print(nrow(expr_df))

# ============================
# 4. Collapse duplicate genes
# ============================

expression_matrix <- expr_df[
  ,
  !colnames(expr_df) %in% "gene"
]

expr_df$variance <- apply(
  expression_matrix,
  1,
  var,
  na.rm = TRUE
)

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
# 5. Final gene matrix
# ============================

expr_gene <- expr_df[
  ,
  !colnames(expr_df) %in% c("gene", "variance")
]

rownames(expr_gene) <- expr_df$gene

expr_gene <- as.matrix(expr_gene)

cat("Final gene matrix dimensions:\n")
print(dim(expr_gene))

# ============================
# 6. Save outputs
# ============================

saveRDS(
  expr_gene,
  "data/processed/expr_GSE23066_gene_symbol.rds"
)

write.csv(
  expr_gene,
  "results/expression/expr_GSE23066_gene_symbol.csv"
)

# ============================
# 7. Log
# ============================

sink("results/logs/04_GSE23066_annotation.log")

cat("GSE23066 annotation completed\n")
cat("Date:", Sys.time(), "\n")
cat("RMA dimensions:", dim(expr_rma), "\n")
cat("Gene matrix dimensions:", dim(expr_gene), "\n")

sink()

cat("Annotation completed successfully.\n")

