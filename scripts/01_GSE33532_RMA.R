############################################################
# 01_RMA_normalization.R
# GSE33532 RMA preprocessing (Affymetrix CEL)
############################################################

source("scripts/00_project_setup.R")

library(affy)

# ============================
# 1. CEL path
# ============================

cel_path <- "data/CEL_GSE33532"

cel_files <- list.files(
  cel_path,
  pattern = "\\.CEL$",
  full.names = TRUE
)

cat("CEL files found:", length(cel_files), "\n")

if(length(cel_files) != 100){
  warning("Expected ~100 CEL files. Check dataset completeness.")
}

# ============================
# 2. Read CEL files
# ============================

raw_affy <- ReadAffy(filenames = cel_files)

cat("CEL files loaded into Affy object\n")

# ============================
# 3. RMA normalization
# ============================

expr_rma <- rma(raw_affy)

# ============================
# 4. Expression matrix
# ============================

expr_matrix <- exprs(expr_rma)

cat("Expression matrix dimensions:\n")
print(dim(expr_matrix))

# ============================
# 5. Quality sanity checks
# ============================

cat("NA count:", sum(is.na(expr_matrix)), "\n")

# ============================
# 6. Save outputs
# ============================

saveRDS(
  expr_matrix,
  "data/processed/expr_GSE33532_RMA.rds"
)

write.csv(
  expr_matrix,
  "results/expression/expr_GSE33532_RMA.csv"
)

cat("RMA normalization completed successfully\n")
