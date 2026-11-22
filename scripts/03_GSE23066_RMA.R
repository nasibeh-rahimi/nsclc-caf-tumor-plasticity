# ============================
# Read CEL files
# ============================

cel_path <- "data/CEL_GSE23066"

cel_files <- list.files(
  cel_path,
  pattern = "\\.CEL$",
  full.names = TRUE
)

cat("Number of CEL files:", length(cel_files), "\n")

# Safety check
stopifnot(length(cel_files) == 10)

# ============================
# Read Affymetrix CEL files
# ============================
library(affy)

raw_affy <- ReadAffy(
  filenames = cel_files
)

# ============================
# RMA normalization
# ============================

eset <- rma(raw_affy)

expr23066 <- exprs(eset)

# ============================
# Summary
# ============================

cat("Expression matrix dimensions:\n")
print(dim(expr23066))

cat("Platform:", annotation(eset), "\n")

# ============================
# Save
# ============================

saveRDS(
  expr23066,
  "data/processed/expr_GSE23066_RMA.rds"
)

write.csv(
  expr23066,
  "results/expression/expr_GSE23066_RMA.csv"
)

# ============================
# Log
dir.create("results/logs", recursive = TRUE)

# ============================

sink("results/logs/03_GSE23066_RMA.log")

cat("GSE23066 RMA completed\n")
cat("Date:", Sys.time(), "\n")
cat("Dimensions:", dim(expr23066), "\n")
cat("Platform:", annotation(eset), "\n")

sink()

cat("Script completed successfully.\n")

dim(expr23066)

annotation(eset)

file.exists("data/processed/expr_GSE23066_RMA.rds")
