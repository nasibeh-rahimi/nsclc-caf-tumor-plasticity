############################################################
# 17_GSE33532_CAF_Signature_Transfer.R
#
# Purpose:
# Transfer the CAF-associated transcriptional program
# identified in GSE149955 to the independent GSE33532
# tumor cohort and construct a signed CAF score.
#
# Input:
#   - GSE149955 CAF-upregulated genes
#   - GSE149955 CAF-downregulated genes
#   - GSE33532 final gene-level expression matrix
#   - GSE33532 sample metadata
#
# Output:
#   - CAF score for 80 GSE33532 tumor samples
#   - CAF score diagnostics
#   - CAF gene coverage table
############################################################

############################
# 1. Setup
############################

project_dir <- getwd()

processed_dir <- file.path(
  project_dir,
  "data",
  "processed"
)

expression_result_dir <- file.path(
  project_dir,
  "results",
  "expression"
)

table_dir <- file.path(
  project_dir,
  "tables"
)

caf_result_dir <- file.path(
  project_dir,
  "results",
  "CAF"
)

dir.create(
  caf_result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################
# 2. Load CAF differential-expression results
############################

caf_up_saved <- readRDS(
  file.path(
    expression_result_dir,
    "GSE149955_DE_upregulated.rds"
  )
)

caf_down_saved <- readRDS(
  file.path(
    expression_result_dir,
    "GSE149955_DE_downregulated.rds"
  )
)


# CAF-associated gene sets
caf_up_genes <- rownames(caf_up_saved)
caf_down_genes <- rownames(caf_down_saved)


cat(
  "Number of CAF-upregulated genes:",
  length(caf_up_genes),
  "\n"
)

cat(
  "Number of CAF-downregulated genes:",
  length(caf_down_genes),
  "\n"
)


############################
# 3. Load GSE33532 expression matrix
############################

expr_GSE33532 <- readRDS(
  file.path(
    processed_dir,
    "expr_GSE33532_gene_final.rds"
  )
)

cat(
  "GSE33532 expression dimensions:",
  nrow(expr_GSE33532),
  "genes x",
  ncol(expr_GSE33532),
  "samples\n"
)


############################
# 4. Load GSE33532 metadata
############################

metadata_GSE33532 <- read.csv(
  file.path(
    table_dir,
    "GSE33532_sample_metadata.csv"
  ),
  stringsAsFactors = FALSE
)


############################
# 5. Check metadata-expression alignment
############################

stopifnot(
  all(metadata_GSE33532$sample %in% colnames(expr_GSE33532))
)

stopifnot(
  nrow(metadata_GSE33532) == ncol(expr_GSE33532)
)

metadata_GSE33532 <- metadata_GSE33532[
  match(
    colnames(expr_GSE33532),
    metadata_GSE33532$sample
  ),
  ,
  drop = FALSE
]

stopifnot(
  all(
    metadata_GSE33532$sample ==
      colnames(expr_GSE33532)
  )
)


############################
# 6. Select tumor samples
############################

tumor_samples <- metadata_GSE33532$sample[
  metadata_GSE33532$group == "Tumor"
]

expr_tumor <- expr_GSE33532[
  ,
  tumor_samples,
  drop = FALSE
]

tumor_meta <- metadata_GSE33532[
  metadata_GSE33532$group == "Tumor",
  ,
  drop = FALSE
]

stopifnot(
  ncol(expr_tumor) == 80
)

stopifnot(
  all(tumor_meta$sample == colnames(expr_tumor))
)

cat(
  "Number of tumor samples:",
  ncol(expr_tumor),
  "\n"
)


############################
# 7. Check CAF gene coverage
############################

caf_up_present <- intersect(
  caf_up_genes,
  rownames(expr_tumor)
)

caf_down_present <- intersect(
  caf_down_genes,
  rownames(expr_tumor)
)

caf_up_missing <- setdiff(
  caf_up_genes,
  rownames(expr_tumor)
)

caf_down_missing <- setdiff(
  caf_down_genes,
  rownames(expr_tumor)
)


coverage_table <- data.frame(
  direction = c("up", "down"),
  total_genes = c(
    length(caf_up_genes),
    length(caf_down_genes)
  ),
  present_genes = c(
    length(caf_up_present),
    length(caf_down_present)
  ),
  missing_genes = c(
    length(caf_up_missing),
    length(caf_down_missing)
  )
)

coverage_table$coverage_percent <- (
  coverage_table$present_genes /
    coverage_table$total_genes
) * 100


print(coverage_table)


write.csv(
  coverage_table,
  file.path(
    caf_result_dir,
    "GSE33532_CAF_gene_coverage.csv"
  ),
  row.names = FALSE
)


############################
# 8. Extract CAF gene expression
############################

caf_up_expr <- expr_tumor[
  caf_up_present,
  ,
  drop = FALSE
]

caf_down_expr <- expr_tumor[
  caf_down_present,
  ,
  drop = FALSE
]


############################
# 9. Gene-wise standardization
############################

caf_up_z <- t(
  scale(
    t(caf_up_expr)
  )
)

caf_down_z <- t(
  scale(
    t(caf_down_expr)
  )
)


############################
# 10. Calculate CAF components
############################

caf_up_component <- colMeans(
  caf_up_z,
  na.rm = TRUE
)

caf_down_component <- colMeans(
  caf_down_z,
  na.rm = TRUE
)


############################
# 11. Construct signed CAF score
############################

CAF_score <- (
  scale(caf_up_component)[, 1] -
    scale(caf_down_component)[, 1]
)

names(CAF_score) <- colnames(expr_tumor)


############################
# 12. Check score alignment
############################

stopifnot(
  identical(
    names(CAF_score),
    colnames(expr_tumor)
  )
)

cat(
  "CAF score alignment:",
  identical(
    names(CAF_score),
    colnames(expr_tumor)
  ),
  "\n"
)


############################
# 13. CAF score summary
############################

print(
  summary(CAF_score)
)

print(
  sd(CAF_score)
)


############################
# 14. Save CAF score
############################

saveRDS(
  CAF_score,
  file.path(
    caf_result_dir,
    "GSE33532_CAF_score.rds"
  )
)


############################
# 15. Save CAF components
############################

CAF_components <- data.frame(
  sample = names(CAF_score),
  CAF_up_component = as.numeric(
    caf_up_component
  ),
  CAF_down_component = as.numeric(
    caf_down_component
  ),
  CAF_score = as.numeric(
    CAF_score
  )
)

write.csv(
  CAF_components,
  file.path(
    caf_result_dir,
    "GSE33532_CAF_score_components.csv"
  ),
  row.names = FALSE
)


############################
# 16. Save tumor metadata
############################

write.csv(
  tumor_meta,
  file.path(
    caf_result_dir,
    "GSE33532_tumor_metadata.csv"
  ),
  row.names = FALSE
)


############################
# 17. Save session information
############################

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    caf_result_dir,
    "17_CAF_Signature_Transfer_sessionInfo.txt"
  )
)


cat(
  "\nCAF signature transfer completed successfully.\n"
)
