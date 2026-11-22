##############################################################
##
## Script:
## 09_GSE149955_annotation.R
##
## Purpose:
## Annotate Affymetrix HG-U133 Plus 2.0 probe sets,
## collapse multiple probes per gene,
## and generate a gene-level expression matrix.
##
## Input:
## data/processed/GSE149955_RMA_ExpressionSet.rds
##
## Output:
## data/processed/GSE149955_gene_expression.rds
## data/processed/GSE149955_probe_annotation_full.rds
##
## Tables:
## tables/GSE149955_probe_annotation.csv
## tables/GSE149955_annotation_statistics.csv
##
##############################################################

##==============================================================
## Load packages
##==============================================================

library(Biobase)
library(AnnotationDbi)
library(hgu133plus2.db)

##==============================================================
## Load normalized ExpressionSet
##==============================================================

eset <- readRDS(
  file.path(
    processed_dir,
    "GSE149955_RMA_ExpressionSet.rds"
  )
)

cat("ExpressionSet loaded successfully.\n")
print(eset)

##==============================================================
## Retrieve probe annotation
##==============================================================

probe2gene <- mapIds(
  hgu133plus2.db,
  keys      = featureNames(eset),
  column    = "SYMBOL",
  keytype   = "PROBEID",
  multiVals = "first"
)

annotation_df <- data.frame(
  probe = names(probe2gene),
  gene  = unname(probe2gene),
  stringsAsFactors = FALSE
)

##==============================================================
## Remove probes without gene symbols
##==============================================================

annotation_df <- annotation_df[
  !is.na(annotation_df$gene),
]

cat("Annotated probes:", nrow(annotation_df), "\n")

##==============================================================
## Extract expression matrix
##==============================================================

expr_matrix <- exprs(eset)

##==============================================================
## Calculate probe variance
##==============================================================

probe_variance <- apply(
  expr_matrix,
  1,
  var
)

annotation_df$variance <-
  probe_variance[annotation_df$probe]

##==============================================================
## Save complete probe annotation
##==============================================================

saveRDS(
  annotation_df,
  file.path(
    processed_dir,
    "GSE149955_probe2gene_mapping.rds"
  )
)

##==============================================================
## Select highest-variance probe for each gene
##==============================================================

annotation_best <- annotation_df[
  order(
    annotation_df$gene,
    -annotation_df$variance
  ),
]

annotation_best <- annotation_best[
  !duplicated(annotation_best$gene),
]

rownames(annotation_best) <- NULL

##==============================================================
## Build gene-level expression matrix
##==============================================================

expr_gene <- expr_matrix[
  annotation_best$probe,
]

rownames(expr_gene) <- annotation_best$gene

##==============================================================
## Save gene-level expression matrix
##==============================================================

saveRDS(
  expr_gene,
  file.path(
    processed_dir,
    "GSE149955_gene_expression.rds"
  )
)

##==============================================================
## Add probe count information
##==============================================================

probe_count <- table(annotation_df$gene)

annotation_best$n_candidate_probes <-
  as.integer(
    probe_count[annotation_best$gene]
  )

##==============================================================
## Save annotation table
##==============================================================

write.csv(
  annotation_best,
  file.path(
    table_dir,
    "GSE149955_probe_annotation.csv"
  ),
  row.names = FALSE
)

##==============================================================
## Annotation statistics
##==============================================================

annotation_stats <- data.frame(
  
  Metric = c(
    "Original probes",
    "Annotated probes",
    "Unique genes",
    "Genes with one probe",
    "Genes with multiple probes"
  ),
  
  Value = c(
    nrow(expr_matrix),
    nrow(annotation_df),
    length(probe_count),
    sum(probe_count == 1),
    sum(probe_count > 1)
  )
  
)

write.csv(
  annotation_stats,
  file.path(
    table_dir,
    "GSE149955_annotation_statistics.csv"
  ),
  row.names = FALSE
)

##==============================================================
## Save sample metadata
##==============================================================

write.csv(
  pData(eset),
  file.path(
    table_dir,
    "GSE149955_metadata_processed.csv"
  ),
  row.names = FALSE
)

##==============================================================
## Summary
##==============================================================

cat("\n")
cat("=========================================\n")
cat(" Script 09 completed successfully\n")
cat("=========================================\n")
cat("Original probes      :", nrow(expr_matrix), "\n")
cat("Annotated probes     :", nrow(annotation_df), "\n")
cat("Unique genes         :", length(probe_count), "\n")
cat("Expression matrix    :", dim(expr_gene)[1], "genes x",
    dim(expr_gene)[2], "samples\n")
cat("=========================================\n")