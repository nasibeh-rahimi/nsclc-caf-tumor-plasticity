############################################################
################################################################
####05_CAF_differential_expression
source("scripts/00_project_setup.R")

# =====================================
# Load expression matrix
# =====================================

expr_gene <- readRDS(
  "data/processed/expr_GSE23066_gene_symbol.rds"
)

cat("Expression matrix:", dim(expr_gene), "\n")

cat("Expression matrix:\n")
print(dim(expr_gene))


sample_info <- data.frame(
  sample = colnames(expr_gene),
  stringsAsFactors = FALSE
)

sample_info$group <- ifelse(
  grepl("_N", sample_info$sample),
  "Normal",
  "CAF"
)

sample_info$group <- factor(
  sample_info$group,
  levels = c("Normal", "CAF")
)

write.csv(
  sample_info,
  file.path(table_dir, "GSE23066_sample_info.csv"),
  row.names = FALSE
)

sample_info
cat("✓ Sample information created\n")

# =====================================
# PCA
# =====================================

pca <- prcomp(
  t(expr_gene),
  scale. = TRUE
)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Group = sample_info$group,
  Sample = sample_info$sample
)

ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 4) +
  geom_text(aes(label = Group), vjust = -1) +
  theme_classic() +
  labs(
    title = "PCA of GSE23066",
    x = paste0(
      "PC1 (",
      round(summary(pca)$importance[2,1] * 100, 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(summary(pca)$importance[2,2] * 100, 1),
      "%)"
    )
  )

ggsave(
  file.path(figure_dir, "PCA_GSE23066.pdf"),
  width = 6,
  height = 5
)

cat("✓ PCA completed\n")

design <- model.matrix(~0 + group, data = sample_info)
colnames(design) <- c("CAF", "Normal")
print(design)

write.csv(
  design,
  file.path(table_dir, "GSE23066_design_matrix.csv")
)

cat("✓ Design matrix created\n")



###Fit the model

library(limma)

fit <- lmFit(
  expr_gene,
  design
)

###Fit the model
contrast_matrix <- makeContrasts(
  CAF_vs_Normal = CAF - Normal,
  levels = design
)

print(contrast_matrix)


###Apply the contrast
fit2 <- contrasts.fit(
  fit,
  contrast_matrix
)

###Empirical Bayes
fit2 <- eBayes(fit2)

saveRDS(
  fit2,
  file.path(
    expression_dir,
    "GSE23066_limma_fit.rds"
  )
)

###Extract the complete results
de_results <- topTable(
  fit2,
  coef = "CAF_vs_Normal",
  number = Inf,
  sort.by = "P"
)


###Save the complete DE table
write.csv(
  de_results,
  file.path(expression_dir, "GSE23066_DE_results.csv"),
  row.names = TRUE
)

###Save the candidate gene set
caf_candidates <- subset(
  de_results,
  P.Value < 0.01
)

write.csv(
  caf_candidates,
  file.path(table_dir, "CAF_candidate_genes_P001.csv"),
  row.names = TRUE
)

cat("Candidate genes:", nrow(caf_candidates), "\n")
cat(
  "Upregulated:",
  sum(caf_candidates$logFC > 0),
  "\n"
)

cat(
  "Downregulated:",
  sum(caf_candidates$logFC < 0),
  "\n"
)



