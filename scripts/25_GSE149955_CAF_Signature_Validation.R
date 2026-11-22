# ============================================================
# Script 25: Validate CAF-associated transcriptional program
# Discovery cohort: GSE149955
#
# Purpose:
# 1. Reconstruct the CAF-associated score in the original
#    GSE149955 CAF vs NF discovery cohort.
# 2. Confirm that the score separates CAF from NF samples.
# 3. Verify transfer of the fixed 360-gene program to GSE33532.
#
# Important:
# This validates a CAF-associated transcriptional program.
# It does NOT establish direct CAF abundance in bulk tumor.
# ============================================================


# ------------------------------------------------------------
# 1. Project setup
# ------------------------------------------------------------

rm(list = ls())

source("scripts/00_project_setup.R")


library(limma)
library(ggplot2)


# ------------------------------------------------------------
# 2. Define paths
# ------------------------------------------------------------

project_dir <- getwd()

results_dir <- file.path(project_dir, "results")
analysis_dir <- file.path(results_dir, "analysis")
figure_dir <- file.path(project_dir, "figures")

dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 3. Load original GSE149955 expression data
# ------------------------------------------------------------

expr_caf <- readRDS(
  file.path(
    project_dir,
    "data",
    "processed",
    "GSE149955_gene_expression.rds"
  )
)

cat("\nGSE149955 expression dimensions:\n")
print(dim(expr_caf))


# ------------------------------------------------------------
# ------------------------------------------------------------
# 4. Load GSE149955 metadata
# ------------------------------------------------------------

metadata_file <- file.path(
  project_dir,
  "tables",
  "GSE149955_sample_info.csv"
)

stopifnot(file.exists(metadata_file))

metadata_caf <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE
)

cat("\nGSE149955 metadata dimensions:\n")
print(dim(metadata_caf))

cat("\nMetadata columns:\n")
print(colnames(metadata_caf))


# ------------------------------------------------------------
# 5. Check expression / metadata alignment
# ------------------------------------------------------------

stopifnot(
  all(colnames(expr_caf) %in% metadata_caf$sample)
)

metadata_caf <- metadata_caf[
  match(colnames(expr_caf), metadata_caf$sample),
  ,
  drop = FALSE
]

stopifnot(
  all(metadata_caf$sample == colnames(expr_caf))
)

cat("\nExpression and metadata successfully aligned.\n")

  
# ------------------------------------------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------
# 4. Load GSE149955 metadata
# ------------------------------------------------------------

metadata_file <- file.path(
  project_dir,
  "tables",
  "GSE149955_sample_info.csv"
)

stopifnot(file.exists(metadata_file))

metadata_caf <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE
)

cat("\nGSE149955 metadata dimensions:\n")
print(dim(metadata_caf))

cat("\nMetadata columns:\n")
print(colnames(metadata_caf))

# ------------------------------------------------------------
# 5. Check expression / metadata alignment
# ------------------------------------------------------------

stopifnot(
  all(colnames(expr_caf) %in% metadata_caf$sample)
)

metadata_caf <- metadata_caf[
  match(colnames(expr_caf), metadata_caf$sample),
  ,
  drop = FALSE
]

stopifnot(
  all(metadata_caf$sample == colnames(expr_caf))
)

cat("\nExpression and metadata successfully aligned.\n")

# ------------------------------------------------------------
# 6. Identify CAF and NF samples
# ------------------------------------------------------------

cat("\nGroup distribution:\n")
print(table(metadata_caf$group))

stopifnot(
  all(c("CAF", "NF") %in% unique(metadata_caf$group))
)

group <- factor(
  metadata_caf$group,
  levels = c("NF", "CAF")
)


# ------------------------------------------------------------
# 7. Load the original CAF differential-expression results
# ------------------------------------------------------------

caf_up <- readRDS(
  file.path(
    results_dir,
    "expression",
    "GSE149955_DE_upregulated.rds"
  )
)

caf_down <- readRDS(
  file.path(
    results_dir,
    "expression",
    "GSE149955_DE_downregulated.rds"
  )
)

caf_up_genes <- rownames(caf_up)
caf_down_genes <- rownames(caf_down)

cat("\nCAF-associated program:\n")
cat("Upregulated genes:", length(caf_up_genes), "\n")
cat("Downregulated genes:", length(caf_down_genes), "\n")
cat("Total genes:", length(unique(c(caf_up_genes, caf_down_genes))), "\n")


# ------------------------------------------------------------
# 8. Check that the expected 360 genes are available
# ------------------------------------------------------------

up_present <- intersect(caf_up_genes, rownames(expr_caf))
down_present <- intersect(caf_down_genes, rownames(expr_caf))

cat("\nGene coverage in GSE149955:\n")
cat(
  "Upregulated:",
  length(up_present), "/",
  length(caf_up_genes), "\n"
)

cat(
  "Downregulated:",
  length(down_present), "/",
  length(caf_down_genes), "\n"
)

stopifnot(
  length(up_present) == length(caf_up_genes),
  length(down_present) == length(caf_down_genes)
)


# ------------------------------------------------------------
# 9. Construct CAF-associated score
#
# Same scoring principle used for GSE33532:
#
#     mean(upregulated gene z-scores)
#     -
#     mean(downregulated gene z-scores)
#
# Then standardize each component before subtraction.
# ------------------------------------------------------------

up_expr <- expr_caf[up_present, , drop = FALSE]
down_expr <- expr_caf[down_present, , drop = FALSE]


# Gene-wise z-score
up_z <- t(scale(t(up_expr)))
down_z <- t(scale(t(down_expr)))


# Mean component scores
up_component <- colMeans(up_z, na.rm = TRUE)
down_component <- colMeans(down_z, na.rm = TRUE)


# Same balanced scoring logic used for GSE33532
CAF_score_discovery <-
  as.numeric(scale(up_component)) -
  as.numeric(scale(down_component))

names(CAF_score_discovery) <- colnames(expr_caf)


# ------------------------------------------------------------
# 10. Basic score summary
# ------------------------------------------------------------

score_table <- data.frame(
  sample = colnames(expr_caf),
  group = group,
  CAF_associated_score = CAF_score_discovery,
  stringsAsFactors = FALSE
)

cat("\nCAF score summary:\n")
print(summary(score_table$CAF_associated_score))


# ------------------------------------------------------------
# 11. CAF vs NF comparison
# ------------------------------------------------------------

caf_scores <- score_table$CAF_associated_score[
  score_table$group == "CAF"
]

nf_scores <- score_table$CAF_associated_score[
  score_table$group == "NF"
]


# Wilcoxon rank-sum test
wilcox_result <- wilcox.test(
  caf_scores,
  nf_scores,
  exact = FALSE
)


# Mean difference
mean_difference <-
  mean(caf_scores) - mean(nf_scores)


# Median difference
median_difference <-
  median(caf_scores) - median(nf_scores)


# Cohen's d
pooled_sd <- sqrt(
  (
    (length(caf_scores) - 1) * var(caf_scores) +
      (length(nf_scores) - 1) * var(nf_scores)
  ) /
    (
      length(caf_scores) +
        length(nf_scores) - 2
    )
)

cohens_d <- mean_difference / pooled_sd


cat("\n============================================\n")
cat("CAF SIGNATURE VALIDATION\n")
cat("============================================\n")

cat("\nCAF samples:", length(caf_scores), "\n")
cat("NF samples:", length(nf_scores), "\n")

cat("\nMean CAF score:\n")
print(mean(caf_scores))

cat("\nMean NF score:\n")
print(mean(nf_scores))

cat("\nMean difference (CAF - NF):\n")
print(mean_difference)

cat("\nMedian difference (CAF - NF):\n")
print(median_difference)

cat("\nCohen's d:\n")
print(cohens_d)

cat("\nWilcoxon test:\n")
print(wilcox_result)


# ------------------------------------------------------------
# 12. Save statistical results
# ------------------------------------------------------------

validation_results <- data.frame(
  comparison = "CAF_vs_NF",
  n_CAF = length(caf_scores),
  n_NF = length(nf_scores),
  mean_CAF = mean(caf_scores),
  mean_NF = mean(nf_scores),
  mean_difference = mean_difference,
  median_CAF = median(caf_scores),
  median_NF = median(nf_scores),
  median_difference = median_difference,
  cohens_d = cohens_d,
  wilcoxon_W = as.numeric(wilcox_result$statistic),
  wilcoxon_p = wilcox_result$p.value
)

write.csv(
  validation_results,
  file.path(
    analysis_dir,
    "GSE149955_CAF_signature_validation_results.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Save sample-level score table
# ------------------------------------------------------------

write.csv(
  score_table,
  file.path(
    analysis_dir,
    "GSE149955_CAF_signature_discovery_scores.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. Plot CAF vs NF
# ------------------------------------------------------------

p <- ggplot(
  score_table,
  aes(
    x = group,
    y = CAF_associated_score
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.10,
    height = 0,
    size = 2
  ) +
  labs(
    title = "CAF-associated transcriptional program in GSE149955",
    subtitle = "Original CAF vs normal fibroblast discovery cohort",
    x = NULL,
    y = "CAF-associated score"
  ) +
  theme_classic(base_size = 12)

ggsave(
  filename = file.path(
    figure_dir,
    "GSE149955_CAF_signature_CAF_vs_NF.png"
  ),
  plot = p,
  width = 7,
  height = 5,
  dpi = 300
)


# ------------------------------------------------------------
# 15. Verify transfer to GSE33532
# ------------------------------------------------------------

expr_tumor <- readRDS(
  file.path(
    project_dir,
    "data",
    "processed",
    "expr_GSE33532_gene_final.rds"
  )
)

caf_up_gse33532 <- intersect(
  caf_up_genes,
  rownames(expr_tumor)
)

caf_down_gse33532 <- intersect(
  caf_down_genes,
  rownames(expr_tumor)
)

coverage_up <-
  length(caf_up_gse33532) / length(caf_up_genes)

coverage_down <-
  length(caf_down_gse33532) / length(caf_down_genes)

coverage_total <-
  length(
    unique(
      c(caf_up_gse33532, caf_down_gse33532)
    )
  ) /
  length(
    unique(
      c(caf_up_genes, caf_down_genes)
    )
  )


cat("\n============================================\n")
cat("TRANSFER CHECK: GSE33532\n")
cat("============================================\n")

cat(
  "\nUpregulated genes:",
  length(caf_up_gse33532),
  "/",
  length(caf_up_genes),
  "\n"
)

cat(
  "Downregulated genes:",
  length(caf_down_gse33532),
  "/",
  length(caf_down_genes),
  "\n"
)

cat(
  "Total program coverage:",
  round(100 * coverage_total, 2),
  "%\n"
)


# ------------------------------------------------------------
# 16. Save transfer validation
# ------------------------------------------------------------

transfer_results <- data.frame(
  dataset = "GSE33532",
  up_genes_total = length(caf_up_genes),
  up_genes_present = length(caf_up_gse33532),
  up_coverage = coverage_up,
  down_genes_total = length(caf_down_genes),
  down_genes_present = length(caf_down_gse33532),
  down_coverage = coverage_down,
  total_genes_total =
    length(unique(c(caf_up_genes, caf_down_genes))),
  total_genes_present =
    length(unique(c(caf_up_gse33532, caf_down_gse33532))),
  total_coverage = coverage_total
)

write.csv(
  transfer_results,
  file.path(
    analysis_dir,
    "GSE33532_CAF_signature_transfer_validation.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. Write concise validation summary
# ------------------------------------------------------------

summary_file <- file.path(
  analysis_dir,
  "GSE149955_CAF_signature_validation_summary.txt"
)

sink(summary_file)

cat("CAF-ASSOCIATED TRANSCRIPTIONAL PROGRAM VALIDATION\n")
cat("=================================================\n\n")

cat("Discovery dataset: GSE149955\n")
cat("CAF samples:", length(caf_scores), "\n")
cat("NF samples:", length(nf_scores), "\n\n")

cat("CAF-associated program:\n")
cat("Upregulated genes:", length(caf_up_genes), "\n")
cat("Downregulated genes:", length(caf_down_genes), "\n")
cat(
  "Total genes:",
  length(unique(c(caf_up_genes, caf_down_genes))),
  "\n\n"
)

cat("Discovery cohort comparison:\n")
cat("Mean CAF score:", mean(caf_scores), "\n")
cat("Mean NF score:", mean(nf_scores), "\n")
cat("Mean difference:", mean_difference, "\n")
cat("Median difference:", median_difference, "\n")
cat("Cohen's d:", cohens_d, "\n")
cat("Wilcoxon p-value:", wilcox_result$p.value, "\n\n")

cat("Transfer to GSE33532:\n")
cat(
  "Upregulated gene coverage:",
  round(100 * coverage_up, 2),
  "%\n"
)
cat(
  "Downregulated gene coverage:",
  round(100 * coverage_down, 2),
  "%\n"
)
cat(
  "Total program coverage:",
  round(100 * coverage_total, 2),
  "%\n\n"
)

cat("Interpretation:\n")
cat(
  "The CAF-associated transcriptional program was derived from\n",
  "the independent GSE149955 CAF-vs-NF comparison. This analysis\n",
  "checks whether the resulting score behaves consistently with\n",
  "the original discovery labels and whether the fixed gene program\n",
  "can be transferred to GSE33532.\n\n"
)

cat(
  "The score should be interpreted as a CAF-associated transcriptional\n",
  "program score rather than a direct measurement of CAF abundance\n",
  "in bulk tumor samples.\n"
)

sink()


# ------------------------------------------------------------
# 18. Final message
# ------------------------------------------------------------

cat("\n============================================\n")
cat("SCRIPT 25 COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

cat(
  "\nResults saved to:\n",
  analysis_dir,
  "\n"
)

cat(
  "\nFigure saved to:\n",
  file.path(
    figure_dir,
    "GSE149955_CAF_signature_CAF_vs_NF.png"
  ),
  "\n"
)

cat("\nDone.\n")

