#####10_GSE149955_differential_expression

project_dir <- "D:/NSCLC_Plasticity2"

data_dir      <- file.path(project_dir, "data")
processed_dir <- file.path(data_dir, "processed")

table_dir     <- file.path(project_dir, "tables")

figure_dir    <- file.path(project_dir, "figures")
de_figure_dir <- file.path(figure_dir, "DE")

result_dir    <- file.path(project_dir, "results")
###Identify genes that are differentially expressed between matched CAF and NF samples using a paired limma model.

library(Biobase)
library(ggplot2)
library(pheatmap)

##==============================================================
## Create output directories
##==============================================================

de_figure_dir <- file.path(
  figure_dir,
  "DE"
)

if (!dir.exists(de_figure_dir)) {
  dir.create(
    de_figure_dir,
    recursive = TRUE
  )
}

############
###verifying that the metadata and expression matrix are perfectly aligned
expr_gene <- readRDS(
  file.path(
    processed_dir,
    "GSE149955_gene_expression.rds"
  )
)

metadata <- read.csv(
  file.path(
    table_dir,
    "GSE149955_metadata_processed.csv"
  ),
  stringsAsFactors = FALSE
)

identical(
  colnames(expr_gene),
  metadata$sample
)

##==============================================================
## Create factors
##==============================================================

metadata$group <- factor(
  metadata$group,
  levels = c("NF", "CAF")
)

metadata$patient <- factor(metadata$patient)

###build the paired design matrix

packageVersion("limma")

library(limma)

design <- model.matrix(
  ~ patient + group,
  data = metadata
)

colnames(design)


##==============================================================
## Fit paired linear model
##==============================================================

fit <- lmFit(
  expr_gene,
  design
)

fit <- eBayes(fit)

###save the fitted model
saveRDS(
  fit,
  file.path(
    result_dir,
    "expression",
    "GSE149955_limma_fit.rds"
  )
)

##==============================================================
## Differential expression analysis
##==============================================================

deg_results <- topTable(
  fit,
  coef = "groupCAF",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

##==============================================================
## Save complete differential expression results
##==============================================================

write.csv(
  deg_results,
  file.path(
    table_dir,
    "GSE149955_DEGs_all.csv"
  ),
  row.names = TRUE
)

###define significant genes
deg_sig <- deg_results[
  deg_results$adj.P.Val < 0.05,
]

###Save significant genes
write.csv(
  deg_sig,
  file.path(
    table_dir,
    "GSE149955_DEGs_significant.csv"
  ),
  row.names = TRUE
)

##==============================================================
## MA Plot:(quality check of the statistical model)
##==============================================================

pdf(
  file.path(
    de_figure_dir,
    "GSE149955_MAplot.pdf"
  ),
  width = 7,
  height = 6
)

plotMA(
  fit,
  coef = "groupCAF",
  main = "CAF vs Normal Fibroblasts",
  xlab = "Average log2 Expression",
  ylab = "Log2 Fold Change"
)

abline(
  h = 0,
  col = "red",
  lty = 2
)

dev.off()

###Volcano plot (global overview of differential expression)
##==============================================================

library(ggplot2)

## Prepare volcano plot data
volcano_data <- deg_results

volcano_data$Significant <-
  volcano_data$adj.P.Val < 0.05

volcano_data$minusLog10FDR <-
  -log10(volcano_data$adj.P.Val)

## Create volcano plot
volcano_plot <- ggplot(
  volcano_data,
  aes(
    x = logFC,
    y = minusLog10FDR
  )
) +
  geom_point(
    aes(color = Significant),
    alpha = 0.6,
    size = 1.5
  ) +
  scale_color_manual(
    values = c(
      "FALSE" = "grey70",
      "TRUE"  = "firebrick"
    ),
    name = "Significant"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "Differential Expression: CAF vs Normal Fibroblasts",
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~Adjusted~P~Value)
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "right"
  )

## Display in RStudio
print(volcano_plot)

## Save as PDF
ggsave(
  filename = file.path(
    de_figure_dir,
    "GSE149955_Volcano.pdf"
  ),
  plot = volcano_plot,
  width = 8,
  height = 7,
  units = "in"
)
##########
###heat map 
##==============================================================
## Heatmap of Top 50 Differentially Expressed Genes
##==============================================================

library(pheatmap)

## Select top 50 genes
top50 <- rownames(deg_results)[1:50]

## Build expression matrix
heatmap_matrix <- expr_gene[top50, ]

## Scale each gene (Z-score)
heatmap_matrix <- t(
  scale(
    t(heatmap_matrix)
  )
)

## Sample annotation
annotation_col <- data.frame(
  Group = metadata$group,
  row.names = metadata$sample
)

## Draw heatmap in RStudio
pheatmap(
  heatmap_matrix,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  fontsize_row = 7,
  main = "Top 50 Differentially Expressed Genes"
)

## Save publication-quality PDF
pheatmap(
  heatmap_matrix,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  fontsize_row = 7,
  main = "Top 50 Differentially Expressed Genes",
  filename = file.path(
    de_figure_dir,
    "GSE149955_Top50_Heatmap.pdf"
  ),
  width = 9,
  height = 10
)
################################################
##check
summary(deg_results$adj.P.Val)
sum(deg_results$adj.P.Val < 0.05)
sum(abs(deg_results$logFC) > 1)
summary(deg_results$logFC)

up_genes <- sum(
  deg_results$adj.P.Val < 0.05 &
    deg_results$logFC > 0
)

down_genes <- sum(
  deg_results$adj.P.Val < 0.05 &
    deg_results$logFC < 0
)

up_genes
down_genes
#################
##==============================================================
## Significant Differentially Expressed Genes
##==============================================================

deg_sig <- deg_results[
  deg_results$adj.P.Val < 0.05,
]

cat("Significant DEGs:", nrow(deg_sig), "\n")

##==============================================================
## Upregulated and Downregulated Genes
##==============================================================

deg_up <- deg_sig[
  deg_sig$logFC > 0,
]

deg_down <- deg_sig[
  deg_sig$logFC < 0,
]

cat("Upregulated genes:", nrow(deg_up), "\n")
cat("Downregulated genes:", nrow(deg_down), "\n")

##==============================================================
## Save R Objects
##==============================================================

saveRDS(
  deg_results,
  file.path(
    expression_dir,
    "GSE149955_DE_results.rds"
  )
)

saveRDS(
  deg_sig,
  file.path(
    expression_dir,
    "GSE149955_DE_significant.rds"
  )
)

saveRDS(
  deg_up,
  file.path(
    expression_dir,
    "GSE149955_DE_upregulated.rds"
  )
)

saveRDS(
  deg_down,
  file.path(
    expression_dir,
    "GSE149955_DE_downregulated.rds"
  )
)

##==============================================================
## Save Tables
##==============================================================

write.csv(
  deg_results,
  file.path(
    table_dir,
    "GSE149955_DE_results.csv"
  ),
  row.names = TRUE
)

write.csv(
  deg_sig,
  file.path(
    table_dir,
    "GSE149955_DE_significant.csv"
  ),
  row.names = TRUE
)

write.csv(
  deg_up,
  file.path(
    table_dir,
    "GSE149955_DE_upregulated.csv"
  ),
  row.names = TRUE
)

write.csv(
  deg_down,
  file.path(
    table_dir,
    "GSE149955_DE_downregulated.csv"
  ),
  row.names = TRUE
)

##==============================================================
## Summary
##==============================================================

cat("\n")
cat("=========================================\n")
cat(" Differential Expression Summary\n")
cat("=========================================\n")
cat("Total genes tested:      ", nrow(deg_results), "\n")
cat("Significant DEGs:        ", nrow(deg_sig), "\n")
cat("Upregulated genes:       ", nrow(deg_up), "\n")
cat("Downregulated genes:     ", nrow(deg_down), "\n")
cat("=========================================\n")

##==============================================================
## Save Session Information
##==============================================================

writeLines(
  capture.output(sessionInfo()),
  file.path(
        result_dir,
        "logs",
        "10_GSE149955_Differential_Expression_sessionInfo.txt"
    )
)








