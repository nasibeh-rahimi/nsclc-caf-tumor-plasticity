###########################################################
######Do tumor and normal samples separate globally? (PCA)
######Are there any outlier samples? (Hierarchical clustering)
######What is the dominant source of variation? (PCA colored by group and patient)
######Is there obvious batch or patient structure that will influence downstream analyses?
###############################################################################
# Script Name : 14_GSE33532_exploratory_QC.R
#
# Project     : Integrative Transcriptomic Analysis of Cancer-Associated
#               Fibroblasts and Tumor Plasticity in NSCLC
#
# Purpose     :
# Perform post-processing exploratory quality assessment of the
# GSE33532 gene expression dataset.
#
# The objectives are to:
#   1. Assess the global transcriptional structure of the dataset.
#   2. Evaluate separation between tumor and normal lung tissues.
#   3. Detect potential sample outliers.
#   4. Assess patient-specific clustering.
#   5. Visualize transcriptional heterogeneity among NSCLC tumors.
#
# Input Files :
#   processed/expr_GSE33532_gene_final.rds
#   tables/GSE33532_sample_metadata.csv
#
# Output Figures :
#   figures/QC/GSE33532_PCA_Tumor_vs_Normal.pdf
#   figures/QC/GSE33532_Hierarchical_Clustering.pdf
#   figures/QC/GSE33532_PCA_By_Patient.pdf
#   figures/QC/GSE33532_TumorOnly_PCA.pdf
#   figures/QC/GSE33532_SampleCorrelation.pdf
#
# Output Logs :
#   results/logs/14_GSE33532_Exploratory_QC_sessionInfo.txt
#
# Author : Nasibeh Rahimi
# Date   : August 2026
###############################################################################
###Load required libraries
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)

###Create QC directory
qc_dir <- file.path(figure_dir, "QC")

if (!dir.exists(qc_dir)) {
  dir.create(qc_dir, recursive = TRUE)
}

###Load objects
expr_GSE33532_gene_final <- readRDS(
  file.path(
    processed_dir,
    "expr_GSE33532_gene_final.rds"
  )
)

metadata_GSE33532 <- read.csv(
  file.path(
    table_dir,
    "GSE33532_sample_metadata.csv"
  ),
  stringsAsFactors = FALSE
)

###Calculate gene variance
gene_variance <- apply(
  expr_GSE33532_gene_final,
  1,
  var
)

###Select top variable genes
top_variable_genes <- names(
  sort(
    gene_variance,
    decreasing = TRUE
  )
)[1:5000]

expr_pca <- expr_GSE33532_gene_final[
  top_variable_genes,
]

###Perform PCA
pca <- prcomp(
  t(expr_pca),
  center = TRUE,
  scale. = TRUE
)

###Create PCA dataframe
pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2]
)

pca_df <- merge(
  pca_df,
  metadata_GSE33532,
  by.x = "Sample",
  by.y = "sample"
)

###Variance explained
percent_var <- round(
  100 *
    summary(pca)$importance[2,1:2],
  1
)

percent_var

###PCA plot
p1 <- ggplot(
  pca_df,
  aes(
    PC1,
    PC2,
    color = group
  )
) +
  geom_point(
    size = 3,
    alpha = 0.9
  ) +
  labs(
    title = "PCA of GSE33532 Gene Expression",
    subtitle = "Top 5,000 Most Variable Genes",
    x = paste0(
      "PC1 (",
      percent_var[1],
      "%)"
    ),
    y = paste0(
      "PC2 (",
      percent_var[2],
      "%)"
    ),
    color = "Group"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

###Save figure
ggsave(
  filename = file.path(
    qc_dir,
    "GSE33532_PCA_Tumor_vs_Normal.pdf"
  ),
  plot = p1,
  width = 8,
  height = 6
)

###Check group distribution in PCA coordinates
aggregate(
  cbind(PC1, PC2) ~ group,
  data = pca_df,
  mean
)

###Check whether patient identity drives PCA
head(pca_df)



###############################################################################
# Hierarchical clustering
###############################################################################

# Use the same top variable genes selected for PCA
expr_hclust <- expr_pca


# Calculate sample-to-sample distance
sample_distance <- dist(
  t(expr_hclust),
  method = "euclidean"
)


# Hierarchical clustering
hc <- hclust(
  sample_distance,
  method = "complete"
)


# Create annotation for samples
annotation_col <- data.frame(
  Group = metadata_GSE33532$group,
  Patient = metadata_GSE33532$patient_id,
  Tissue = metadata_GSE33532$tissue_code
)


rownames(annotation_col) <- metadata_GSE33532$sample


# Confirm annotation matches expression matrix
all(
  rownames(annotation_col) ==
    colnames(expr_hclust)
)



# Save hierarchical clustering heatmap
###############################################################################

pdf(
  file.path(
    qc_dir,
    "GSE33532_Hierarchical_Clustering.pdf"
  ),
  width = 10,
  height = 8
)

heatmap_result <- pheatmap(
  as.matrix(sample_distance),
  annotation_col = annotation_col,
  cluster_rows = hc,
  cluster_cols = hc,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "GSE33532 Sample Hierarchical Clustering"
)

dev.off()


###############################################################################
# Tumor-only PCA
###############################################################################

# Select tumor samples
tumor_metadata <- metadata_GSE33532[
  metadata_GSE33532$group == "Tumor",
]

# Expression matrix (tumors only)
expr_tumor <- expr_GSE33532_gene_final[
  ,
  tumor_metadata$sample
]

# Check dimensions
dim(expr_tumor)

###Confirm sample order
all(
  colnames(expr_tumor) ==
    tumor_metadata$sample
)

###Select highly variable genes
gene_var_tumor <- apply(
  expr_tumor,
  1,
  var
)

top5000_tumor <- names(
  sort(
    gene_var_tumor,
    decreasing = TRUE
  )[1:5000]
)

expr_tumor_pca <- expr_tumor[
  top5000_tumor,
]

###Run PCA
pca_tumor <- prcomp(
  t(expr_tumor_pca),
  center = TRUE,
  scale. = TRUE
)

###Variance explained
percent_var_tumor <-
  round(
    100 *
      pca_tumor$sdev^2 /
      sum(pca_tumor$sdev^2),
    1
  )

percent_var_tumor[1:10]

###Build plotting dataframe
tumor_pca_df <- data.frame(
  Sample = tumor_metadata$sample,
  PC1 = pca_tumor$x[,1],
  PC2 = pca_tumor$x[,2],
  Patient = tumor_metadata$patient_id,
  Tissue = tumor_metadata$tissue_code
)

###plot
library(ggplot2)

p <- ggplot(
  tumor_pca_df,
  aes(
    PC1,
    PC2,
    color = Patient,
    shape = Tissue
  )
) +
  geom_point(size = 3) +
  labs(
    title = "Tumor-only PCA (GSE33532)",
    x = paste0(
      "PC1 (",
      percent_var_tumor[1],
      "%)"
    ),
    y = paste0(
      "PC2 (",
      percent_var_tumor[2],
      "%)"
    )
  ) +
  theme_bw()

ggsave(
  filename = file.path(
    qc_dir,
    "GSE33532_TumorOnly_PCA.pdf"
  ),
  plot = p,
  width = 8,
  height = 6
)

print(p)
#########################################
##How do you know patient identity is not driving the observed heterogeneity?
###Patient clustering in PCA
library(ggplot2)

tumor_pca_df$Patient <- factor(tumor_pca_df$Patient)

p_patient <- ggplot(
  tumor_pca_df,
  aes(
    PC1,
    PC2,
    color = Patient
  )
) +
  geom_point(size = 3) +
  theme_bw() +
  labs(
    title = "Tumor-only PCA colored by patient",
    x = paste0("PC1 (", percent_var_tumor[1], "%)"),
    y = paste0("PC2 (", percent_var_tumor[2], "%)")
  )
ggsave(
  filename = file.path(
    qc_dir,
    "GSE33532_TumorOnly_PCA_byPatient.pdf"
  ),
  plot = p_patient,
  width = 8,
  height = 6
)

print(p_patient)






























