##########08_GSE149955_RMA
########
project_dir <- "D:/NSCLC_Plasticity2"

data_dir      <- file.path(project_dir, "data")
script_dir    <- file.path(project_dir, "scripts")
table_dir     <- file.path(project_dir, "tables")
result_dir    <- file.path(project_dir, "results")

cel_dir       <- file.path(data_dir, "CEL_GSE149955")
processed_dir <- file.path(data_dir, "processed")
qc_dir        <- file.path(result_dir, "figures", "QC")

metadata_file <- file.path(
  table_dir,
  "GSE149955_sample_info.csv"
)

sample_info <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE
)

colnames(sample_info)

head(sample_info)

#########
###Read metadata 
library(stringr)
library(affy)
library(Biobase)

## Locate CEL files
cel_files <- list.files(
  path = cel_dir,
  pattern = "\\.CEL(\\.gz)?$",
  full.names = TRUE
)

cat("CEL files found:", length(cel_files), "\n")

stopifnot(length(cel_files) == 50)
## Read CEL files
raw_data <- ReadAffy(
  filenames = cel_files
)

print(raw_data)
#################################%%%%%%%%%%%
gsm_from_files <- str_extract(
  basename(sampleNames(raw_data)),
  "GSM\\d+"
)

stopifnot(!any(is.na(gsm_from_files)))

match_index <- match(
  sample_info$sample,
  gsm_from_files
)

stopifnot(!any(is.na(match_index)))
stopifnot(length(unique(match_index)) == nrow(sample_info))

raw_data <- raw_data[, match_index]

gsm_after_reordering <- str_extract(
  basename(sampleNames(raw_data)),
  "GSM\\d+"
)

stopifnot(
  all(gsm_after_reordering == sample_info$sample)
)


#################################%%%%%%%%%%%%%%
## Extract GSM IDs
cel_names <- basename(cel_files)
gsm_from_files <- str_extract(cel_names, "GSM\\d+")

## Compare with metadata
head(gsm_from_files)
head(sample_info$sample)

####
cel_files <- list.files(
  path = cel_dir,
  pattern = "\\.CEL(\\.gz)?$",
  full.names = TRUE
)

length(cel_files)

head(basename(cel_files))

###extract the GSM IDs (Validate GSM mapping)
library(stringr)

gsm_from_files <- str_extract(
  basename(cel_files),
  "GSM\\d+"
)

head(gsm_from_files)

###Read the CEL files
raw_data <- ReadAffy(
  filenames = cel_files
)
print(raw_data)
#########
figure_dir <- file.path(project_dir, "figures")

qc_dir <- file.path(
  figure_dir,
  "QC"
)
####
dataset_id <- "GSE149955"


###Raw Quality Control
##==============================================================
## Raw QC - Boxplot


pdf(
  file.path(
    qc_dir,
    paste0(dataset_id, "_Raw_Boxplot.pdf")
  ),
  width = 10,
  height = 6
)

boxplot(
  raw_data,
  col = "lightblue",
  las = 2,
  main = "Raw Affymetrix Intensities",
  ylab = "Raw Intensity"
)

dev.off()

cat("Raw boxplot saved.\n")

##==============================================================
## Raw QC - Density Plot
##==============================================================

pdf(
  file.path(
    qc_dir,
    paste0(dataset_id, "_Raw_Density.pdf")
  ),
  width = 10,
  height = 6
)

hist(
  raw_data,
  lwd = 2,
  main = "Raw Intensity Distributions",
  xlab = "Raw Intensity"
)

dev.off()

cat("Raw density plot saved.\n")

##==============================================================
library(affy)

## RMA Normalization
##==============================================================

cat("Starting RMA normalization...\n")

eset <- rma(raw_data)

cat("RMA normalization completed.\n\n")

print(eset)

##==============================================================
##Create colors
## Sample Colors
##==============================================================

sample_colors <- ifelse(
  sample_info$group == "CAF",
  "salmon",
  "skyblue"
)

##==============================================================
## RMA QC - Boxplot
##==============================================================

pdf(
  file.path(
    qc_dir,
    "GSE149955_RMA_Boxplot.pdf"
  ),
  width = 12,
  height = 7
)

boxplot(
  exprs(eset),
  col = sample_colors,
  las = 2,
  outline = FALSE,
  cex.axis = 0.7,
  main = "GSE149955 - RMA Normalized Expression",
  ylab = "Log2 Expression"
)

legend(
  "topright",
  legend = c("Normal Fibroblast", "Cancer-Associated Fibroblast"),
  fill = c("skyblue", "salmon"),
  bty = "n"
)

dev.off()

cat("Normalized boxplot saved.\n")

##==============================================================
## RMA QC - Density Plot
##==============================================================

pdf(
  file.path(
    qc_dir,
    "GSE149955_RMA_Density.pdf"
  ),
  width = 10,
  height = 7
)

plot(
  density(exprs(eset)[,1]),
  col = sample_colors[1],
  lwd = 2,
  main = "GSE149955 - RMA Density Plot",
  xlab = "Log2 Expression",
  ylim = c(0, 0.6)
)

for(i in 2:ncol(exprs(eset))){
  lines(
    density(exprs(eset)[,i]),
    col = sample_colors[i]
  )
}

legend(
  "topright",
  legend = c("NF", "CAF"),
  col = c("skyblue", "salmon"),
  lwd = 2,
  bty = "n"
)

dev.off()

cat("Normalized density plot saved.\n")

##==============================================================
## PCA of RMA-normalized expression
##==============================================================

# PCA expects samples in rows and variables in columns
pca <- prcomp(
  t(exprs(eset)),
  center = TRUE,
  scale. = FALSE
)

# Percentage of variance explained
percent_var <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
##plot
pdf(
  file.path(
    qc_dir,
    "GSE149955_PCA.pdf"
  ),
  width = 8,
  height = 7
)

plot(
  pca$x[,1],
  pca$x[,2],
  col = sample_colors,
  pch = 19,
  cex = 1.3,
  xlab = paste0("PC1 (", round(percent_var[1],1), "%)"),
  ylab = paste0("PC2 (", round(percent_var[2],1), "%)"),
  main = "GSE149955 - PCA"
)

legend(
  "topright",
  legend = c("NF","CAF"),
  col = c("skyblue","salmon"),
  pch = 19,
  bty = "n"
)

dev.off()

cat("PCA plot saved.\n")

##==============================================================
## Hierarchical Clustering
##==============================================================

sample_dist <- dist(
  t(exprs(eset)),
  method = "euclidean"
)

hc <- hclust(
  sample_dist,
  method = "complete"
)

pdf(
  file.path(
    qc_dir,
    "GSE149955_Hierarchical_Clustering.pdf"
  ),
  width = 12,
  height = 7
)

plot(
  hc,
  labels = sample_info$sample,
  main = "GSE149955 - Hierarchical Clustering",
  xlab = "",
  sub = "",
  cex = 0.7,
  hang = -1
)

dev.off()

cat("Hierarchical clustering saved.\n")

######cleaning the ExpressionSet.
###Step 1 — Make sample names consistent
sampleNames(eset)
###
##==============================================================
## Update sample names
##==============================================================

library(stringr)

sampleNames(eset) <- str_extract(
  sampleNames(eset),
  "GSM\\d+"
)

sampleNames(eset)
###Step 2 — Verify the order
all(sampleNames(eset) == sample_info$sample)
##Step 3 — Attach metadata
##==============================================================
## Attach sample metadata
##==============================================================

library(Biobase)
pd <- AnnotatedDataFrame(data = sample_info)

phenoData(eset) <- pd

###Step 4 — Check it
varLabels(eset)
pData(eset)[1:5, ]

##==============================================================
## Save normalized ExpressionSet
##==============================================================

saveRDS(
  eset,
  file = file.path(
    processed_dir,
    "GSE149955_RMA_ExpressionSet.rds"
  )
)

cat("Normalized ExpressionSet saved successfully.\n")

##==============================================================
## Script completed successfully
##
## Output:
## data/processed/GSE149955_RMA_ExpressionSet.rds
##
## Next script:
## 09_GSE149955_annotation.R
##==============================================================

