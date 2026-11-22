#########################################
###create GSE33532 metadata
sample_names <- colnames(expr_GSE33532_gene_final)

head(sample_names, 20)


metadata_test <- data.frame(
  sample = sample_names
)

head(metadata_test)



metadata_test$patient_id <- sub(
  ".*_(\\d+)_.*",
  "\\1",
  metadata_test$sample
)

head(metadata_test)


metadata_test$tissue_code <- sub(
  ".*_([A-Z])_TK\\.CEL",
  "\\1",
  metadata_test$sample
)

head(metadata_test)
################

metadata_GSE33532 <- metadata_test

metadata_GSE33532$group <- ifelse(
  metadata_GSE33532$tissue_code == "N",
  "Normal",
  "Tumor"
)

metadata_GSE33532$group <- factor(
  metadata_GSE33532$group,
  levels = c("Normal", "Tumor")
)

head(metadata_GSE33532)


###save
write.csv(
  metadata_GSE33532,
  file.path(
    table_dir,
    "GSE33532_sample_metadata.csv"
  ),
  row.names = FALSE
)





