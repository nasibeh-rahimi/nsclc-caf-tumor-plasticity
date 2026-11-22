

rm(list = ls())

source("scripts/00_project_setup.R")

###Load GSE33532 objects
expr_GSE33532_RMA <- readRDS(
  file.path(processed_dir,
            "expr_GSE33532_RMA.rds")
)

expr_GSE33532_gene_symbol <- readRDS(
  file.path(processed_dir,
            "expr_GSE33532_gene_symbol.rds")
)

expr_GSE33532_gene_level_partial <- readRDS(
  file.path(processed_dir,
            "expr_GSE33532_gene_level_partial.rds")
)

###Check object classes
class(expr_GSE33532_RMA)

class(expr_GSE33532_gene_symbol)

class(expr_GSE33532_gene_level_partial)

###Check dimensions
dim(expr_GSE33532_RMA)

dim(expr_GSE33532_gene_symbol)

dim(expr_GSE33532_gene_level_partial)


###Check sample names
head(colnames(expr_GSE33532_gene_level_partial))

tail(colnames(expr_GSE33532_gene_level_partial))

##Check gene identifiers
head(rownames(expr_GSE33532_gene_level_partial))

###Check duplicated genes
sum(duplicated(
  rownames(expr_GSE33532_gene_level_partial)
))

###Check missing values
sum(is.na(expr_GSE33532_gene_level_partial))

###Check expression distribution
summary(
  as.numeric(expr_GSE33532_gene_level_partial[1:100,1:5])
)

#########
genes <- rownames(expr_GSE33532_gene_symbol)

head(
  genes[grepl("\\.", genes)]
)

sum(grepl("\\.[0-9]+$", genes))


###save a clean final object
expr_GSE33532_gene_final <- expr_GSE33532_gene_symbol

saveRDS(
  expr_GSE33532_gene_final,
  file.path(
    processed_dir,
    "expr_GSE33532_gene_final.rds"
  )
)





