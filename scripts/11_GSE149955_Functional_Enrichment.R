############################################
##### 11_GSE149955_Functional_Enrichment
####### its goal: Take the 360 DEGs and answer one biological question
#### What biological processes distinguish CAFs from normal fibroblasts?

##==============================================================
## Load packages
##==============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(DOSE)
library(ggplot2)
library(AnnotationDbi)

##==============================================================
## Project directories
##==============================================================

project_dir <- getwd()

processed_dir <- file.path(
  project_dir,
  "data",
  "processed"
)

result_dir <- file.path(
  project_dir,
  "results"
)

expression_dir <- file.path(
  result_dir,
  "expression"
)

enrichment_dir <- file.path(
  result_dir,
  "enrichment"
)

table_dir <- file.path(
  project_dir,
  "tables"
)

figure_dir <- file.path(
  project_dir,
  "figures"
)

enrichment_figure_dir <-
  file.path(
    figure_dir,
    "Enrichment"
  )

dir.create(
  enrichment_figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

##==============================================================
## Load differential expression results
##==============================================================

deg_results <- readRDS(
  file.path(
    expression_dir,
    "GSE149955_DE_results.rds"
  )
)

deg_sig <- readRDS(
  file.path(
    expression_dir,
    "GSE149955_DE_significant.rds"
  )
)

deg_up <- readRDS(
  file.path(
    expression_dir,
    "GSE149955_DE_upregulated.rds"
  )
)

deg_down <- readRDS(
  file.path(
    expression_dir,
    "GSE149955_DE_downregulated.rds"
  )
)

##==============================================================
## Load gene-level expression matrix
## Used to define the enrichment background universe
##==============================================================

expr_gene <- readRDS(
  file.path(
    processed_dir,
    "GSE149955_gene_expression.rds"
  )
)

##==============================================================
## Check imported results
##==============================================================

cat(
  "Total genes:",
  nrow(deg_results),
  "\n"
)

cat(
  "Significant genes:",
  nrow(deg_sig),
  "\n"
)

cat(
  "Upregulated genes:",
  nrow(deg_up),
  "\n"
)

cat(
  "Downregulated genes:",
  nrow(deg_down),
  "\n"
)

cat(
  "Genes in expression matrix:",
  nrow(expr_gene),
  "\n"
)

##==============================================================
## Extract gene symbols
##==============================================================

genes_up <- rownames(deg_up)

genes_down <- rownames(deg_down)

genes_all <- rownames(deg_results)

length(genes_up)

length(genes_down)

##==============================================================
## Define enrichment background universe
##
## The universe is all genes represented in the
## GSE149955 gene-level expression matrix.
##
## This is preferable to using the entire human genome
## because only genes measurable on this experiment's
## platform could have entered the DE analysis.
##==============================================================

universe_symbols <- rownames(expr_gene)

universe_entrez_df <- bitr(
  universe_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

universe_entrez <- unique(
  universe_entrez_df$ENTREZID
)

cat(
  "Background gene symbols:",
  length(universe_symbols),
  "\n"
)

cat(
  "Background Entrez IDs:",
  length(universe_entrez),
  "\n"
)

##==============================================================
## Convert DE gene symbols to Entrez IDs
##==============================================================

up_entrez <- bitr(
  genes_up,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

down_entrez <- bitr(
  genes_down,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

## Remove duplicated Entrez IDs
up_entrez <- up_entrez[
  !duplicated(up_entrez$ENTREZID),
]

down_entrez <- down_entrez[
  !duplicated(down_entrez$ENTREZID),
]

##==============================================================
## Check ID conversion
##==============================================================

head(up_entrez)

head(down_entrez)

cat(
  "Upregulated Entrez IDs:",
  nrow(up_entrez),
  "\n"
)

cat(
  "Downregulated Entrez IDs:",
  nrow(down_entrez),
  "\n"
)

##==============================================================
## GO Biological Process Enrichment
## Upregulated genes
##==============================================================

ego_up_BP <- enrichGO(
  gene          = up_entrez$ENTREZID,
  universe      = universe_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

head(ego_up_BP)

##==============================================================
## GO Biological Process Enrichment
## Downregulated genes
##==============================================================

ego_down_BP <- enrichGO(
  gene          = down_entrez$ENTREZID,
  universe      = universe_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

head(ego_down_BP)

cat(
  "GO BP Upregulated terms:",
  nrow(as.data.frame(ego_up_BP)),
  "\n"
)

cat(
  "GO BP Downregulated terms:",
  nrow(as.data.frame(ego_down_BP)),
  "\n"
)

##==============================================================
## Which upregulated genes participated in significant GO BP?
##==============================================================

ego_up_BP_df <- as.data.frame(
  ego_up_BP
)

if (nrow(ego_up_BP_df) > 0) {
  
  ## Because readable = TRUE was used in enrichGO(),
  ## geneID already contains gene SYMBOLs.
  
  genes_in_GO_BP_symbols <- unique(
    unlist(
      strsplit(
        ego_up_BP_df$geneID,
        "/"
      )
    )
  )
  
  genes_in_GO_BP_symbols <- trimws(
    genes_in_GO_BP_symbols
  )
  
  genes_in_GO_BP_symbols <- unique(
    genes_in_GO_BP_symbols[
      !is.na(genes_in_GO_BP_symbols) &
        genes_in_GO_BP_symbols != ""
    ]
  )
  
  ## Upregulated genes that did not appear
  ## in any significant GO BP term
  genes_not_in_GO_BP <- setdiff(
    genes_up,
    genes_in_GO_BP_symbols
  )
  
} else {
  
  ## No significant GO BP terms
  genes_not_in_GO_BP <- genes_up
  
}

##==============================================================
## Save upregulated genes not represented in significant GO BP
##==============================================================

write.csv(
  data.frame(
    Gene = genes_not_in_GO_BP
  ),
  file.path(
    table_dir,
    "GSE149955_Upregulated_Genes_Not_in_GO_BP.csv"
  ),
  row.names = FALSE
)

cat(
  "Upregulated genes not represented in significant GO BP:",
  length(genes_not_in_GO_BP),
  "\n"
)

##==============================================================
## GO Cellular Component (CC)
##==============================================================

##-------------------------
## Run enrichment
##-------------------------

ego_up_CC <- enrichGO(
  gene          = up_entrez$ENTREZID,
  universe      = universe_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_down_CC <- enrichGO(
  gene          = down_entrez$ENTREZID,
  universe      = universe_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

##-------------------------
## Inspect results
##-------------------------

cc_up_df <- as.data.frame(
  ego_up_CC
)

cc_down_df <- as.data.frame(
  ego_down_CC
)

cat(
  "GO CC Upregulated terms:",
  nrow(cc_up_df),
  "\n"
)

cat(
  "GO CC Downregulated terms:",
  nrow(cc_down_df),
  "\n"
)

head(cc_up_df)

head(cc_down_df)

##-------------------------
## Save tables
##-------------------------

if (nrow(cc_up_df) > 0) {
  
  write.csv(
    cc_up_df,
    file.path(
      table_dir,
      "GSE149955_GO_CC_Upregulated.csv"
    ),
    row.names = FALSE
  )
  
} else {
  
  message(
    "No significant GO CC terms for upregulated genes."
  )
  
}

if (nrow(cc_down_df) > 0) {
  
  write.csv(
    cc_down_df,
    file.path(
      table_dir,
      "GSE149955_GO_CC_Downregulated.csv"
    ),
    row.names = FALSE
  )
  
} else {
  
  message(
    "No significant GO CC terms for downregulated genes."
  )
  
}

##-------------------------
## Dotplots
##-------------------------

if (nrow(cc_up_df) > 0) {
  
  cc_up_plot <- dotplot(
    ego_up_CC,
    showCategory = 20
  ) +
    ggtitle(
      "GO Cellular Component: Upregulated Genes"
    )
  
  print(cc_up_plot)
  
  ggsave(
    filename = file.path(
      enrichment_figure_dir,
      "GSE149955_GO_CC_Upregulated_Dotplot.pdf"
    ),
    plot = cc_up_plot,
    width = 10,
    height = 8
  )
  
} else {
  
  message(
    "Skipping GO CC Upregulated dotplot."
  )
  
}

if (nrow(cc_down_df) > 0) {
  
  cc_down_plot <- dotplot(
    ego_down_CC,
    showCategory = 20
  ) +
    ggtitle(
      "GO Cellular Component: Downregulated Genes"
    )
  
  print(cc_down_plot)
  
  ggsave(
    filename = file.path(
      enrichment_figure_dir,
      "GSE149955_GO_CC_Downregulated_Dotplot.pdf"
    ),
    plot = cc_down_plot,
    width = 10,
    height = 8
  )
  
} else {
  
  message(
    "Skipping GO CC Downregulated dotplot."
  )
  
}

##==============================================================
## GO Molecular Function (MF)
##==============================================================

##-------------------------
## Run enrichment
##-------------------------

ego_up_MF <- enrichGO(
  gene          = up_entrez$ENTREZID,
  universe      = universe_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_down_MF <- enrichGO(
  gene          = down_entrez$ENTREZID,
  universe      = universe_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

##-------------------------
## Inspect results
##-------------------------

mf_up_df <- as.data.frame(
  ego_up_MF
)

mf_down_df <- as.data.frame(
  ego_down_MF
)

cat(
  "GO MF Upregulated terms:",
  nrow(mf_up_df),
  "\n"
)

cat(
  "GO MF Downregulated terms:",
  nrow(mf_down_df),
  "\n"
)

head(mf_up_df)

head(mf_down_df)

##==============================================================
## Reactome Pathway Enrichment
##==============================================================

reactome_up <- enrichPathway(
  gene          = up_entrez$ENTREZID,
  universe      = universe_entrez,
  organism      = "human",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

reactome_down <- enrichPathway(
  gene          = down_entrez$ENTREZID,
  universe      = universe_entrez,
  organism      = "human",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

##==============================================================
## Inspect Reactome results
##==============================================================

reactome_up_df <- as.data.frame(
  reactome_up
)

reactome_down_df <- as.data.frame(
  reactome_down
)

cat(
  "Reactome Upregulated:",
  nrow(reactome_up_df),
  "\n"
)

cat(
  "Reactome Downregulated:",
  nrow(reactome_down_df),
  "\n"
)

head(reactome_up_df)

head(reactome_down_df)

##==============================================================
## Save Reactome tables
##==============================================================

if (nrow(reactome_up_df) > 0) {
  
  write.csv(
    reactome_up_df,
    file.path(
      table_dir,
      "GSE149955_Reactome_Upregulated.csv"
    ),
    row.names = FALSE
  )
  
} else {
  
  message(
    "No significant Reactome pathways for upregulated genes."
  )
}

if (nrow(reactome_down_df) > 0) {
  
  write.csv(
    reactome_down_df,
    file.path(
      table_dir,
      "GSE149955_Reactome_Downregulated.csv"
    ),
    row.names = FALSE
  )
  
} else {
  
  message(
    "No significant Reactome pathways for downregulated genes."
  )
}

##==============================================================
## Reactome Downregulated Dotplot
##==============================================================

if (nrow(reactome_down_df) > 0) {
  
  reactome_down_plot <- dotplot(
    reactome_down,
    showCategory = 20
  ) +
    ggtitle(
      "Reactome: Downregulated Genes"
    )
  
  print(reactome_down_plot)
  
  ggsave(
    filename = file.path(
      enrichment_figure_dir,
      "GSE149955_Reactome_Downregulated_Dotplot.pdf"
    ),
    plot = reactome_down_plot,
    width = 10,
    height = 8
  )
  
} else {
  
  message(
    "Skipping Reactome Downregulated dotplot."
  )
}

##==============================================================
## End of Script
##==============================================================
