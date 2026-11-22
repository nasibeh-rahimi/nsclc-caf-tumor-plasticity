############################################################
# 15_GSE33532_Tumor_Plasticity.R
#
# Project:
# Integrative Transcriptomic Analysis of CAF Activation
# and Tumor Plasticity in NSCLC
#
# Purpose:
# Develop and evaluate a tumor-level transcriptional axis
# associated with cancer-cell plasticity in NSCLC tumors.
#
# Biological Question:
# Do NSCLC tumors exhibit a continuous transcriptional
# state associated with adaptive tumor plasticity?
#
# Strategy:
# 1. Restrict analysis to tumor samples
# 2. Define literature-informed candidate plasticity programs
# 3. Evaluate their behavior in GSE33532
# 4. Derive a tumor-level plasticity-associated score
# 5. Identify genes positively and negatively associated
#    with the inferred axis
# 6. Validate the resulting state using independent
#    biological programs
#
# Important:
# The plasticity axis is derived independently of the CAF
# signature to avoid circularity.
#
# Input:
# expr_GSE33532_gene_final.rds
# GSE33532_sample_metadata.csv
#
# Output:
# Tumor-only expression matrix
# Candidate plasticity program scores
# Tumor plasticity-associated score
# Gene-wise correlations with the score
# QC and validation figures
#
# Interpretation:
# The score represents a transcriptional state associated
# with tumor plasticity. It does not directly measure
# dynamic phenotypic plasticity or prove temporal cell-state
# transitions because the dataset is cross-sectional bulk
# transcriptomics.
#
# Date:
# August 2026
############################################################
############################################################
# 1. Load project setup
############################################################

rm(list = ls())

source("scripts/00_project_setup.R")

cat("\n============================================\n")
cat("SCRIPT 15: GSE33532 Tumor Plasticity\n")
cat("============================================\n")

############################################################
# 2. Load GSE33532 expression and metadata
############################################################

expr_GSE33532 <- readRDS(
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

cat("\nFull expression matrix:\n")
cat(
  nrow(expr_GSE33532),
  "genes x",
  ncol(expr_GSE33532),
  "samples\n"
)

cat("\nMetadata rows:", nrow(metadata_GSE33532), "\n")

############################################################
# 3. Validate sample order
############################################################

sample_order_ok <- all(
  colnames(expr_GSE33532) ==
    metadata_GSE33532$sample
)

cat(
  "\nSample order matches metadata:",
  sample_order_ok,
  "\n"
)

stopifnot(sample_order_ok)

############################################################
# 4. Extract tumor samples
############################################################

tumor_idx <- metadata_GSE33532$group == "Tumor"

expr_tumor <- expr_GSE33532[
  ,
  tumor_idx,
  drop = FALSE
]

metadata_tumor <- metadata_GSE33532[
  tumor_idx,
  ,
  drop = FALSE
]

cat(
  "\nTumor expression matrix:",
  nrow(expr_tumor),
  "genes x",
  ncol(expr_tumor),
  "tumors\n"
)

############################################################
#5. Install/load msigdbr
if (!requireNamespace("msigdbr", quietly = TRUE)) {
  install.packages("msigdbr")
}

library(msigdbr)
library(dplyr)
library(tidyr)
library(purrr)

packageVersion("msigdbr")

msigdbr_version <- as.character(
  packageVersion("msigdbr")
)

cat(
  "\nmsigdbr version:",
  msigdbr_version,
  "\n"
)

############################################################
# 6. Candidate plasticity programs
############################################################

plasticity_program_info <- data.frame(
  
  program = c(
    "EMT_Byers76",
    "Stemness_Wong",
    "TGF_beta_Hallmark",
    "Hypoxia_Hallmark",
    "Therapy_Resistance_DTC_UP"
  ),
  
  role = c(
    "Core",
    "Core",
    "Supporting",
    "Supporting",
    "Independent validation"
  ),
  
  biological_rationale = c(
    "NSCLC-specific epithelial-mesenchymal state and phenotypic transition",
    "ESC-like stemness and dedifferentiation-associated tumor state",
    "Transcriptional response associated with TGF-beta signaling and microenvironmental cues",
    "Transcriptional adaptation associated with hypoxic stress",
    "Drug-tolerant state associated with adaptive survival under therapeutic pressure"
  ),
  
  expected_direction = c(
    "Mesenchymal direction",
    "Positive",
    "Positive",
    "Positive",
    "Positive"
  ),
  
  primary_plasticity_component = c(
    TRUE,
    TRUE,
    FALSE,
    FALSE,
    FALSE
  ),
  
  stringsAsFactors = FALSE
)

print(plasticity_program_info)


############################################################
# 7. Inspect tumor gene universe
############################################################

tumor_genes <- rownames(expr_tumor)

cat(
  "\nNumber of genes in tumor expression matrix:",
  length(tumor_genes),
  "\n"
)

cat(
  "\nFirst 20 genes:\n"
)

print(
  head(
    tumor_genes,
    20
  )
)

cat(
  "\nNumber of unique genes:",
  length(unique(tumor_genes)),
  "\n"
)

cat(
  "\nNumber of duplicated gene identifiers:",
  sum(duplicated(tumor_genes)),
  "\n"
)

############################################################
# 8. Check expression matrix properties
############################################################

summary(
  as.numeric(
    expr_tumor[1:100, 1:10]
  )
)

cat(
  "\nExpression matrix dimensions:\n"
)

print(
  dim(expr_tumor)
)

############################################################
############################################################
# 8B. Enforce tumor matrix QC
############################################################

stopifnot(
  nrow(expr_tumor) == 21355
)

stopifnot(
  ncol(expr_tumor) == 80
)

stopifnot(
  length(unique(tumor_genes)) == 21355
)

stopifnot(
  sum(duplicated(tumor_genes)) == 0
)

stopifnot(
  sum(is.na(expr_tumor)) == 0
)

cat(
  "\nTumor matrix QC passed.\n"
)
############################################################
# 9. Literature-derived gene-set framework
############################################################

plasticity_gene_sets_raw <- list(
  
  EMT_Byers76 = NULL,
  
  Stemness_Wong = NULL,
  
  TGF_beta_Hallmark = NULL,
  
  Hypoxia_Hallmark = NULL,
  
  Therapy_Resistance_DTC_UP = NULL
  
)

cat(
  "\nLiterature-derived gene-set framework initialized.\n"
)

print(
  names(plasticity_gene_sets_raw)
)


############################################################
# 10. Gene-set provenance
############################################################

plasticity_gene_set_provenance <- data.frame(
  
  program = c(
    "EMT_Byers76",
    "Stemness_Wong",
    "TGF_beta_Hallmark",
    "Hypoxia_Hallmark",
    "Therapy_Resistance_DTC_UP"
  ),
  
  source = c(
    "Byers et al. NSCLC EMT signature",
    "MSigDB: WONG_EMBRYONIC_STEM_CELL_CORE",
    "MSigDB: HALLMARK_TGF_BETA_SIGNALING",
    "MSigDB: HALLMARK_HYPOXIA",
    "Published NSCLC drug-tolerant-cell signature"
  ),
  
  source_type = c(
    "Published NSCLC signature",
    "MSigDB",
    "MSigDB",
    "MSigDB",
    "Published signature"
  ),
  
  role = c(
    "Core",
    "Core",
    "Supporting",
    "Supporting",
    "Independent validation"
  ),
  
  stringsAsFactors = FALSE
)

print(
  plasticity_gene_set_provenance
)

############################################################
# 11. Retrieve MSigDB gene sets
############################################################

cat("\n============================================\n")
cat("SECTION 11: RETRIEVE MSigDB GENE SETS\n")
cat("============================================\n")

############################################################
# 11A. Retrieve Hallmark collection
############################################################

msigdb_H <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

cat(
  "\nNumber of Hallmark gene sets:",
  length(unique(msigdb_H$gs_name)),
  "\n"
)

############################################################
# 11B. Retrieve C2 collection
############################################################

msigdb_C2 <- msigdbr(
  species = "Homo sapiens",
  collection = "C2"
)

cat(
  "\nNumber of C2 gene sets:",
  length(unique(msigdb_C2$gs_name)),
  "\n"
)

############################################################
# 12. Verify exact MSigDB gene-set names
############################################################

target_msigdbr_sets <- c(
  "WONG_EMBRYONIC_STEM_CELL_CORE",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_HYPOXIA"
)

cat(
  "\nChecking requested MSigDB gene sets:\n\n"
)

for (gs in target_msigdbr_sets) {
  
  found_H <- gs %in% unique(msigdb_H$gs_name)
  found_C2 <- gs %in% unique(msigdb_C2$gs_name)
  
  cat(
    gs,
    "\n  Hallmark collection:",
    found_H,
    "\n  C2 collection:",
    found_C2,
    "\n\n"
  )
}

############################################################
# 13. Extract MSigDB gene sets
############################################################

Stemness_Wong <- msigdb_C2 %>%
  filter(
    gs_name == "WONG_EMBRYONIC_STEM_CELL_CORE"
  ) %>%
  pull(gene_symbol) %>%
  unique()

TGF_beta_Hallmark <- msigdb_H %>%
  filter(
    gs_name == "HALLMARK_TGF_BETA_SIGNALING"
  ) %>%
  pull(gene_symbol) %>%
  unique()

Hypoxia_Hallmark <- msigdb_H %>%
  filter(
    gs_name == "HALLMARK_HYPOXIA"
  ) %>%
  pull(gene_symbol) %>%
  unique()

############################################################
# 14. MSigDB gene-set QC
############################################################

msigdbr_gene_set_sizes <- data.frame(
  
  program = c(
    "Stemness_Wong",
    "TGF_beta_Hallmark",
    "Hypoxia_Hallmark"
  ),
  
  gene_count = c(
    length(Stemness_Wong),
    length(TGF_beta_Hallmark),
    length(Hypoxia_Hallmark)
  ),
  
  stringsAsFactors = FALSE
)

print(msigdbr_gene_set_sizes)

############################################################
# Check that none are empty
############################################################

stopifnot(
  length(Stemness_Wong) > 0,
  length(TGF_beta_Hallmark) > 0,
  length(Hypoxia_Hallmark) > 0
)

############################################################
# Check for duplicated genes
############################################################

cat("\nDuplicates:\n")

cat(
  "Stemness:",
  length(Stemness_Wong) -
    length(unique(Stemness_Wong)),
  "\n"
)

cat(
  "TGF-beta:",
  length(TGF_beta_Hallmark) -
    length(unique(TGF_beta_Hallmark)),
  "\n"
)

cat(
  "Hypoxia:",
  length(Hypoxia_Hallmark) -
    length(unique(Hypoxia_Hallmark)),
  "\n"
)

############################################################
# 15. Harmonize MSigDB gene sets to tumor gene universe
############################################################

harmonize_gene_set <- function(
    gene_set,
    tumor_genes
) {
  
  gene_set <- unique(gene_set)
  
  present <- intersect(
    gene_set,
    tumor_genes
  )
  
  missing <- setdiff(
    gene_set,
    tumor_genes
  )
  
  return(
    list(
      original = gene_set,
      present = present,
      missing = missing
    )
  )
}


Stemness_Wong_harmonized <- harmonize_gene_set(
  Stemness_Wong,
  tumor_genes
)

TGF_beta_Hallmark_harmonized <- harmonize_gene_set(
  TGF_beta_Hallmark,
  tumor_genes
)

Hypoxia_Hallmark_harmonized <- harmonize_gene_set(
  Hypoxia_Hallmark,
  tumor_genes
)


############################################################
# 16. Gene-set coverage
############################################################

msigdbr_coverage <- data.frame(
  
  program = c(
    "Stemness_Wong",
    "TGF_beta_Hallmark",
    "Hypoxia_Hallmark"
  ),
  
  original_n = c(
    
    length(
      Stemness_Wong_harmonized$original
    ),
    
    length(
      TGF_beta_Hallmark_harmonized$original
    ),
    
    length(
      Hypoxia_Hallmark_harmonized$original
    )
  ),
  
  present_n = c(
    
    length(
      Stemness_Wong_harmonized$present
    ),
    
    length(
      TGF_beta_Hallmark_harmonized$present
    ),
    
    length(
      Hypoxia_Hallmark_harmonized$present
    )
  ),
  
  missing_n = c(
    
    length(
      Stemness_Wong_harmonized$missing
    ),
    
    length(
      TGF_beta_Hallmark_harmonized$missing
    ),
    
    length(
      Hypoxia_Hallmark_harmonized$missing
    )
  ),
  
  stringsAsFactors = FALSE
)

msigdbr_coverage$coverage_pct <-
  round(
    100 *
      msigdbr_coverage$present_n /
      msigdbr_coverage$original_n,
    2
  )

print(msigdbr_coverage)

############################################################
# 16B. Inspect missing genes
############################################################

cat("\nMissing Stemness genes:\n")
print(
  Stemness_Wong_harmonized$missing
)

cat("\nMissing TGF-beta genes:\n")
print(
  TGF_beta_Hallmark_harmonized$missing
)

cat("\nMissing Hypoxia genes:\n")
print(
  Hypoxia_Hallmark_harmonized$missing
)

############################################################
# 17. Store harmonized MSigDB gene sets
############################################################

plasticity_gene_sets <- list(
  
  Stemness_Wong =
    Stemness_Wong_harmonized$present,
  
  TGF_beta_Hallmark =
    TGF_beta_Hallmark_harmonized$present,
  
  Hypoxia_Hallmark =
    Hypoxia_Hallmark_harmonized$present
)

############################################################
# Verify harmonized set sizes
############################################################

cat("\nHarmonized MSigDB gene-set sizes:\n")

print(
  sapply(
    plasticity_gene_sets,
    length
  )
)

############################################################
# Save harmonized MSigDB sets
############################################################

saveRDS(
  plasticity_gene_sets,
  file.path(
    plasticity_dir,
    "plasticity_gene_sets_MSigDB_harmonized.rds"
  )
)

############################################################
# Save coverage table
############################################################

write.csv(
  msigdbr_coverage,
  file.path(
    table_dir,
    "plasticity_MSigDB_gene_set_coverage.csv"
  ),
  row.names = FALSE
)

cat(
  "\nHarmonized MSigDB gene sets saved.\n"
)

############################################################
############################################################
# 18. Byers 76-gene NSCLC EMT signature
############################################################
#
# Byers et al. developed a 76-gene EMT signature in NSCLC.
#
# The original Supplemental Table 1 contains:
#   96 Affymetrix probes
#   representing 76 unique genes
#
# The published 76GS score is a weighted sum of expression,
# where weights are based on correlation with CDH1.
#
# IMPORTANT:
# At this stage we only define and harmonize the gene set.
# We are NOT calculating the EMT score yet.
############################################################

EMT_Byers76 <- c(
  
  "ANKRD22",
  "ANTXR2",
  "AP1M2",
  "AXL",
  "BSPRY",
  "C1orf116",
  "C1orf172",
  "C3orf21",
  "CARD6",
  "CDH1",
  "CDH3",
  "CDS1",
  "CLDN4",
  "CLDN7",
  "CRB3",
  "DSP",
  "ELMO3",
  "ENPP5",
  "EPB41L5",
  "EPHA1",
  "EPN3",
  "EPPK1",
  "ERBB3",
  "EVPL",
  "F11R",
  "FN1",
  "FXYD3",
  "GALNT3",
  "GALNT5",
  "GPR110",
  "GPR56",
  "GRHL1",
  "GRHL2",
  "HNMT",
  "INADL",
  "ITGB6",
  "KLC3",
  "KRT19",
  "KRTCAP3",
  "LIX1L",
  "LRRC54",
  "MAL2",
  "MAPK13",
  "MMP2",
  "MPP7",
  "MPZL2",
  "MTAC2D1",
  "MUC1",
  "NRP1",
  "PPARG",
  "PRR5",
  "PRSS22",
  "PRSS8",
  "RAB25",
  "RBM35A",
  "RBPMS",
  "S100A14",
  "SCNN1A",
  "SERINC2",
  "SH3YL1",
  "SHROOM3",
  "SPINT2",
  "SSH3",
  "ST14",
  "STAP2",
  "TACSTD1",
  "TACSTD2",
  "TGFBI",
  "TJP3",
  "TMC4",
  "TMEM125",
  "TMEM30B",
  "TMEM45B",
  "TNFRSF21",
  "VIM",
  "ZEB1"
)

############################################################
# Validate the signature
############################################################

cat(
  "\nNumber of Byers genes:",
  length(EMT_Byers76),
  "\n"
)

cat(
  "Number of unique Byers genes:",
  length(unique(EMT_Byers76)),
  "\n"
)

stopifnot(
  length(EMT_Byers76) == 76,
  length(unique(EMT_Byers76)) == 76
)

cat(
  "\nByers 76-gene signature validated.\n"
)

############################################################
# 19. Harmonize historical Byers gene symbols
############################################################

Byers_symbol_map <- c(
  
  "C1orf172" = "KDF1",
  "C3orf21"  = "XXYLT1",
  "GPR110"   = "ADGRF1",
  "GPR56"    = "ADGRG1",
  "INADL"    = "PATJ",
  "LRRC54"   = "TSKU",
  "MTAC2D1"  = "TC2N",
  "RBM35A"   = "ESRP1",
  "TACSTD1"  = "EPCAM"
  
)

############################################################
# Check whether the replacement symbols exist in our matrix
############################################################

Byers_mapping_check <- data.frame(
  
  original_symbol = names(Byers_symbol_map),
  
  harmonized_symbol = unname(Byers_symbol_map),
  
  present_in_array =
    unname(Byers_symbol_map) %in% tumor_genes,
  
  stringsAsFactors = FALSE
  
)

print(
  Byers_mapping_check
)

############################################################
# Check whether any replacement symbols are duplicated
############################################################

cat(
  "\nNumber of duplicated harmonized symbols:",
  sum(
    duplicated(
      unname(Byers_symbol_map)
    )
  ),
  "\n"
)


############################################################
# 20. Create harmonized Byers 76-gene signature
############################################################

EMT_Byers76_harmonized <- EMT_Byers76

############################################################
# Replace historical symbols with current HGNC symbols
############################################################

EMT_Byers76_harmonized[
  EMT_Byers76_harmonized %in%
    names(Byers_symbol_map)
] <-
  Byers_symbol_map[
    EMT_Byers76_harmonized[
      EMT_Byers76_harmonized %in%
        names(Byers_symbol_map)
    ]
  ]

############################################################
# Validate
############################################################

cat(
  "\nOriginal Byers gene count:",
  length(EMT_Byers76),
  "\n"
)

cat(
  "Harmonized Byers gene count:",
  length(EMT_Byers76_harmonized),
  "\n"
)

cat(
  "Unique harmonized genes:",
  length(
    unique(
      EMT_Byers76_harmonized
    )
  ),
  "\n"
)

stopifnot(
  length(EMT_Byers76_harmonized) == 76,
  length(
    unique(
      EMT_Byers76_harmonized
    )
  ) == 76
)

############################################################
# 21. Final Byers coverage after symbol harmonization
############################################################

EMT_Byers76_present <- intersect(
  EMT_Byers76_harmonized,
  tumor_genes
)

EMT_Byers76_missing <- setdiff(
  EMT_Byers76_harmonized,
  tumor_genes
)

EMT_Byers76_coverage_final <- data.frame(
  
  program = "EMT_Byers76",
  
  original_n = 76,
  
  present_n =
    length(
      EMT_Byers76_present
    ),
  
  missing_n =
    length(
      EMT_Byers76_missing
    ),
  
  coverage_pct =
    round(
      100 *
        length(
          EMT_Byers76_present
        ) /
        76,
      1
    ),
  
  stringsAsFactors = FALSE
)

print(
  EMT_Byers76_coverage_final
)

cat(
  "\nRemaining missing Byers genes after harmonization:\n"
)

print(
  EMT_Byers76_missing
)


############################################################
# 22. Update harmonized gene-set collection
############################################################

# Create the harmonized collection from the already
# harmonized MSigDB gene sets

plasticity_gene_sets_harmonized <-
  plasticity_gene_sets

# Add the harmonized Byers EMT signature

plasticity_gene_sets_harmonized$EMT_Byers76 <-
  EMT_Byers76_present

############################################################
# Reorder programs for clarity
############################################################

plasticity_gene_sets_harmonized <-
  plasticity_gene_sets_harmonized[
    c(
      "EMT_Byers76",
      "Stemness_Wong",
      "TGF_beta_Hallmark",
      "Hypoxia_Hallmark"
    )
  ]

############################################################
# Check final harmonized gene-set sizes
############################################################

cat(
  "\nFinal harmonized gene-set sizes:\n"
)

print(
  sapply(
    plasticity_gene_sets_harmonized,
    length
  )
)

############################################################
# Final coverage check
############################################################

cat(
  "\nFinal Byers EMT coverage:\n"
)

print(
  EMT_Byers76_coverage_final
)

############################################################
# Verify that all final gene sets contain
# only genes present in the tumor expression matrix
############################################################

for (program in names(plasticity_gene_sets_harmonized)) {
  
  stopifnot(
    all(
      plasticity_gene_sets_harmonized[[program]] %in%
        tumor_genes
    )
  )
  
}

cat(
  "\nAll harmonized gene sets are present in the tumor gene universe.\n"
)


############################################################
# 23. Finalize harmonized plasticity gene-set collection
############################################################

plasticity_gene_sets_harmonized <- list(
  
  EMT_Byers76 =
    EMT_Byers76_present,
  
  Stemness_Wong =
    Stemness_Wong_harmonized$present,
  
  TGF_beta_Hallmark =
    TGF_beta_Hallmark_harmonized$present,
  
  Hypoxia_Hallmark =
    Hypoxia_Hallmark_harmonized$present
  
)

cat("\nFinal harmonized gene-set sizes:\n")

print(
  sapply(
    plasticity_gene_sets_harmonized,
    length
  )
)

############################################################
# QC
############################################################

stopifnot(
  length(plasticity_gene_sets_harmonized$EMT_Byers76) == 76,
  length(plasticity_gene_sets_harmonized$Stemness_Wong) == 334,
  length(plasticity_gene_sets_harmonized$TGF_beta_Hallmark) == 54,
  length(plasticity_gene_sets_harmonized$Hypoxia_Hallmark) == 199
)

cat("\nFinal plasticity gene-set QC passed.\n")


############################################################
# 24. Prepare Byers EMT expression matrix
############################################################

EMT_Byers76_expression <- expr_tumor[
  plasticity_gene_sets_harmonized$EMT_Byers76,
  ,
  drop = FALSE
]

cat(
  "\nByers EMT expression matrix:",
  nrow(EMT_Byers76_expression),
  "genes x",
  ncol(EMT_Byers76_expression),
  "tumors\n"
)

stopifnot(
  nrow(EMT_Byers76_expression) == 76,
  ncol(EMT_Byers76_expression) == 80
)


############################################################
# 25. Check expression of CDH1
############################################################

stopifnot(
  "CDH1" %in% rownames(
    EMT_Byers76_expression
  )
)

CDH1_expression <- as.numeric(
  expr_tumor["CDH1", ]
)

############################################################
# 26. Calculate dataset-specific Byers weights
############################################################

Byers76_weights <- sapply(
  
  rownames(
    EMT_Byers76_expression
  ),
  
  function(gene) {
    
    cor(
      as.numeric(
        EMT_Byers76_expression[gene, ]
      ),
      CDH1_expression,
      method = "pearson",
      use = "pairwise.complete.obs"
    )
    
  }
)

Byers76_weight_table <- data.frame(
  
  gene = names(Byers76_weights),
  
  weight_CDH1 =
    as.numeric(Byers76_weights),
  
  stringsAsFactors = FALSE
)

Byers76_weight_table <- Byers76_weight_table[
  order(
    Byers76_weight_table$weight_CDH1,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

cat("\nByers weight summary:\n")

print(
  summary(
    Byers76_weight_table$weight_CDH1
  )
)

############################################################
# 27. Check weights
############################################################

cat(
  "\nPositive weights:",
  sum(
    Byers76_weight_table$weight_CDH1 > 0
  ),
  "\n"
)

cat(
  "Negative weights:",
  sum(
    Byers76_weight_table$weight_CDH1 < 0
  ),
  "\n"
)

cat(
  "NA weights:",
  sum(
    is.na(
      Byers76_weight_table$weight_CDH1
    )
  ),
  "\n"
)

stopifnot(
  sum(
    is.na(
      Byers76_weight_table$weight_CDH1
    )
  ) == 0
)


############################################################
# 28. Calculate epithelial-oriented 76GS
############################################################

Byers76_gene_order <- rownames(
  EMT_Byers76_expression
)

Byers76_weights_ordered <- Byers76_weights[
  Byers76_gene_order
]

stopifnot(
  identical(
    names(Byers76_weights_ordered),
    Byers76_gene_order
  )
)

Byers76_score_epithelial <- colSums(
  
  sweep(
    EMT_Byers76_expression,
    MARGIN = 1,
    STATS = Byers76_weights_ordered,
    FUN = "*"
  )
  
)

############################################################
# 29. Center the epithelial score
############################################################

Byers76_score_epithelial_centered <-
  Byers76_score_epithelial -
  mean(
    Byers76_score_epithelial
  )

cat(
  "\nMean centered epithelial score:",
  mean(
    Byers76_score_epithelial_centered
  ),
  "\n"
)


############################################################
# 30. Convert to mesenchymal-oriented EMT score
############################################################

EMT_Byers76_score <-
  -Byers76_score_epithelial_centered

cat(
  "\nSummary of mesenchymal-oriented EMT score:\n"
)

print(
  summary(
    EMT_Byers76_score
  )
)

cat(
  "\nMean EMT score:",
  mean(
    EMT_Byers76_score
  ),
  "\n"
)

############################################################
# 31. Check score direction
############################################################

cat(
  "\nCorrelation between EMT score and CDH1:",
  cor(
    EMT_Byers76_score,
    CDH1_expression,
    method = "pearson"
  ),
  "\n"
)

cat(
  "\nExpected: negative correlation because higher EMT\n",
  "should correspond to lower epithelial CDH1 expression.\n"
)


############################################################
# 32. Prepare GSVA/ssGSEA scoring
############################################################

cat("\n============================================\n")
cat("SECTION 32: SAMPLE-LEVEL PROGRAM SCORING\n")
cat("============================================\n")

############################################################
# Check GSVA installation
############################################################

if (!requireNamespace("GSVA", quietly = TRUE)) {
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  BiocManager::install(
    "GSVA",
    ask = FALSE,
    update = FALSE
  )
}

library(GSVA)

cat(
  "\nGSVA version:",
  as.character(
    packageVersion("GSVA")
  ),
  "\n"
)

############################################################
# 33. Prepare MSigDB gene sets for scoring
############################################################

program_gene_sets <- list(
  
  Stemness_Wong =
    plasticity_gene_sets_harmonized$Stemness_Wong,
  
  TGF_beta_Hallmark =
    plasticity_gene_sets_harmonized$TGF_beta_Hallmark,
  
  Hypoxia_Hallmark =
    plasticity_gene_sets_harmonized$Hypoxia_Hallmark
  
)

############################################################
# Verify gene-set sizes
############################################################

cat("\nGene-set sizes for scoring:\n")

print(
  sapply(
    program_gene_sets,
    length
  )
)

############################################################
# Verify every gene exists in expression matrix
############################################################

program_coverage_check <- data.frame(
  
  program = names(program_gene_sets),
  
  gene_set_size =
    sapply(
      program_gene_sets,
      length
    ),
  
  genes_present =
    sapply(
      program_gene_sets,
      function(gs) {
        sum(
          gs %in%
            rownames(expr_tumor)
        )
      }
    ),
  
  stringsAsFactors = FALSE
)

program_coverage_check$coverage_pct <-
  round(
    100 *
      program_coverage_check$genes_present /
      program_coverage_check$gene_set_size,
    1
  )

print(
  program_coverage_check
)

stopifnot(
  all(
    program_coverage_check$genes_present ==
      program_coverage_check$gene_set_size
  )
)

############################################################
# 34. Run ssGSEA
############################################################
#
# ssGSEA is used here for sample-level scoring of
# transcriptional programs.
#
# Important:
# ssgseaParam() is used specifically for the ssGSEA method.
############################################################

cat(
  "\nRunning ssGSEA for Stemness, TGF-beta and Hypoxia...\n"
)

ssgsea_param <- ssgseaParam(
  
  expr_tumor,
  
  program_gene_sets
  
)

program_scores_ssgsea <- gsva(
  
  ssgsea_param
  
)

############################################################
# Check result
############################################################

cat(
  "\nssGSEA score matrix dimensions:\n"
)

print(
  dim(
    program_scores_ssgsea
  )
)

cat(
  "\nPrograms:\n"
)

print(
  rownames(
    program_scores_ssgsea
  )
)

cat(
  "\nTumors:\n"
)

print(
  head(
    colnames(
      program_scores_ssgsea
    )
  )
)

############################################################
# Validate dimensions
############################################################

stopifnot(
  nrow(program_scores_ssgsea) == 3,
  ncol(program_scores_ssgsea) == 80
)

stopifnot(
  identical(
    colnames(program_scores_ssgsea),
    colnames(expr_tumor)
  )
)

cat(
  "\nssGSEA QC passed.\n"
)

############################################################
# 35. Summarize program scores
############################################################

program_score_summary <- data.frame(
  
  program =
    rownames(
      program_scores_ssgsea
    ),
  
  min =
    apply(
      program_scores_ssgsea,
      1,
      min
    ),
  
  Q1 =
    apply(
      program_scores_ssgsea,
      1,
      quantile,
      probs = 0.25
    ),
  
  median =
    apply(
      program_scores_ssgsea,
      1,
      median
    ),
  
  mean =
    apply(
      program_scores_ssgsea,
      1,
      mean
    ),
  
  Q3 =
    apply(
      program_scores_ssgsea,
      1,
      quantile,
      probs = 0.75
    ),
  
  max =
    apply(
      program_scores_ssgsea,
      1,
      max
    ),
  
  SD =
    apply(
      program_scores_ssgsea,
      1,
      sd
    ),
  
  stringsAsFactors = FALSE
)

print(
  program_score_summary
)

###########################################
############################################################
# 36. Correlation among candidate plasticity programs
############################################################

cat("\n============================================\n")
cat("SECTION 36: CORRELATION AMONG PLASTICITY PROGRAMS\n")
cat("============================================\n")

############################################################
# 36A. Construct program score matrix
############################################################

candidate_program_scores <- rbind(
  
  EMT_Byers76 =
    EMT_Byers76_score,
  
  Stemness_Wong =
    as.numeric(
      program_scores_ssgsea["Stemness_Wong", ]
    ),
  
  TGF_beta_Hallmark =
    as.numeric(
      program_scores_ssgsea["TGF_beta_Hallmark", ]
    ),
  
  Hypoxia_Hallmark =
    as.numeric(
      program_scores_ssgsea["Hypoxia_Hallmark", ]
    )
)
############################################################
# Extract named EMT score for Script 16
############################################################

EMT_score <- as.numeric(
  candidate_program_scores["EMT_Byers76", ]
)

names(EMT_score) <- colnames(candidate_program_scores)

############################################################
# Validate dimensions
############################################################

cat(
  "\nCandidate program score matrix:\n"
)

print(
  dim(candidate_program_scores)
)

print(
  rownames(candidate_program_scores)
)

stopifnot(
  nrow(candidate_program_scores) == 4,
  ncol(candidate_program_scores) == 80
)

############################################################
# Check sample order
############################################################

stopifnot(
  identical(
    colnames(candidate_program_scores),
    colnames(expr_tumor)
  )
)

cat(
  "\nProgram score matrix QC passed.\n"
)

############################################################
# 36B. Pearson correlation matrix
############################################################

candidate_cor_pearson <- cor(
  t(candidate_program_scores),
  method = "pearson",
  use = "pairwise.complete.obs"
)

cat(
  "\nPearson correlation matrix:\n"
)

print(
  round(
    candidate_cor_pearson,
    3
  )
)

############################################################
# 36C. Spearman correlation matrix
############################################################

candidate_cor_spearman <- cor(
  t(candidate_program_scores),
  method = "spearman",
  use = "pairwise.complete.obs"
)

cat(
  "\nSpearman correlation matrix:\n"
)

print(
  round(
    candidate_cor_spearman,
    3
  )
)

############################################################
# 36D. Compare Pearson and Spearman
############################################################

cat(
  "\nDifference between Pearson and Spearman correlations:\n"
)

print(
  round(
    candidate_cor_pearson -
      candidate_cor_spearman,
    3
  )
)

############################################################
# 36E. Check score distributions
############################################################

cat(
  "\nProgram score summaries:\n"
)

for (program in rownames(candidate_program_scores)) {
  
  cat(
    "\n---",
    program,
    "---\n"
  )
  
  print(
    summary(
      candidate_program_scores[
        program,
      ]
    )
  )
}

############################################################
# 36F. Save correlation matrices
############################################################

write.csv(
  candidate_cor_pearson,
  file.path(
    table_dir,
    "plasticity_candidate_correlations_Pearson.csv"
  ),
  row.names = TRUE
)

write.csv(
  candidate_cor_spearman,
  file.path(
    table_dir,
    "plasticity_candidate_correlations_Spearman.csv"
  ),
  row.names = TRUE
)

cat(
  "\nCorrelation matrices saved.\n"
)

############################################################
# 37. EMT–Stemness relationship analysis
############################################################
#
# Purpose:
# Examine the relationship between the two candidate
# core plasticity programs:
#
#   1. EMT_Byers76
#   2. Stemness_Wong
#
# We do NOT construct a composite plasticity score yet.
#
# Questions:
#   - Is the relationship statistically significant?
#   - Is it approximately linear?
#   - Is it robust to rank-based correlation?
#   - Are there obvious outliers?
#   - Does the cohort appear to contain distinct states?
############################################################

cat("\n============================================\n")
cat("SECTION 37: EMT-STEMNESS RELATIONSHIP\n")
cat("============================================\n")

############################################################
# 37A. Extract scores
############################################################

EMT_score <- as.numeric(
  candidate_program_scores[
    "EMT_Byers76",
  ]
)
names(EMT_score) <- colnames(candidate_program_scores)


Stemness_score <- as.numeric(
  candidate_program_scores[
    "Stemness_Wong",
  ]
)

############################################################
# Validate
############################################################

stopifnot(
  length(EMT_score) == 80,
  length(Stemness_score) == 80
)

stopifnot(
  !any(is.na(EMT_score)),
  !any(is.na(Stemness_score))
)

############################################################
# 37B. Pearson correlation test
############################################################

EMT_Stemness_pearson <- cor.test(
  EMT_score,
  Stemness_score,
  method = "pearson"
)

cat(
  "\nPearson correlation:\n"
)

print(
  EMT_Stemness_pearson
)

############################################################
# 37C. Spearman correlation test
############################################################

EMT_Stemness_spearman <- cor.test(
  EMT_score,
  Stemness_score,
  method = "spearman",
  exact = FALSE
)

cat(
  "\nSpearman correlation:\n"
)

print(
  EMT_Stemness_spearman
)

############################################################
# 37D. Extract key statistics
############################################################

EMT_Stemness_statistics <- data.frame(
  
  comparison = "EMT_Byers76_vs_Stemness_Wong",
  
  pearson_r =
    unname(
      EMT_Stemness_pearson$estimate
    ),
  
  pearson_p =
    EMT_Stemness_pearson$p.value,
  
  spearman_rho =
    unname(
      EMT_Stemness_spearman$estimate
    ),
  
  spearman_p =
    EMT_Stemness_spearman$p.value,
  
  stringsAsFactors = FALSE
)

print(
  EMT_Stemness_statistics
)

############################################################
# 37E. Linear regression
############################################################
#
# Model:
#
# Stemness = beta0 + beta1 * EMT + error
#
############################################################

EMT_Stemness_lm <- lm(
  Stemness_score ~ EMT_score
)

cat(
  "\nLinear regression:\n"
)

print(
  summary(
    EMT_Stemness_lm
  )
)

############################################################
# Extract regression statistics
############################################################

EMT_Stemness_lm_summary <- summary(
  EMT_Stemness_lm
)

EMT_Stemness_regression_statistics <- data.frame(
  
  intercept =
    coef(
      EMT_Stemness_lm
    )[1],
  
  EMT_beta =
    coef(
      EMT_Stemness_lm
    )[2],
  
  R_squared =
    EMT_Stemness_lm_summary$r.squared,
  
  adjusted_R_squared =
    EMT_Stemness_lm_summary$adj.r.squared,
  
  regression_p =
    coef(
      EMT_Stemness_lm_summary
    )[2, 4],
  
  stringsAsFactors = FALSE
)

print(
  EMT_Stemness_regression_statistics
)

############################################################
# 37F. Create analysis data frame
############################################################

EMT_Stemness_df <- data.frame(
  
  sample =
    colnames(
      expr_tumor
    ),
  
  EMT_Byers76 =
    EMT_score,
  
  Stemness_Wong =
    Stemness_score,
  
  stringsAsFactors = FALSE
)

############################################################
# Add patient identifier if available
############################################################

if (
  "patient" %in%
  colnames(metadata_tumor)
) {
  
  EMT_Stemness_df$patient <-
    metadata_tumor$patient
  
}

############################################################
# Inspect
############################################################

cat(
  "\nFirst rows of EMT-Stemness dataset:\n"
)

print(
  head(
    EMT_Stemness_df
  )
)

############################################################
# 37G. Identify potential outliers
############################################################
#
# Standardized scores are used only for QC here.
# They are NOT used to construct the final axis.
############################################################

EMT_Stemness_df$EMT_z <- as.numeric(
  scale(
    EMT_Stemness_df$EMT_Byers76
  )
)

EMT_Stemness_df$Stemness_z <- as.numeric(
  scale(
    EMT_Stemness_df$Stemness_Wong
  )
)

############################################################
# Define potential extreme observations
#
# |z| > 3 is a conservative QC threshold.
############################################################

potential_outliers <- EMT_Stemness_df[
  abs(
    EMT_Stemness_df$EMT_z
  ) > 3 |
    abs(
      EMT_Stemness_df$Stemness_z
    ) > 3,
  ,
  drop = FALSE
]

cat(
  "\nPotential extreme observations (|z| > 3):\n"
)

print(
  potential_outliers
)

cat(
  "\nNumber of potential extreme observations:",
  nrow(
    potential_outliers
  ),
  "\n"
)

############################################################
# 37H. Save analysis tables
############################################################

write.csv(
  EMT_Stemness_statistics,
  file.path(
    table_dir,
    "plasticity_EMT_Stemness_correlation_statistics.csv"
  ),
  row.names = FALSE
)

write.csv(
  EMT_Stemness_regression_statistics,
  file.path(
    table_dir,
    "plasticity_EMT_Stemness_regression_statistics.csv"
  ),
  row.names = FALSE
)

write.csv(
  EMT_Stemness_df,
  file.path(
    table_dir,
    "plasticity_EMT_Stemness_tumor_scores.csv"
  ),
  row.names = FALSE
)

############################################################
# 37I. Scatterplot
############################################################

library(ggplot2)

EMT_Stemness_plot <- ggplot(
  EMT_Stemness_df,
  aes(
    x = EMT_Byers76,
    y = Stemness_Wong
  )
) +
  
  geom_point(
    size = 2.5,
    alpha = 0.8
  ) +
  
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  
  labs(
    title = "Relationship between EMT and Wong Stemness",
    x = "Byers 76GS EMT score (mesenchymal-oriented)",
    y = "Wong embryonic stemness score"
  ) +
  
  theme_classic()

print(
  EMT_Stemness_plot
)

############################################################
# 37J. Save plot
############################################################

ggsave(
  filename = file.path(
    figure_dir,
    "plasticity_EMT_vs_Stemness.png"
  ),
  plot = EMT_Stemness_plot,
  width = 7,
  height = 6,
  dpi = 300
)

cat(
  "\nEMT-Stemness plot saved.\n"
)

############################################################
# 37K. QC conclusion
############################################################

cat("\n============================================\n")
cat("SECTION 37 COMPLETE\n")
cat("============================================\n")

cat(
  "\nPearson r:",
  round(
    EMT_Stemness_pearson$estimate,
    3
  ),
  "\n"
)

cat(
  "Pearson p:",
  signif(
    EMT_Stemness_pearson$p.value,
    3
  ),
  "\n"
)

cat(
  "Spearman rho:",
  round(
    EMT_Stemness_spearman$estimate,
    3
  ),
  "\n"
)

cat(
  "Spearman p:",
  signif(
    EMT_Stemness_spearman$p.value,
    3
  ),
  "\n"
)

cat(
  "R-squared:",
  round(
    EMT_Stemness_lm_summary$r.squared,
    3
  ),
  "\n"
)

cat(
  "Potential |z| > 3 observations:",
  nrow(
    potential_outliers
  ),
  "\n"
)


############################################################
# 38. Exploratory PCA of core plasticity programs
############################################################
#
# Purpose:
# Explore whether EMT and stemness occupy a shared
# transcriptional axis.
#
# IMPORTANT:
# This PCA is exploratory.
# PC1 is NOT yet defined as the final tumor plasticity score.
#
# The PCA will be revisited after independent validation
# of stemness using the Liu et al. stemness signatures.
############################################################

cat("\n============================================\n")
cat("SECTION 38: EXPLORATORY CORE PROGRAM PCA\n")
cat("============================================\n")

############################################################
# 38A. Construct core program matrix
############################################################

core_program_scores_df <- as.data.frame(
  t(
    candidate_program_scores[
      c(
        "EMT_Byers76",
        "Stemness_Wong"
      ),
      ,
      drop = FALSE
    ]
  )
)

############################################################
# Validate dimensions
############################################################

cat(
  "\nCore program score data frame:\n"
)

print(
  dim(
    core_program_scores_df
  )
)

stopifnot(
  nrow(core_program_scores_df) == 80,
  ncol(core_program_scores_df) == 2
)

############################################################
# 38B. Standardize core programs
############################################################

core_program_scores_scaled <- scale(
  core_program_scores_df
)

############################################################
# 38C. PCA
############################################################

core_pca <- prcomp(
  
  core_program_scores_scaled,
  
  center = FALSE,
  
  scale. = FALSE
  
)

############################################################
# PCA summary
############################################################

cat(
  "\nPCA summary:\n"
)

print(
  summary(
    core_pca
  )
)

############################################################
# 38D. PCA rotation
############################################################

cat(
  "\nPCA rotation:\n"
)

print(
  core_pca$rotation
)

############################################################
# 38E. PCA scores
############################################################

core_pca_scores <- as.data.frame(
  core_pca$x
)

core_pca_scores$sample <- rownames(
  core_pca_scores
)

############################################################
# Check sample order
############################################################

stopifnot(
  identical(
    core_pca_scores$sample,
    colnames(expr_tumor)
  )
)

############################################################
# 38F. Variance explained
############################################################

core_pca_variance <- (
  core_pca$sdev^2 /
    sum(
      core_pca$sdev^2
    )
)

core_pca_variance_df <- data.frame(
  
  PC = paste0(
    "PC",
    seq_along(
      core_pca_variance
    )
  ),
  
  variance_explained =
    core_pca_variance,
  
  variance_explained_pct =
    100 *
    core_pca_variance,
  
  stringsAsFactors = FALSE
  
)

print(
  core_pca_variance_df
)

############################################################
# 38G. Correlations with biological reference markers
############################################################

core_pca_vs_CDH1 <- cor.test(
  
  core_pca_scores$PC1,
  
  CDH1_expression,
  
  method = "pearson"
  
)

cat(
  "\nPC1 vs CDH1:\n"
)

print(
  core_pca_vs_CDH1
)

############################################################
# 38H. Correlations with supporting programs
############################################################

core_pca_vs_supporting <- data.frame(
  
  variable = c(
    "CDH1",
    "TGF_beta_Hallmark",
    "Hypoxia_Hallmark"
  ),
  
  PC1_pearson_r = c(
    
    cor(
      core_pca_scores$PC1,
      CDH1_expression,
      method = "pearson"
    ),
    
    cor(
      core_pca_scores$PC1,
      candidate_program_scores[
        "TGF_beta_Hallmark",
      ]
    ),
    
    cor(
      core_pca_scores$PC1,
      candidate_program_scores[
        "Hypoxia_Hallmark",
      ]
    )
    
  ),
  
  PC2_pearson_r = c(
    
    cor(
      core_pca_scores$PC2,
      CDH1_expression,
      method = "pearson"
    ),
    
    cor(
      core_pca_scores$PC2,
      candidate_program_scores[
        "TGF_beta_Hallmark",
      ]
    ),
    
    cor(
      core_pca_scores$PC2,
      candidate_program_scores[
        "Hypoxia_Hallmark",
      ]
    )
    
  ),
  
  stringsAsFactors = FALSE
  
)

print(
  core_pca_vs_supporting
)

############################################################
# 38I. Save exploratory PCA results
############################################################

write.csv(
  
  core_pca_variance_df,
  
  file.path(
    table_dir,
    "plasticity_core_PCA_variance_explained.csv"
  ),
  
  row.names = FALSE
  
)

write.csv(
  
  core_pca_vs_supporting,
  
  file.path(
    table_dir,
    "plasticity_core_PCA_biological_correlations.csv"
  ),
  
  row.names = FALSE
  
)

write.csv(
  
  core_pca_scores,
  
  file.path(
    table_dir,
    "plasticity_core_PCA_scores.csv"
  ),
  
  row.names = FALSE
  
)

cat(
  "\nExploratory PCA results saved.\n"
)

############################################################
# 38J. Interpretation note
############################################################

cat(
  "\nIMPORTANT:\n",
  "This PCA is exploratory and does not yet define\n",
  "the final tumor plasticity axis.\n",
  "Independent stemness validation is required first.\n"
)


############################################################
############################################################
# 39. Independent stemness signature validation
#     Liu et al. 2022 — Supplementary Table S2
############################################################

library(readxl)

liu_s2 <- read_excel(
  file.path(
    project_dir,
    "data",
    "supplementary",
    "mmc4.xlsx"
  ),
  sheet = "Table S2",
  col_names = FALSE
)

# ---------------------------------------------------------
# Inspect imported structure
# ---------------------------------------------------------

print(dim(liu_s2))
print(liu_s2[1:10, 1:12])


# ---------------------------------------------------------
# Extract the 12 stemness signature names
# Row 5 contains the signature names
# ---------------------------------------------------------

stemness_signature_names <- as.character(
  unlist(liu_s2[5, 1:12])
)

print(stemness_signature_names)


# ---------------------------------------------------------
# Extract genes
# Genes begin in Row 6
# ---------------------------------------------------------

stemness_gene_sets <- lapply(
  1:12,
  function(i) {
    
    genes <- as.character(
      unlist(liu_s2[6:nrow(liu_s2), i])
    )
    
    genes <- genes[
      !is.na(genes) &
        trimws(genes) != ""
    ]
    
    unique(trimws(genes))
  }
)

names(stemness_gene_sets) <- stemness_signature_names


# ---------------------------------------------------------
# Check number of genes in each signature
# ---------------------------------------------------------

stemness_original_counts <- data.frame(
  signature = names(stemness_gene_sets),
  original_n = sapply(
    stemness_gene_sets,
    length
  )
)

print(stemness_original_counts)


# ---------------------------------------------------------
# Map each signature to genes present in GSE33532
# ---------------------------------------------------------

stemness_gene_sets_mapped <- lapply(
  stemness_gene_sets,
  function(gs) intersect(gs, tumor_genes)
)


# ---------------------------------------------------------
# Coverage table
# ---------------------------------------------------------

stemness_coverage <- data.frame(
  signature = names(stemness_gene_sets),
  
  original_n = sapply(
    stemness_gene_sets,
    length
  ),
  
  mapped_n = sapply(
    stemness_gene_sets_mapped,
    length
  )
)

stemness_coverage$missing_n <-
  stemness_coverage$original_n -
  stemness_coverage$mapped_n

stemness_coverage$coverage_pct <-
  round(
    100 *
      stemness_coverage$mapped_n /
      stemness_coverage$original_n,
    2
  )

stemness_coverage <-
  stemness_coverage[
    order(
      -stemness_coverage$coverage_pct
    ),
    ,
    drop = FALSE
  ]

print(stemness_coverage)


# ---------------------------------------------------------
# Save results
# ---------------------------------------------------------


write.csv(
  stemness_coverage,
  file.path(
    results_dir,
    "15_Liu2022_stemness_coverage.csv"
  ),
  row.names = FALSE
)

saveRDS(
  stemness_gene_sets_mapped,
  file.path(
    results_dir,
    "15_Liu2022_stemness_gene_sets_mapped.rds"
  )
)


############################################################
# 40. Independent stemness signature scoring
#     Liu et al. 2022 — 12 stemness signatures
############################################################

# ---------------------------------------------------------
# Prepare mapped Liu stemness gene sets
# ---------------------------------------------------------

liu_stemness_gene_sets <- stemness_gene_sets_mapped

# Remove signatures with fewer than 10 genes
# This is a technical safeguard, not a biological exclusion
liu_stemness_gene_sets_usable <- liu_stemness_gene_sets[
  sapply(liu_stemness_gene_sets, length) >= 10
]

cat(
  "Number of usable Liu stemness signatures:",
  length(liu_stemness_gene_sets_usable),
  "\n"
)

print(
  sapply(
    liu_stemness_gene_sets_usable,
    length
  )
)


# ---------------------------------------------------------
# Run ssGSEA
# ---------------------------------------------------------

liu_ssgsea_param <- ssgseaParam(
  expr_tumor,
  liu_stemness_gene_sets_usable
)

liu_stemness_scores <- gsva(
  liu_ssgsea_param
)

# Expected:
# rows = signatures
# columns = tumor samples

dim(liu_stemness_scores)


# ---------------------------------------------------------
# Save score matrix
# ---------------------------------------------------------

saveRDS(
  liu_stemness_scores,
  file.path(
    results_dir,
    "15_Liu2022_stemness_ssGSEA_scores.rds"
  )
)

write.csv(
  liu_stemness_scores,
  file.path(
    results_dir,
    "15_Liu2022_stemness_ssGSEA_scores.csv"
  ),
  row.names = TRUE
)


############################################################
# 41. Agreement among independent stemness signatures
############################################################

# Convert to samples × signatures
liu_stemness_scores_df <- as.data.frame(
  t(liu_stemness_scores)
)

# ---------------------------------------------------------
# Pearson correlation
# ---------------------------------------------------------

liu_stemness_cor_pearson <- cor(
  liu_stemness_scores_df,
  method = "pearson",
  use = "pairwise.complete.obs"
)

print(
  round(
    liu_stemness_cor_pearson,
    3
  )
)


# ---------------------------------------------------------
# Spearman correlation
# ---------------------------------------------------------

liu_stemness_cor_spearman <- cor(
  liu_stemness_scores_df,
  method = "spearman",
  use = "pairwise.complete.obs"
)

print(
  round(
    liu_stemness_cor_spearman,
    3
  )
)


# ---------------------------------------------------------
# Save correlation matrices
# ---------------------------------------------------------

# Save correlation matrices
# ---------------------------------------------------------

write.csv(
  liu_stemness_cor_pearson,
  file.path(
    results_dir,
    "15_Liu2022_stemness_cor_pearson.csv"
  ),
  row.names = TRUE
)

write.csv(
  liu_stemness_cor_spearman,
  file.path(
    results_dir,
    "15_Liu2022_stemness_cor_spearman.csv"
  ),
  row.names = TRUE
)
############################################################
# 42. Independent validation of Wong stemness score
############################################################

# Extract primary Wong stemness score
wong_score <- as.numeric(
  program_scores_ssgsea["Stemness_Wong", ]
)

names(wong_score) <- colnames(
  program_scores_ssgsea 
)


############################################################
# Correlation of Wong with each Liu signature
############################################################

wong_validation <- do.call(
  rbind,
  lapply(
    colnames(liu_stemness_scores_df),
    function(sig) {
      
      x <- wong_score
      y <- liu_stemness_scores_df[, sig]
      
      pearson_test <- cor.test(
        x,
        y,
        method = "pearson"
      )
      
      spearman_test <- cor.test(
        x,
        y,
        method = "spearman"
      )
      
      data.frame(
        signature = sig,
        
        pearson_r =
          unname(
            pearson_test$estimate
          ),
        
        pearson_p =
          pearson_test$p.value,
        
        spearman_rho =
          unname(
            spearman_test$estimate
          ),
        
        spearman_p =
          spearman_test$p.value
      )
    }
  )
)


############################################################
# Multiple-testing correction
############################################################

wong_validation$pearson_FDR <-
  p.adjust(
    wong_validation$pearson_p,
    method = "BH"
  )

wong_validation$spearman_FDR <-
  p.adjust(
    wong_validation$spearman_p,
    method = "BH"
  )


############################################################
# Display results
############################################################

print(
  wong_validation
)


############################################################
# Save results
############################################################

write.csv(
  wong_validation,
  file.path(
    results_dir,
    "15_Wong_independent_stemness_validation.csv"
  ),
  row.names = FALSE
)

############################################################
# 43. PCA of independent Liu stemness signatures
############################################################

cat("\nRunning PCA across the 12 Liu stemness signatures...\n")

# Samples x signatures
liu_stemness_scores_df <- as.data.frame(
  t(liu_stemness_scores)
)

# Standardize each signature
liu_stemness_scaled <- scale(
  liu_stemness_scores_df
)

# PCA
liu_stemness_pca <- prcomp(
  liu_stemness_scaled,
  center = FALSE,
  scale. = FALSE
)

cat("\nPCA variance explained:\n")
print(
  summary(liu_stemness_pca)$importance
)

# Sample-level PC scores
liu_stemness_PC_scores <- liu_stemness_pca$x

# Orient PC1 so that higher values correspond to higher Wong stemness
wong_score <- as.numeric(
  program_scores_ssgsea["Stemness_Wong", ]
)

names(wong_score) <- colnames(
  program_scores_ssgsea
)

if (
  cor(
    liu_stemness_PC_scores[, 1],
    wong_score,
    method = "pearson"
  ) < 0
) {
  liu_stemness_PC_scores[, 1] <-
    -liu_stemness_PC_scores[, 1]
  
  liu_stemness_pca$rotation[, 1] <-
    -liu_stemness_pca$rotation[, 1]
}

cat("\nCorrelation of oriented PC1 with Wong score:\n")
print(
  cor.test(
    liu_stemness_PC_scores[, 1],
    wong_score,
    method = "pearson"
  )
)

saveRDS(
  liu_stemness_pca,
  file.path(
    results_dir,
    "15_Liu_stemness_PCA.rds"
  )
)

write.csv(
  liu_stemness_PC_scores,
  file.path(
    results_dir,
    "15_Liu_stemness_PCA_scores.csv"
  ),
  row.names = TRUE
)

cat("\nPCA completed and saved.\n")
############################################################
# 44. Inspect Liu stemness PCA loadings
############################################################

pc1_loadings <- liu_stemness_pca$rotation[, 1]
pc2_loadings <- liu_stemness_pca$rotation[, 2]

pc_loadings <- data.frame(
  signature = rownames(liu_stemness_pca$rotation),
  PC1 = pc1_loadings,
  PC2 = pc2_loadings
)

pc_loadings <- pc_loadings[
  order(-abs(pc_loadings$PC1)),
]

print(pc_loadings)

############################################################
# 45. Correlation of PC1 with each Liu signature
############################################################

pc1_signature_correlations <- sapply(
  colnames(liu_stemness_scores_df),
  function(sig) {
    cor(
      liu_stemness_PC_scores[, 1],
      liu_stemness_scores_df[, sig],
      method = "pearson"
    )
  }
)

print(
  sort(
    pc1_signature_correlations,
    decreasing = TRUE
  )
)

############################################################
# 46A. EMT–Stemness relationship using corrected scores
############################################################

wong_score <- as.numeric(
  program_scores_ssgsea["Stemness_Wong", ]
)

names(wong_score) <- colnames(
  program_scores_ssgsea
)

length(wong_score)
head(wong_score)
head(names(wong_score))

###################
emt_score <- as.numeric(EMT_Byers76_score)
names(emt_score) <- names(EMT_Byers76_score)

stemness_score <- as.numeric(wong_score)
names(stemness_score) <- names(wong_score)

# Check sample alignment
stopifnot(
  identical(names(emt_score), names(stemness_score))
)

emt_stemness_pearson <- cor.test(
  emt_score,
  stemness_score,
  method = "pearson"
)

emt_stemness_spearman <- cor.test(
  emt_score,
  stemness_score,
  method = "spearman"
)

cat("\nPearson correlation: EMT vs Wong stemness\n")
print(emt_stemness_pearson)

cat("\nSpearman correlation: EMT vs Wong stemness\n")
print(emt_stemness_spearman)

############################################################
# 46B. EMT–Consensus Stemness Relationship
############################################################

# Ensure EMT score is a named vector
emt_score <- as.numeric(EMT_Byers76_score)
names(emt_score) <- names(EMT_Byers76_score)

# Extract Liu stemness PC1
pc1_score <- liu_stemness_PC_scores[, 1]
names(pc1_score) <- rownames(liu_stemness_PC_scores)

# Check sample names before correlation
cat("EMT samples:", length(emt_score), "\n")
cat("PC1 samples:", length(pc1_score), "\n")

cat("\nEMT and PC1 sample names identical:\n")
print(identical(names(emt_score), names(pc1_score)))

############################################################
# Align explicitly by sample ID
############################################################

common_samples <- intersect(
  names(emt_score),
  names(pc1_score)
)

cat("\nCommon samples:", length(common_samples), "\n")

emt_score_aligned <- emt_score[common_samples]
pc1_score_aligned <- pc1_score[common_samples]

# Final safety check
stopifnot(
  identical(
    names(emt_score_aligned),
    names(pc1_score_aligned)
  )
)

############################################################
# Pearson correlation
############################################################

emt_PC1_pearson <- cor.test(
  emt_score_aligned,
  pc1_score_aligned,
  method = "pearson"
)

############################################################
# Spearman correlation
############################################################

emt_PC1_spearman <- cor.test(
  emt_score_aligned,
  pc1_score_aligned,
  method = "spearman"
)

############################################################
# Print results
############################################################

cat("\n============================================\n")
cat("EMT vs Liu Stemness PC1\n")
cat("============================================\n")

cat("\nPearson correlation:\n")
print(emt_PC1_pearson)

cat("\nSpearman correlation:\n")
print(emt_PC1_spearman)

############################################################
# Save tumor expression matrix for downstream scripts
############################################################

saveRDS(
  expr_tumor,
  file.path(
    processed_dir,
    "GSE33532_expr_tumor.rds"
  )
)

############################################################
############################################################
# 47. Save final plasticity scores
############################################################

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n============================================\n")
cat("SECTION 47: SAVE FINAL PLASTICITY SCORES\n")
cat("============================================\n")

# ---- Byers EMT score ----
EMT_score_final <- as.numeric(EMT_Byers76_score)
names(EMT_score_final) <- colnames(expr_tumor)

stopifnot(length(EMT_score_final) == 80)
stopifnot(identical(names(EMT_score_final), colnames(expr_tumor)))


# ---- Wong stemness score ----
Wong_score_final <- as.numeric(
  program_scores_ssgsea["Stemness_Wong", ]
)

names(Wong_score_final) <- colnames(expr_tumor)

stopifnot(length(Wong_score_final) == 80)
stopifnot(identical(names(Wong_score_final), colnames(expr_tumor)))


# ---- Liu PC1 score ----
Liu_PC1_score_final <- as.numeric(
  liu_stemness_PC_scores[, 1]
)

names(Liu_PC1_score_final) <- rownames(liu_stemness_PC_scores)

stopifnot(length(Liu_PC1_score_final) == 80)
stopifnot(identical(names(Liu_PC1_score_final), colnames(expr_tumor)))


# ---- Save individual score objects ----

saveRDS(
  EMT_score_final,
  file.path(
    results_dir,
    "GSE33532_Byers76_EMT_score.rds"
  )
)

saveRDS(
  Wong_score_final,
  file.path(
    results_dir,
    "GSE33532_Wong_stemness_score.rds"
  )
)

saveRDS(
  Liu_PC1_score_final,
  file.path(
    results_dir,
    "GSE33532_Liu_PC1_stemness_score.rds"
  )
)


# ---- Combined table ----

plasticity_score_table <- data.frame(
  sample = colnames(expr_tumor),
  EMT_Byers76 = EMT_score_final,
  Wong_stemness = Wong_score_final,
  Liu_PC1_stemness = Liu_PC1_score_final,
  stringsAsFactors = FALSE
)

write.csv(
  plasticity_score_table,
  file.path(
    results_dir,
    "GSE33532_core_plasticity_scores.csv"
  ),
  row.names = FALSE
)


cat("\nSaved final plasticity scores:\n")
cat("Byers EMT:", length(EMT_score_final), "samples\n")
cat("Wong stemness:", length(Wong_score_final), "samples\n")
cat("Liu PC1 stemness:", length(Liu_PC1_score_final), "samples\n")

cat("\nSection 47 completed successfully.\n")

############################################################
# 48. Save dependencies for downstream Script 16
############################################################
# Check required objects before saving Script 16 dependencies
############################################################

required_objects <- c(
  "EMT_score",
  "wong_score",
  "Byers76_gene_order",
  "stemness_gene_sets_mapped",
  "liu_stemness_gene_sets_usable"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    envir = .GlobalEnv,
    inherits = FALSE
  )
]

if (length(missing_objects) > 0) {
  stop(
    "Missing objects required by Script 16: ",
    paste(missing_objects, collapse = ", ")
  )
}

cat("All required Script 16 objects are available.\n")


script16_dependencies <- list(
  
  EMT_score = EMT_score,
  
  wong_score = wong_score,
  
  Byers76_gene_order = Byers76_gene_order,
  
  stemness_gene_sets_mapped =
    stemness_gene_sets_mapped,
  
  liu_stemness_gene_sets_usable =
    liu_stemness_gene_sets_usable
  
)

saveRDS(
  script16_dependencies,
  file.path(
    processed_dir,
    "Script16_dependencies.rds"
  )
)

cat(
  "\nScript 16 dependencies saved successfully.\n"
)

##########################
sessionInfo()


