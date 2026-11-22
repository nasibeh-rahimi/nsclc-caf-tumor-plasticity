############################################################
# 18_GSE33532_Integrated_Analysis_Table.R
#
# Purpose:
#
# Build one analysis-ready table for the 80 GSE33532 tumor
# samples by integrating:
#
#   1. CAF-associated transcriptional program score
#      - CAF signature DISCOVERED in GSE149955
#      - transferred/scored in GSE33532 tumors
#
#   2. Byers76 EMT/plasticity score
#
#   3. Wong broad stemness score
#
#   4. Zakaria-derived putative CSC-associated state score
#
#   5. GSE33532 patient and subsample information
#
#   6. Histology
#
# IMPORTANT DATASET ROLES:
#
# GSE149955:
#   25 matched NF + 25 CAF samples
#   Used for CAF-associated differential expression.
#
# GSE33532:
#   20 patients
#   4 tumor subsamples per patient = 80 tumor samples
#   Used as the tumor analysis cohort.
#
# The CAF score in this script therefore does NOT mean that
# GSE33532 contains CAF samples. It represents the expression
# of the CAF-associated transcriptional program derived from
# GSE149955 within GSE33532 tumor samples.
############################################################

############################################################
# 1. Project paths
############################################################

project_dir <- getwd()

# CAF-derived signature was generated from GSE149955,
# but the resulting score is calculated in GSE33532.
caf_result_dir <- file.path(
  project_dir,
  "results",
  "CAF"
)

# The three saved plasticity scores are stored directly
# under results/.
score_result_dir <- file.path(
  project_dir,
  "results"
)

# Integrated analysis outputs
analysis_dir <- file.path(
  project_dir,
  "results",
  "analysis"
)

dir.create(
  analysis_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 2. Load GSE33532 CAF-associated score
############################################################

# IMPORTANT:
#
# The CAF-associated gene program was DISCOVERED in GSE149955.
# It was then TRANSFERRED to GSE33532.
#
# Therefore this object contains CAF-associated scores for
# the 80 GSE33532 tumor samples.

CAF_score <- readRDS(
  file.path(
    caf_result_dir,
    "GSE33532_CAF_score.rds"
  )
)

cat(
  "GSE33532 CAF-associated score samples:",
  length(CAF_score),
  "\n"
)

stopifnot(
  length(CAF_score) == 80
)


############################################################
# 3. Load GSE33532 tumor metadata
############################################################

# This is NOT CAF-sample metadata.
#
# It describes the 80 GSE33532 tumor samples to which the
# GSE149955-derived CAF-associated program was transferred.

tumor_meta <- read.csv(
  file.path(
    caf_result_dir,
    "GSE33532_tumor_metadata.csv"
  ),
  stringsAsFactors = FALSE
)

# Basic metadata checks

stopifnot(
  nrow(tumor_meta) == 80
)

stopifnot(
  all(
    c(
      "sample",
      "patient_id",
      "tissue_code"
    ) %in% names(tumor_meta)
  )
)

stopifnot(
  !anyDuplicated(tumor_meta$sample)
)

cat(
  "GSE33532 tumor metadata samples:",
  nrow(tumor_meta),
  "\n"
)


############################################################
# 4. Load GSE33532 plasticity scores
############################################################

# ----------------------------------------------------------
# 4A. Byers76 EMT score
# ----------------------------------------------------------

EMT_score <- readRDS(
  file.path(
    score_result_dir,
    "GSE33532_Byers76_EMT_score.rds"
  )
)


# ----------------------------------------------------------
# 4B. Wong broad stemness score
# ----------------------------------------------------------

Wong_score <- readRDS(
  file.path(
    score_result_dir,
    "GSE33532_Wong_stemness_score.rds"
  )
)


# ----------------------------------------------------------
# 4C. Zakaria-derived putative CSC-associated state score
# ----------------------------------------------------------

Zakaria_score <- readRDS(
  file.path(
    score_result_dir,
    "GSE33532_Zakaria_CSC_balanced_score.rds"
  )
)


############################################################
# 5. Check score dimensions
############################################################

cat(
  "Byers76 EMT score samples:",
  length(EMT_score),
  "\n"
)

cat(
  "Wong stemness score samples:",
  length(Wong_score),
  "\n"
)

cat(
  "Zakaria CSC-associated state score samples:",
  length(Zakaria_score),
  "\n"
)

stopifnot(
  length(EMT_score) == 80
)

stopifnot(
  length(Wong_score) == 80
)

stopifnot(
  length(Zakaria_score) == 80
)


############################################################
# 6. Check sample-name overlap across all score objects
############################################################

score_sample_sets <- list(
  CAF_associated = names(CAF_score),
  Byers76_EMT = names(EMT_score),
  Wong_stemness = names(Wong_score),
  Zakaria_CSC_state = names(Zakaria_score)
)

common_samples <- Reduce(
  intersect,
  score_sample_sets
)

cat(
  "\nNumber of samples common to all four scores:",
  length(common_samples),
  "\n"
)

stopifnot(
  length(common_samples) == 80
)


############################################################
# 7. Confirm metadata contains the same 80 samples
############################################################

metadata_samples <- tumor_meta$sample

stopifnot(
  length(intersect(
    common_samples,
    metadata_samples
  )) == 80
)

stopifnot(
  setequal(
    common_samples,
    metadata_samples
  )
)


############################################################
# 8. Create integrated score table
############################################################

analysis_table <- data.frame(
  
  sample = common_samples,
  
  CAF_associated_score =
    as.numeric(
      CAF_score[common_samples]
    ),
  
  Byers76_EMT =
    as.numeric(
      EMT_score[common_samples]
    ),
  
  Wong_stemness =
    as.numeric(
      Wong_score[common_samples]
    ),
  
  Zakaria_CSC_state =
    as.numeric(
      Zakaria_score[common_samples]
    ),
  
  stringsAsFactors = FALSE
)


############################################################
# 9. Add GSE33532 patient and subsample information
############################################################

metadata_match <- match(
  analysis_table$sample,
  tumor_meta$sample
)

stopifnot(
  !any(is.na(metadata_match))
)

analysis_table$patient_id <-
  tumor_meta$patient_id[metadata_match]

analysis_table$subsample <-
  tumor_meta$tissue_code[metadata_match]


############################################################
# 10. Add histology
############################################################

# Histology is assigned at the patient level.
#
# Each patient contributes four tumor subsamples:
# A, B, C, and D.
#
# Therefore all four tumor subsamples from a given patient
# receive the same histology.

histology_map <- c(
  
  "02" = "Adenocarcinoma",
  "07" = "Adenocarcinoma",
  "08" = "Mixed",
  "10" = "Mixed",
  "11" = "Mixed",
  "17" = "Adenocarcinoma",
  "22" = "Adenocarcinoma",
  "27" = "Adenocarcinoma",
  "28" = "Adenocarcinoma",
  
  "34" = "Squamous",
  
  "41" = "Adenocarcinoma",
  "42" = "Squamous",
  "43" = "Mixed",
  "57" = "Adenocarcinoma",
  "58" = "Adenocarcinoma",
  "59" = "Mixed",
  "62" = "Adenocarcinoma",
  "85" = "Squamous",
  "93" = "Squamous",
  "98" = "Mixed"
)

analysis_table$histology <- unname(
  histology_map[
    formatC(
      as.numeric(analysis_table$patient_id),
      width = 2,
      format = "d",
      flag = "0"
    )
  ]
)


############################################################
# 11. Verify histology assignment
############################################################

stopifnot(
  !anyNA(analysis_table$histology)
)

# Every patient must have only one histology.
histology_check <- tapply(
  analysis_table$histology,
  analysis_table$patient_id,
  function(x) length(unique(x))
)

stopifnot(
  all(histology_check == 1)
)


############################################################
# 12. Convert grouping variables to factors
############################################################

analysis_table$patient_id <-
  factor(
    analysis_table$patient_id
  )

analysis_table$subsample <-
  factor(
    analysis_table$subsample,
    levels = c(
      "A",
      "B",
      "C",
      "D"
    )
  )

analysis_table$histology <-
  factor(
    analysis_table$histology,
    levels = c(
      "Adenocarcinoma",
      "Mixed",
      "Squamous"
    )
  )


############################################################
# 13. Validate histology distribution
############################################################

cat("\n============================================\n")
cat("HISTOLOGY DISTRIBUTION\n")
cat("============================================\n")

print(
  table(
    analysis_table$histology
  )
)


############################################################
# 14. Validate patient structure
############################################################

cat("\n============================================\n")
cat("PATIENT DISTRIBUTION\n")
cat("============================================\n")

patient_counts <- table(
  analysis_table$patient_id
)

print(patient_counts)

# 20 patients, four tumor subsamples each
stopifnot(
  length(patient_counts) == 20
)

stopifnot(
  all(patient_counts == 4)
)


############################################################
# 15. Validate subsample structure
############################################################

cat("\n============================================\n")
cat("SUBSAMPLE DISTRIBUTION\n")
cat("============================================\n")

print(
  table(
    analysis_table$subsample
  )
)

stopifnot(
  all(
    table(
      analysis_table$subsample
    ) == 20
  )
)


############################################################
# 16. Check missing values
############################################################

cat("\n============================================\n")
cat("MISSING VALUES\n")
cat("============================================\n")

missing_values <- colSums(
  is.na(analysis_table)
)

print(missing_values)

stopifnot(
  !anyNA(analysis_table)
)


############################################################
# 17. Check score distributions
############################################################

cat("\n============================================\n")
cat("SCORE SUMMARIES\n")
cat("============================================\n")

print(
  summary(
    analysis_table[
      ,
      c(
        "CAF_associated_score",
        "Byers76_EMT",
        "Wong_stemness",
        "Zakaria_CSC_state"
      )
    ]
  )
)


############################################################
# 18. Display first rows
############################################################

cat("\n============================================\n")
cat("FIRST 10 ROWS\n")
cat("============================================\n")

print(
  head(
    analysis_table,
    10
  )
)


############################################################
# 19. Save integrated analysis table
############################################################

write.csv(
  analysis_table,
  file.path(
    analysis_dir,
    "GSE33532_integrated_analysis_table.csv"
  ),
  row.names = FALSE
)

saveRDS(
  analysis_table,
  file.path(
    analysis_dir,
    "GSE33532_integrated_analysis_table.rds"
  )
)


############################################################
# 20. Save session information
############################################################

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    analysis_dir,
    "18_integrated_analysis_table_sessionInfo.txt"
  )
)


############################################################
# 21. Final confirmation
############################################################

cat("\n============================================\n")
cat("INTEGRATED ANALYSIS TABLE CREATED\n")
cat("============================================\n")

cat(
  "Samples:",
  nrow(analysis_table),
  "\n"
)

cat(
  "Patients:",
  length(
    unique(
      analysis_table$patient_id
    )
  ),
  "\n"
)

cat(
  "Tumor subsamples per patient: 4\n"
)

cat(
  "CAF-associated score: GSE149955-derived signature\n"
)

cat(
  "Applied to: GSE33532 tumors\n"
)

cat(
  "Byers76 EMT: loaded\n"
)

cat(
  "Wong stemness: loaded\n"
)

cat(
  "Zakaria CSC-associated state: loaded\n"
)

cat(
  "\nFiles saved to:\n",
  analysis_dir,
  "\n"
)

cat(
  "\nScript 18 completed successfully.\n"
)