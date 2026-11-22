############################################################
# 19_GSE33532_Patient_Level_Aggregation.R
#
# Purpose:
#   Aggregate the four tumor subsamples (A/B/C/D) from each
#   GSE33532 patient into one patient-level observation.
#
# Input:
#   GSE33532_integrated_analysis_table.csv
#
# Output:
#   Patient-level analysis table with n = 20 patients
#
# Primary downstream use:
#   CAF-associated score vs Byers76 EMT
#
# Secondary:
#   CAF-associated score vs Wong stemness
#
# Exploratory:
#   CAF-associated score vs Zakaria CSC-associated state
############################################################

############################################################
# 1. Setup
############################################################

rm(list = ls())

project_dir <- getwd()

analysis_dir <- file.path(
  project_dir,
  "results",
  "analysis"
)

output_csv <- file.path(
  analysis_dir,
  "GSE33532_patient_level_analysis_table.csv"
)

output_rds <- file.path(
  analysis_dir,
  "GSE33532_patient_level_analysis_table.rds"
)


############################################################
# 2. Load integrated 80-sample table
############################################################

integrated_file <- file.path(
  analysis_dir,
  "GSE33532_integrated_analysis_table.csv"
)

stopifnot(file.exists(integrated_file))

analysis_table <- read.csv(
  integrated_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Integrated table loaded.\n")
cat("Samples:", nrow(analysis_table), "\n")


############################################################
# 3. Basic structural validation
############################################################

required_columns <- c(
  "sample",
  "patient_id",
  "subsample",
  "histology",
  "CAF_associated_score",
  "Byers76_EMT",
  "Wong_stemness",
  "Zakaria_CSC_state"
)

stopifnot(
  all(required_columns %in% colnames(analysis_table))
)

stopifnot(
  nrow(analysis_table) == 80
)

stopifnot(
  length(unique(analysis_table$patient_id)) == 20
)

stopifnot(
  all(table(analysis_table$patient_id) == 4)
)

stopifnot(
  setequal(
    unique(analysis_table$subsample),
    c("A", "B", "C", "D")
  )
)

stopifnot(
  !anyNA(analysis_table$histology)
)


############################################################
# 4. Check that each patient has one histology
############################################################

histology_check <- aggregate(
  histology ~ patient_id,
  data = analysis_table,
  FUN = function(x) length(unique(x))
)

stopifnot(
  all(histology_check$histology == 1)
)

cat("Histology is consistent within patients.\n")


############################################################
# 5. Check score missingness before aggregation
############################################################

score_columns <- c(
  "CAF_associated_score",
  "Byers76_EMT",
  "Wong_stemness",
  "Zakaria_CSC_state"
)

missing_counts <- sapply(
  analysis_table[, score_columns, drop = FALSE],
  function(x) sum(is.na(x))
)

print(missing_counts)

stopifnot(
  all(missing_counts == 0)
)


############################################################
# 6. Aggregate scores to patient level
############################################################

patient_level <- aggregate(
  analysis_table[, score_columns],
  by = list(
    patient_id = analysis_table$patient_id
  ),
  FUN = mean
)


############################################################
# 7. Add patient-level histology
############################################################

patient_histology <- unique(
  analysis_table[, c("patient_id", "histology")]
)

stopifnot(
  nrow(patient_histology) == 20
)

patient_level <- merge(
  patient_level,
  patient_histology,
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE
)


############################################################
# 8. Reorder columns
############################################################

patient_level <- patient_level[
  ,
  c(
    "patient_id",
    "histology",
    "CAF_associated_score",
    "Byers76_EMT",
    "Wong_stemness",
    "Zakaria_CSC_state"
  )
]


############################################################
# 9. Sort by patient ID
############################################################

patient_level <- patient_level[
  order(patient_level$patient_id),
  ,
  drop = FALSE
]

rownames(patient_level) <- NULL


############################################################
# 10. Convert histology to factor
############################################################

patient_level$histology <- factor(
  patient_level$histology,
  levels = c(
    "Adenocarcinoma",
    "Mixed",
    "Squamous"
  )
)


############################################################
# 11. Final validation
############################################################

cat("\n============================================\n")
cat("PATIENT-LEVEL VALIDATION\n")
cat("============================================\n")

cat(
  "Patients:",
  nrow(patient_level),
  "\n"
)

cat(
  "Histology groups:\n"
)

print(
  table(patient_level$histology)
)

cat(
  "\nMissing values:\n"
)

print(
  sapply(
    patient_level,
    function(x) sum(is.na(x))
  )
)

stopifnot(
  nrow(patient_level) == 20
)

stopifnot(
  !anyNA(patient_level)
)

stopifnot(
  length(unique(patient_level$patient_id)) == 20
)


############################################################
# 12. Display patient-level table
############################################################

cat("\n============================================\n")
cat("PATIENT-LEVEL ANALYSIS TABLE\n")
cat("============================================\n")

print(patient_level)


############################################################
# 13. Score summaries
############################################################

cat("\n============================================\n")
cat("PATIENT-LEVEL SCORE SUMMARIES\n")
cat("============================================\n")

print(
  summary(
    patient_level[
      ,
      score_columns,
      drop = FALSE
    ]
  )
)


############################################################
# 14. Save patient-level table
############################################################

write.csv(
  patient_level,
  output_csv,
  row.names = FALSE
)

saveRDS(
  patient_level,
  output_rds
)


############################################################
# 15. Final confirmation
############################################################

cat("\n============================================\n")
cat("SCRIPT 19 COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

cat(
  "Patient-level observations:",
  nrow(patient_level),
  "\n"
)

cat(
  "Output CSV:\n",
  output_csv,
  "\n"
)

cat(
  "Output RDS:\n",
  output_rds,
  "\n"
)



