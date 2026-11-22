############################################################
# Script 26
# GSE33532: Zakaria-derived putative CSC-associated state
# Exploratory analysis
#
# Purpose:
#   Evaluate the relationship between a Zakaria-derived
#   putative CSC-associated transcriptional program and
#   independent tumor-state scores in GSE33532.
#
# Primary inferential unit:
#   Patient-level averages across tumor subsamples A-D.
#
# Important:
#   The Zakaria score represents a transcriptional state,
#   not direct CSC abundance.
############################################################

# 1. Packages
# ============================================================
library(dplyr)
library(tidyr)
library(ggplot2)


# 2. Project directories
# ============================================================
project_dir <- getwd()

results_dir <- file.path(
  project_dir,
  "results",
  "analysis"
)

figure_dir <- file.path(
  project_dir,
  "figures"
)

# ============================================================
# 3. Load existing scores and Zakaria gene program
# ============================================================

Zakaria_program <- readRDS(
  file.path(
    project_dir,
    "results",
    "GSE33532_Zakaria_gene_program.rds"
  )
)

Zakaria_score <- readRDS(
  file.path(
    project_dir,
    "results",
    "GSE33532_Zakaria_CSC_balanced_score.rds"
  )
)

EMT_score <- readRDS(
  file.path(
    project_dir,
    "results",
    "GSE33532_Byers76_EMT_score.rds"
  )
)

Wong_score <- readRDS(
  file.path(
    project_dir,
    "results",
    "GSE33532_Wong_stemness_score.rds"
  )
)

CAF_score <- readRDS(
  file.path(
    project_dir,
    "results",
    "CAF",
    "GSE33532_CAF_score.rds"
  )
)

tumor_meta <- read.csv(
  file.path(
    project_dir,
    "results",
    "CAF",
    "GSE33532_tumor_metadata.csv"
  ),
  stringsAsFactors = FALSE
)

# ============================================================
# 4. Verify loaded objects
# ============================================================

cat("\n============================================\n")
cat("ZAKARIA EXPLORATORY ANALYSIS\n")
cat("============================================\n\n")

cat("Zakaria program object:\n")
print(class(Zakaria_program))
print(dim(Zakaria_program))

cat("\nZakaria score:\n")
print(length(Zakaria_score))

cat("\nByers EMT score:\n")
print(length(EMT_score))

cat("\nWong stemness score:\n")
print(length(Wong_score))

cat("\nCAF score:\n")
print(length(CAF_score))

cat("\nTumor metadata:\n")
print(dim(tumor_meta))
print(names(tumor_meta))

# ============================================================
# 5. Inspect Zakaria gene program structure
# ============================================================

cat("\n============================================\n")
cat("ZAKARIA GENE PROGRAM STRUCTURE\n")
cat("============================================\n\n")

cat("Zakaria program names:\n")
print(names(Zakaria_program))

cat("\nStructure:\n")
str(Zakaria_program, max.level = 1)

# ============================================================
# 6. Verify sample-name alignment
# ============================================================

cat("\n============================================\n")
cat("SAMPLE ALIGNMENT CHECK\n")
cat("============================================\n\n")

sample_names <- names(Zakaria_score)

cat("Number of Zakaria samples:", length(sample_names), "\n")

cat("\nFirst 10 Zakaria sample names:\n")
print(head(sample_names, 10))

cat("\nSample names in EMT score identical to Zakaria:\n")
print(identical(names(EMT_score), sample_names))

cat("\nSample names in Wong score identical to Zakaria:\n")
print(identical(names(Wong_score), sample_names))

cat("\nSample names in CAF score identical to Zakaria:\n")
print(identical(names(CAF_score), sample_names))

cat("\nSample names in metadata identical to Zakaria:\n")
print(identical(tumor_meta$sample, sample_names))

# ============================================================
# ============================================================
# 7. Build sample-level integrated table
# ============================================================

patient_histology <- c(
  "02" = "adenocarcinoma",
  "07" = "adenocarcinoma",
  "08" = "mixed",
  "10" = "mixed",
  "11" = "mixed",
  "17" = "adenocarcinoma",
  "22" = "adenocarcinoma",
  "27" = "adenocarcinoma",
  "28" = "adenocarcinoma",
  "34" = "squamous",
  "41" = "adenocarcinoma",
  "42" = "squamous",
  "43" = "mixed",
  "57" = "adenocarcinoma",
  "58" = "adenocarcinoma",
  "59" = "mixed",
  "62" = "adenocarcinoma",
  "85" = "squamous",
  "93" = "squamous",
  "98" = "mixed"
)


integrated_sample <- tumor_meta %>%
  mutate(
    CAF_associated_score = as.numeric(CAF_score),
    Byers76_EMT = as.numeric(EMT_score),
    Wong_stemness = as.numeric(Wong_score),
    Zakaria_CSC_state = as.numeric(Zakaria_score),
    histology = unname(
      patient_histology[
        sprintf("%02d", patient_id)
      ]
    )
  )


cat("\n============================================\n")
cat("SAMPLE-LEVEL INTEGRATED TABLE\n")
cat("============================================\n\n")

print(dim(integrated_sample))

cat("\nHistology counts:\n")
print(table(integrated_sample$histology))

cat("\nMissing values:\n")
print(colSums(is.na(integrated_sample)))

# ============================================================
# 8. Aggregate tumor subsamples to patient level
# ============================================================

patient_level <- integrated_sample %>%
  group_by(patient_id) %>%
  summarise(
    CAF_associated_score = mean(
      CAF_associated_score,
      na.rm = TRUE
    ),
    Byers76_EMT = mean(
      Byers76_EMT,
      na.rm = TRUE
    ),
    Wong_stemness = mean(
      Wong_stemness,
      na.rm = TRUE
    ),
    Zakaria_CSC_state = mean(
      Zakaria_CSC_state,
      na.rm = TRUE
    ),
    histology = first(histology),
    n_subsamples = n(),
    .groups = "drop"
  )


cat("\n============================================\n")
cat("PATIENT-LEVEL TABLE\n")
cat("============================================\n\n")

print(dim(patient_level))

cat("\nNumber of patients:\n")
print(nrow(patient_level))

cat("\nSubsamples per patient:\n")
print(table(patient_level$n_subsamples))

cat("\nPatient-level histology:\n")
print(table(patient_level$histology))

cat("\nFirst rows:\n")
print(head(patient_level))

# ============================================================
# 9. Exploratory patient-level correlations
# ============================================================

correlation_results <- data.frame(
  comparison = character(),
  method = character(),
  n = integer(),
  estimate = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)


run_correlation <- function(
    x,
    y,
    x_name,
    y_name
) {
  
  pearson <- cor.test(
    x,
    y,
    method = "pearson"
  )
  
  spearman <- cor.test(
    x,
    y,
    method = "spearman",
    exact = FALSE
  )
  
  rbind(
    data.frame(
      comparison = paste(
        x_name,
        "vs",
        y_name
      ),
      method = "Pearson",
      n = sum(
        complete.cases(x, y)
      ),
      estimate = unname(
        pearson$estimate
      ),
      p_value = pearson$p.value
    ),
    
    data.frame(
      comparison = paste(
        x_name,
        "vs",
        y_name
      ),
      method = "Spearman",
      n = sum(
        complete.cases(x, y)
      ),
      estimate = unname(
        spearman$estimate
      ),
      p_value = spearman$p.value
    )
  )
}


correlation_results <- rbind(
  
  run_correlation(
    patient_level$Zakaria_CSC_state,
    patient_level$Byers76_EMT,
    "Zakaria_CSC_state",
    "Byers76_EMT"
  ),
  
  run_correlation(
    patient_level$Zakaria_CSC_state,
    patient_level$Wong_stemness,
    "Zakaria_CSC_state",
    "Wong_stemness"
  ),
  
  run_correlation(
    patient_level$Zakaria_CSC_state,
    patient_level$CAF_associated_score,
    "Zakaria_CSC_state",
    "CAF_associated_score"
  )
)


# BH correction across the six correlation tests
correlation_results$BH_adjusted_p <- p.adjust(
  correlation_results$p_value,
  method = "BH"
)


cat("\n============================================\n")
cat("EXPLORATORY CORRELATION RESULTS\n")
cat("============================================\n\n")

print(
  correlation_results,
  row.names = FALSE
)

# ============================================================
# 10. Save exploratory Zakaria results
# ============================================================

write.csv(
  correlation_results,
  file.path(
    results_dir,
    "GSE33532_Zakaria_exploratory_correlations.csv"
  ),
  row.names = FALSE
)


write.csv(
  patient_level,
  file.path(
    results_dir,
    "GSE33532_Zakaria_patient_level_analysis.csv"
  ),
  row.names = FALSE
)


write.csv(
  integrated_sample,
  file.path(
    results_dir,
    "GSE33532_Zakaria_sample_level_analysis.csv"
  ),
  row.names = FALSE
)


# Save a concise interpretation table

Zakaria_interpretation <- data.frame(
  analysis = c(
    "Zakaria vs Byers EMT",
    "Zakaria vs Wong stemness",
    "Zakaria vs CAF-associated program"
  ),
  interpretation = c(
    "Positive but non-significant association at the patient level.",
    "Negative but non-significant association at the patient level.",
    "Weak negative but non-significant association at the patient level."
  ),
  n_patients = 20,
  stringsAsFactors = FALSE
)


write.csv(
  Zakaria_interpretation,
  file.path(
    results_dir,
    "GSE33532_Zakaria_exploratory_interpretation.csv"
  ),
  row.names = FALSE
)


cat("\n============================================\n")
cat("OUTPUTS SAVED\n")
cat("============================================\n\n")

cat(
  "Correlation results:\n",
  file.path(
    results_dir,
    "GSE33532_Zakaria_exploratory_correlations.csv"
  ),
  "\n\n"
)

cat(
  "Patient-level analysis:\n",
  file.path(
    results_dir,
    "GSE33532_Zakaria_patient_level_analysis.csv"
  ),
  "\n\n"
)

cat(
  "Sample-level analysis:\n",
  file.path(
    results_dir,
    "GSE33532_Zakaria_sample_level_analysis.csv"
  ),
  "\n\n"
)

cat(
  "Interpretation:\n",
  file.path(
    results_dir,
    "GSE33532_Zakaria_exploratory_interpretation.csv"
  ),
  "\n"
)





