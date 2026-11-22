# ============================================================
# Script 24
# GSE33532 Integrated Statistical Results
#
# Purpose:
# Integrate the primary and secondary CAF–plasticity analyses
# into one reproducible statistical results table.
#
# Core statistics are recalculated directly from the saved
# patient-level analysis table.
#
# Sensitivity analyses are loaded from Scripts 21 and 23.
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 0. Project setup
# ------------------------------------------------------------

source("scripts/00_project_setup.R")

library(dplyr)

project_dir <- getwd()

analysis_dir <- file.path(
  project_dir,
  "results",
  "analysis"
)

# ------------------------------------------------------------
# 1. Load patient-level analysis table
# ------------------------------------------------------------

patient_table_file <- file.path(
  analysis_dir,
  "GSE33532_patient_level_analysis_table.rds"
)

stopifnot(file.exists(patient_table_file))

patient_data <- readRDS(
  patient_table_file
)

# ------------------------------------------------------------
# 2. Validate patient-level dataset
# ------------------------------------------------------------

required_columns <- c(
  "patient_id",
  "histology",
  "CAF_associated_score",
  "Byers76_EMT",
  "Wong_stemness"
)

missing_columns <- setdiff(
  required_columns,
  colnames(patient_data)
)

if (length(missing_columns) > 0) {
  
  stop(
    paste(
      "Missing required columns:",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}

# Expected independent patient sample size
stopifnot(nrow(patient_data) == 20)

# Each patient must occur once
stopifnot(
  length(unique(patient_data$patient_id)) ==
    nrow(patient_data)
)

# ------------------------------------------------------------
# 3. Basic cohort description
# ------------------------------------------------------------

n_patients <- nrow(patient_data)

histology_counts <- table(
  patient_data$histology
)

histology_text <- paste(
  names(histology_counts),
  histology_counts,
  sep = " = ",
  collapse = "; "
)

cat("\n============================================================\n")
cat("COHORT DESCRIPTION\n")
cat("============================================================\n\n")

cat(
  "Independent patients:",
  n_patients,
  "\n"
)

cat(
  "Histology distribution:",
  histology_text,
  "\n"
)

# ------------------------------------------------------------
# 4. Primary analysis: CAF vs Byers76 EMT
# ------------------------------------------------------------

emt_pearson_test <- cor.test(
  patient_data$CAF_associated_score,
  patient_data$Byers76_EMT,
  method = "pearson"
)

emt_spearman_test <- cor.test(
  patient_data$CAF_associated_score,
  patient_data$Byers76_EMT,
  method = "spearman",
  exact = FALSE
)

emt_model <- lm(
  Byers76_EMT ~ CAF_associated_score,
  data = patient_data
)

emt_adjusted_model <- lm(
  Byers76_EMT ~ CAF_associated_score + histology,
  data = patient_data
)

# ------------------------------------------------------------
# 5. Secondary analysis: CAF vs Wong stemness
# ------------------------------------------------------------

wong_pearson_test <- cor.test(
  patient_data$CAF_associated_score,
  patient_data$Wong_stemness,
  method = "pearson"
)

wong_spearman_test <- cor.test(
  patient_data$CAF_associated_score,
  patient_data$Wong_stemness,
  method = "spearman",
  exact = FALSE
)

wong_model <- lm(
  Wong_stemness ~ CAF_associated_score,
  data = patient_data
)

wong_adjusted_model <- lm(
  Wong_stemness ~ CAF_associated_score + histology,
  data = patient_data
)

# ------------------------------------------------------------
# 6. Helper function for regression coefficient extraction
# ------------------------------------------------------------

extract_model_coefficient <- function(
    model,
    term
) {
  
  coef_table <- summary(model)$coefficients
  
  estimate <- coef_table[
    term,
    "Estimate"
  ]
  
  se <- coef_table[
    term,
    "Std. Error"
  ]
  
  t_value <- coef_table[
    term,
    "t value"
  ]
  
  p_value <- coef_table[
    term,
    "Pr(>|t|)"
  ]
  
  df <- df.residual(model)
  
  critical_t <- qt(
    0.975,
    df = df
  )
  
  ci_low <- estimate -
    critical_t * se
  
  ci_high <- estimate +
    critical_t * se
  
  data.frame(
    estimate = estimate,
    SE = se,
    t = t_value,
    p = p_value,
    CI_low = ci_low,
    CI_high = ci_high
  )
}

# ------------------------------------------------------------
# 7. Extract regression results
# ------------------------------------------------------------

emt_unadjusted <- extract_model_coefficient(
  emt_model,
  "CAF_associated_score"
)

emt_adjusted <- extract_model_coefficient(
  emt_adjusted_model,
  "CAF_associated_score"
)

wong_unadjusted <- extract_model_coefficient(
  wong_model,
  "CAF_associated_score"
)

wong_adjusted <- extract_model_coefficient(
  wong_adjusted_model,
  "CAF_associated_score"
)

# ------------------------------------------------------------
# 8. Model summaries
# ------------------------------------------------------------

emt_unadjusted_summary <- summary(
  emt_model
)

emt_adjusted_summary <- summary(
  emt_adjusted_model
)

wong_unadjusted_summary <- summary(
  wong_model
)

wong_adjusted_summary <- summary(
  wong_adjusted_model
)

# ------------------------------------------------------------
# 9. Load sensitivity-analysis outputs
# ------------------------------------------------------------

emt_loo_file <- file.path(
  analysis_dir,
  "GSE33532_CAF_EMT_leave_one_out.csv"
)

emt_influence_file <- file.path(
  analysis_dir,
  "GSE33532_CAF_EMT_influence_diagnostics.csv"
)

wong_loo_file <- file.path(
  analysis_dir,
  "GSE33532_CAF_Wong_leave_one_out.csv"
)

wong_influence_file <- file.path(
  analysis_dir,
  "GSE33532_CAF_Wong_influence_diagnostics.csv"
)

sensitivity_files <- c(
  emt_loo_file,
  emt_influence_file,
  wong_loo_file,
  wong_influence_file
)

missing_sensitivity <- sensitivity_files[
  !file.exists(sensitivity_files)
]

if (length(missing_sensitivity) > 0) {
  
  stop(
    paste(
      "Missing sensitivity-analysis files:\n",
      paste(
        missing_sensitivity,
        collapse = "\n"
      )
    )
  )
}

# Load sensitivity results
emt_loo <- read.csv(
  emt_loo_file,
  stringsAsFactors = FALSE
)

emt_influence <- read.csv(
  emt_influence_file,
  stringsAsFactors = FALSE
)

wong_loo <- read.csv(
  wong_loo_file,
  stringsAsFactors = FALSE
)

wong_influence <- read.csv(
  wong_influence_file,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 10. Influence summaries
# ------------------------------------------------------------
# EMT influence summary
# -----------------------------------------------------

emt_influence <- read.csv(
  emt_influence_file,
  stringsAsFactors = FALSE
)

emt_flagged <- emt_influence %>%
  dplyr::filter(Any_diagnostic_flag)

emt_influence_summary <- data.frame(
  relationship = "CAF–Byers EMT",
  n_influence_checked = nrow(emt_influence),
  n_flagged = nrow(emt_flagged),
  flagged_patients = if (nrow(emt_flagged) > 0) {
    paste(emt_flagged$patient_id, collapse = "; ")
  } else {
    NA_character_
  },
  cook_flag = if (any(emt_influence$Cook_flag)) {
    paste(
      emt_influence$patient_id[emt_influence$Cook_flag],
      collapse = "; "
    )
  } else {
    NA_character_
  },
  leverage_flag = if (any(emt_influence$Leverage_flag)) {
    paste(
      emt_influence$patient_id[emt_influence$Leverage_flag],
      collapse = "; "
    )
  } else {
    NA_character_
  },
  residual_flag = if (any(emt_influence$Residual_flag)) {
    paste(
      emt_influence$patient_id[emt_influence$Residual_flag],
      collapse = "; "
    )
  } else {
    NA_character_
  },
  DFBETA_flag = NA_character_,
  stringsAsFactors = FALSE
)


# -----------------------------------------------------
# Wong influence summary
# -----------------------------------------------------

wong_influence <- read.csv(
  wong_influence_file,
  stringsAsFactors = FALSE
)

wong_influence$Any_diagnostic_flag <-
  wong_influence$flag_cook |
  wong_influence$flag_leverage |
  wong_influence$flag_residual |
  wong_influence$flag_DFBETA

wong_flagged <- wong_influence %>%
  dplyr::filter(Any_diagnostic_flag)

wong_influence_summary <- data.frame(
  relationship = "CAF–Wong",
  n_influence_checked = nrow(wong_influence),
  n_flagged = nrow(wong_flagged),
  flagged_patients = if (nrow(wong_flagged) > 0) {
    paste(wong_flagged$patient_id, collapse = "; ")
  } else {
    NA_character_
  },
  cook_flag = if (any(wong_influence$flag_cook)) {
    paste(
      wong_influence$patient_id[wong_influence$flag_cook],
      collapse = "; "
    )
  } else {
    NA_character_
  },
  leverage_flag = if (any(wong_influence$flag_leverage)) {
    paste(
      wong_influence$patient_id[wong_influence$flag_leverage],
      collapse = "; "
    )
  } else {
    NA_character_
  },
  residual_flag = if (any(wong_influence$flag_residual)) {
    paste(
      wong_influence$patient_id[wong_influence$flag_residual],
      collapse = "; "
    )
  } else {
    NA_character_
  },
  DFBETA_flag = if (any(wong_influence$flag_DFBETA)) {
    paste(
      wong_influence$patient_id[wong_influence$flag_DFBETA],
      collapse = "; "
    )
  } else {
    NA_character_
  },
  stringsAsFactors = FALSE
)


# -----------------------------------------------------
# Combine
# -----------------------------------------------------

influence_summary <- rbind(
  emt_influence_summary,
  wong_influence_summary
)


# -----------------------------------------------------
# Save
# -----------------------------------------------------

write.csv(
  influence_summary,
  file.path(
    analysis_dir,
    "GSE33532_integrated_influence_summary.csv"
  ),
  row.names = FALSE
)

print(influence_summary)


# ------------------------------------------------------------
# 11. Leave-one-out ranges
# ------------------------------------------------------------

emt_loo_min <- min(
  emt_loo$Pearson_r,
  na.rm = TRUE
)

emt_loo_max <- max(
  emt_loo$Pearson_r,
  na.rm = TRUE
)

wong_loo_min <- min(
  wong_loo$Pearson_r,
  na.rm = TRUE
)

wong_loo_max <- max(
  wong_loo$Pearson_r,
  na.rm = TRUE
)

# ------------------------------------------------------------
# 12. Influential patient labels
# ------------------------------------------------------------

emt_influential_text <- if (
  nrow(emt_flagged) == 0
) {
  
  "None"
  
} else {
  
  paste(
    emt_flagged$patient_id,
    collapse = ", "
  )
}

wong_influential_text <- if (
  nrow(wong_flagged) == 0
) {
  
  "None"
  
} else {
  
  paste(
    wong_flagged$patient_id,
    collapse = ", "
  )
}

# ------------------------------------------------------------
# 13. Integrated primary + secondary results
# ------------------------------------------------------------

results_table <- data.frame(
  
  analysis = c(
    "CAF–Byers76 EMT",
    "CAF–Byers76 EMT",
    "CAF–Byers76 EMT",
    "CAF–Wong stemness",
    "CAF–Wong stemness",
    "CAF–Wong stemness"
  ),
  
  model = c(
    "Pearson correlation",
    "Spearman correlation",
    "Unadjusted linear regression",
    "Pearson correlation",
    "Spearman correlation",
    "Unadjusted linear regression"
  ),
  
  sample_size = rep(
    n_patients,
    6
  ),
  
  estimate = c(
    
    unname(
      emt_pearson_test$estimate
    ),
    
    unname(
      emt_spearman_test$estimate
    ),
    
    emt_unadjusted$estimate,
    
    unname(
      wong_pearson_test$estimate
    ),
    
    unname(
      wong_spearman_test$estimate
    ),
    
    wong_unadjusted$estimate
  ),
  
  estimate_type = c(
    "Pearson r",
    "Spearman rho",
    "CAF beta",
    "Pearson r",
    "Spearman rho",
    "CAF beta"
  ),
  
  p_value = c(
    
    emt_pearson_test$p.value,
    
    emt_spearman_test$p.value,
    
    emt_unadjusted$p,
    
    wong_pearson_test$p.value,
    
    wong_spearman_test$p.value,
    
    wong_unadjusted$p
  ),
  
  CI_low = c(
    
    emt_pearson_test$conf.int[1],
    
    NA,
    
    emt_unadjusted$CI_low,
    
    wong_pearson_test$conf.int[1],
    
    NA,
    
    wong_unadjusted$CI_low
  ),
  
  CI_high = c(
    
    emt_pearson_test$conf.int[2],
    
    NA,
    
    emt_unadjusted$CI_high,
    
    wong_pearson_test$conf.int[2],
    
    NA,
    
    wong_unadjusted$CI_high
  ),
  
  R_squared = c(
    
    NA,
    NA,
    emt_unadjusted_summary$r.squared,
    
    NA,
    NA,
    wong_unadjusted_summary$r.squared
  ),
  
  adjusted_R_squared = c(
    
    NA,
    NA,
    emt_unadjusted_summary$adj.r.squared,
    
    NA,
    NA,
    wong_unadjusted_summary$adj.r.squared
  ),
  
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 14. Histology-adjusted results
# ------------------------------------------------------------

adjusted_table <- data.frame(
  
  analysis = c(
    "CAF–Byers76 EMT",
    "CAF–Wong stemness"
  ),
  
  model = c(
    "Histology-adjusted linear regression",
    "Histology-adjusted linear regression"
  ),
  
  sample_size = c(
    n_patients,
    n_patients
  ),
  
  predictor = c(
    "CAF_associated_score",
    "CAF_associated_score"
  ),
  
  estimate = c(
    emt_adjusted$estimate,
    wong_adjusted$estimate
  ),
  
  SE = c(
    emt_adjusted$SE,
    wong_adjusted$SE
  ),
  
  t = c(
    emt_adjusted$t,
    wong_adjusted$t
  ),
  
  p_value = c(
    emt_adjusted$p,
    wong_adjusted$p
  ),
  
  CI_low = c(
    emt_adjusted$CI_low,
    wong_adjusted$CI_low
  ),
  
  CI_high = c(
    emt_adjusted$CI_high,
    wong_adjusted$CI_high
  ),
  
  R_squared = c(
    emt_adjusted_summary$r.squared,
    wong_adjusted_summary$r.squared
  ),
  
  adjusted_R_squared = c(
    emt_adjusted_summary$adj.r.squared,
    wong_adjusted_summary$adj.r.squared
  ),
  
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 15. Sensitivity results
# ------------------------------------------------------------

sensitivity_table <- data.frame(
  
  analysis = c(
    "CAF–Byers76 EMT",
    "CAF–Wong stemness"
  ),
  
  sample_size = c(
    n_patients,
    n_patients
  ),
  
  full_Pearson_r = c(
    
    unname(
      emt_pearson_test$estimate
    ),
    
    unname(
      wong_pearson_test$estimate
    )
  ),
  
  full_Pearson_p = c(
    
    emt_pearson_test$p.value,
    
    wong_pearson_test$p.value
  ),
  
  LOO_Pearson_min = c(
    emt_loo_min,
    wong_loo_min
  ),
  
  LOO_Pearson_max = c(
    emt_loo_max,
    wong_loo_max
  ),
  
  LOO_direction_consistent = c(
    
    emt_loo_min < 0 &
      emt_loo_max < 0,
    
    wong_loo_min < 0 &
      wong_loo_max < 0
  ),
  
  n_influential_flagged = c(
    nrow(emt_flagged),
    nrow(wong_flagged)
  ),
  
  influential_patients = c(
    emt_influential_text,
    wong_influential_text
  ),
  
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 16. Print results
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("INTEGRATED STATISTICAL RESULTS\n")
cat("============================================================\n\n")

cat(
  "Independent patients:",
  n_patients,
  "\n"
)

cat(
  "Histology:",
  histology_text,
  "\n\n"
)

cat("------------------------------------------------------------\n")
cat("PRIMARY + SECONDARY UNADJUSTED ANALYSES\n")
cat("------------------------------------------------------------\n\n")

print(
  results_table,
  row.names = FALSE
)

cat("\n------------------------------------------------------------\n")
cat("HISTOLOGY-ADJUSTED MODELS\n")
cat("------------------------------------------------------------\n\n")

print(
  adjusted_table,
  row.names = FALSE
)

cat("\n------------------------------------------------------------\n")
cat("SENSITIVITY ANALYSES\n")
cat("------------------------------------------------------------\n\n")

print(
  sensitivity_table,
  row.names = FALSE
)

# ------------------------------------------------------------
# 17. Interpretation table
# ------------------------------------------------------------

interpretation_table <- data.frame(
  
  analysis = c(
    "CAF–Byers76 EMT",
    "CAF–Wong stemness"
  ),
  
  primary_interpretation = c(
    
    paste0(
      "No evidence of a patient-level association ",
      "between the CAF-associated transcriptional ",
      "score and the Byers76 EMT score."
    ),
    
    paste0(
      "Negative direction of association with the ",
      "Wong broad stemness score, but not statistically ",
      "supported."
    )
  ),
  
  histology_adjusted_interpretation = c(
    
    paste0(
      "The CAF coefficient remained close to zero ",
      "after adjustment for histology."
    ),
    
    paste0(
      "The CAF coefficient remained negative after ",
      "adjustment for histology, but was not statistically ",
      "supported."
    )
  ),
  
  sensitivity_interpretation = c(
    
    paste0(
      "LOO Pearson correlations remained negative ",
      "and close to zero; the result was not dependent ",
      "on a single patient."
    ),
    
    paste0(
      "LOO Pearson correlations remained negative ",
      "after omission of every individual patient, ",
      "although influential observations were detected."
    )
  ),
  
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 18. Study limitations
# ------------------------------------------------------------

limitations <- c(
  
  "1. The patient-level analysis contains 20 independent patients.",
  
  "2. GSE33532 contains four tumor subsamples per patient; these were aggregated before patient-level inferential analysis.",
  
  "3. The cohort contains adenocarcinoma, mixed, and squamous histologies, introducing biological heterogeneity.",
  
  "4. The CAF score represents a CAF-associated transcriptional program transferred from GSE149955; it is not a direct measurement of CAF abundance.",
  
  "5. The Byers76 score represents an EMT-associated transcriptional state and should not be interpreted as a direct measurement of all forms of tumor plasticity.",
  
  "6. The Wong score represents a broad stemness/pluripotency-like transcriptional program and should not be interpreted as a direct measurement of CSC abundance.",
  
  "7. These analyses are observational associations and cannot establish causality.",
  
  "8. Bulk transcriptomic measurements mix signals from multiple cellular compartments and may obscure tumor-cell-specific relationships.",
  
  "9. The small patient-level sample size limits precision and statistical power.",
  
  "10. Influence diagnostics identified potentially influential observations in the CAF–Wong analysis; these observations were retained rather than removed from the primary analysis."
)

# ------------------------------------------------------------
# 19. Save output tables
# ------------------------------------------------------------

results_table_file <- file.path(
  analysis_dir,
  "GSE33532_integrated_statistical_results.csv"
)

adjusted_table_file <- file.path(
  analysis_dir,
  "GSE33532_integrated_histology_adjusted_results.csv"
)

sensitivity_table_file <- file.path(
  analysis_dir,
  "GSE33532_integrated_sensitivity_results.csv"
)

interpretation_file <- file.path(
  analysis_dir,
  "GSE33532_integrated_interpretation.csv"
)

write.csv(
  results_table,
  results_table_file,
  row.names = FALSE
)

write.csv(
  adjusted_table,
  adjusted_table_file,
  row.names = FALSE
)

write.csv(
  sensitivity_table,
  sensitivity_table_file,
  row.names = FALSE
)

write.csv(
  interpretation_table,
  interpretation_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 20. Human-readable report
# ------------------------------------------------------------

report_file <- file.path(
  analysis_dir,
  "GSE33532_integrated_statistical_results.txt"
)

report_lines <- c(
  
  "GSE33532 INTEGRATED STATISTICAL RESULTS",
  "=======================================",
  "",
  
  "STUDY COHORT",
  "------------",
  
  paste(
    "Independent patients:",
    n_patients
  ),
  
  paste(
    "Histology distribution:",
    histology_text
  ),
  
  "",
  
  "PRIMARY ANALYSIS: CAF-associated score vs Byers76 EMT",
  "-----------------------------------------------------",
  
  paste(
    "Pearson r =",
    round(
      unname(emt_pearson_test$estimate),
      4
    ),
    "; p =",
    signif(
      emt_pearson_test$p.value,
      4
    )
  ),
  
  paste(
    "Pearson 95% CI = [",
    round(
      emt_pearson_test$conf.int[1],
      4
    ),
    ",",
    round(
      emt_pearson_test$conf.int[2],
      4
    ),
    "]"
  ),
  
  paste(
    "Spearman rho =",
    round(
      unname(emt_spearman_test$estimate),
      4
    ),
    "; p =",
    signif(
      emt_spearman_test$p.value,
      4
    )
  ),
  
  paste(
    "Unadjusted CAF beta =",
    round(
      emt_unadjusted$estimate,
      4
    ),
    "; 95% CI [",
    round(
      emt_unadjusted$CI_low,
      4
    ),
    ",",
    round(
      emt_unadjusted$CI_high,
      4
    ),
    "]; p =",
    signif(
      emt_unadjusted$p,
      4
    )
  ),
  
  paste(
    "Unadjusted R-squared =",
    round(
      emt_unadjusted_summary$r.squared,
      4
    )
  ),
  
  paste(
    "Adjusted R-squared =",
    round(
      emt_unadjusted_summary$adj.r.squared,
      4
    )
  ),
  
  "",
  
  "Histology-adjusted CAF beta:",
  
  paste(
    round(
      emt_adjusted$estimate,
      4
    ),
    "; 95% CI [",
    round(
      emt_adjusted$CI_low,
      4
    ),
    ",",
    round(
      emt_adjusted$CI_high,
      4
    ),
    "]; p =",
    signif(
      emt_adjusted$p,
      4
    )
  ),
  
  "",
  
  "CAF–EMT SENSITIVITY",
  
  paste(
    "LOO Pearson range:",
    round(
      emt_loo_min,
      4
    ),
    "to",
    round(
      emt_loo_max,
      4
    )
  ),
  
  paste(
    "Flagged influential observations:",
    emt_influential_text
  ),
  
  "",
  
  "SECONDARY ANALYSIS: CAF-associated score vs Wong stemness",
  "---------------------------------------------------------",
  
  paste(
    "Pearson r =",
    round(
      unname(wong_pearson_test$estimate),
      4
    ),
    "; p =",
    signif(
      wong_pearson_test$p.value,
      4
    )
  ),
  
  paste(
    "Pearson 95% CI = [",
    round(
      wong_pearson_test$conf.int[1],
      4
    ),
    ",",
    round(
      wong_pearson_test$conf.int[2],
      4
    ),
    "]"
  ),
  
  paste(
    "Spearman rho =",
    round(
      unname(wong_spearman_test$estimate),
      4
    ),
    "; p =",
    signif(
      wong_spearman_test$p.value,
      4
    )
  ),
  
  paste(
    "Unadjusted CAF beta =",
    round(
      wong_unadjusted$estimate,
      4
    ),
    "; 95% CI [",
    round(
      wong_unadjusted$CI_low,
      4
    ),
    ",",
    round(
      wong_unadjusted$CI_high,
      4
    ),
    "]; p =",
    signif(
      wong_unadjusted$p,
      4
    )
  ),
  
  paste(
    "Unadjusted R-squared =",
    round(
      wong_unadjusted_summary$r.squared,
      4
    )
  ),
  
  paste(
    "Adjusted R-squared =",
    round(
      wong_unadjusted_summary$adj.r.squared,
      4
    )
  ),
  
  "",
  
  "Histology-adjusted CAF beta:",
  
  paste(
    round(
      wong_adjusted$estimate,
      4
    ),
    "; 95% CI [",
    round(
      wong_adjusted$CI_low,
      4
    ),
    ",",
    round(
      wong_adjusted$CI_high,
      4
    ),
    "]; p =",
    signif(
      wong_adjusted$p,
      4
    )
  ),
  
  "",
  
  "CAF–WONG SENSITIVITY",
  
  paste(
    "LOO Pearson range:",
    round(
      wong_loo_min,
      4
    ),
    "to",
    round(
      wong_loo_max,
      4
    )
  ),
  
  paste(
    "Flagged influential observations:",
    wong_influential_text
  ),
  
  "",
  
  "STUDY LIMITATIONS",
  "-----------------",
  
  limitations,
  
  "",
  
  "OVERALL STATISTICAL INTERPRETATION",
  "-----------------------------------",
  
  "The patient-level analysis did not provide statistical",
  "support for the primary CAF–Byers76 EMT association.",
  
  "The CAF–Wong analysis showed a negative direction but",
  "was not statistically supported.",
  
  "The Wong relationship remained negative under leave-one-out",
  "analysis, although influence diagnostics identified",
  "potentially influential observations.",
  
  "These analyses are observational and should not be interpreted",
  "as evidence of causality."
)

writeLines(
  report_lines,
  report_file
)

# ------------------------------------------------------------
# 21. Final verification
# ------------------------------------------------------------

output_files <- c(
  
  results_table_file,
  
  adjusted_table_file,
  
  sensitivity_table_file,
  
  interpretation_file,
  
  report_file
)

cat("\n============================================================\n")
cat("SCRIPT 24 OUTPUT VERIFICATION\n")
cat("============================================================\n\n")

for (f in output_files) {
  
  cat(
    ifelse(
      file.exists(f),
      "OK: ",
      "MISSING: "
    ),
    f,
    "\n"
  )
}

cat("\n============================================================\n")
cat("Script 24 completed successfully.\n")
cat("============================================================\n")

