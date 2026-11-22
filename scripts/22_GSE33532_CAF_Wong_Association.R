# ============================================================
# Script 22
# GSE33532: CAF-associated activity and Wong stemness
#
# Secondary analysis:
# Is CAF-associated transcriptional activity associated
# with the broad Wong stemness/pluripotency-like program?
#
# Unit of analysis:
# Patient level
#
# n = 20 patients
#
# Primary analysis:
#   Wong_stemness ~ CAF_associated_score
#
# Sensitivity analysis:
#   Wong_stemness ~ CAF_associated_score + histology
#
# This is a secondary analysis.
# ============================================================


# ============================================================
# 1. Setup
# ============================================================

rm(list = ls())

source("scripts/00_project_setup.R")

project_dir <- getwd()

analysis_dir <- file.path(
  project_dir,
  "results",
  "analysis"
)

figure_dir <- file.path(
  project_dir,
  "figures"
)

dir.create(
  analysis_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. Load patient-level analysis table
# ============================================================

patient_file <- file.path(
  analysis_dir,
  "GSE33532_patient_level_analysis_table.csv"
)

patient_table <- read.csv(
  patient_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# 3. Validate input
# ============================================================

required_columns <- c(
  "patient_id",
  "histology",
  "CAF_associated_score",
  "Wong_stemness"
)

missing_columns <- setdiff(
  required_columns,
  colnames(patient_table)
)

if (length(missing_columns) > 0) {
  stop(
    paste(
      "Missing required columns:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

if (nrow(patient_table) != 20) {
  stop(
    paste(
      "Expected 20 patients, found",
      nrow(patient_table)
    )
  )
}

if (anyDuplicated(patient_table$patient_id) > 0) {
  stop("Duplicate patient IDs detected.")
}

if (
  anyNA(patient_table$CAF_associated_score) ||
  anyNA(patient_table$Wong_stemness) ||
  anyNA(patient_table$histology)
) {
  stop("Missing values detected in analysis variables.")
}


# ============================================================
# 4. Define variables
# ============================================================

x <- patient_table$CAF_associated_score
y <- patient_table$Wong_stemness


# ============================================================
# 5. Histology distribution
# ============================================================

histology_counts <- table(
  patient_table$histology
)

cat("\n============================================\n")
cat("PATIENT-LEVEL COHORT\n")
cat("============================================\n\n")

cat(
  "Number of patients:",
  nrow(patient_table),
  "\n\n"
)

print(histology_counts)


# ============================================================
# 6. Descriptive summaries
# ============================================================

cat("\n============================================\n")
cat("VARIABLE SUMMARIES\n")
cat("============================================\n\n")

cat("CAF-associated score:\n")
print(summary(x))

cat("\nWong stemness score:\n")
print(summary(y))


# ============================================================
# 7. Pearson correlation
# ============================================================

pearson_test <- cor.test(
  x,
  y,
  method = "pearson"
)


# ============================================================
# 8. Spearman correlation
# ============================================================

spearman_test <- cor.test(
  x,
  y,
  method = "spearman",
  exact = FALSE
)


# ============================================================
# 9. Unadjusted linear model
# ============================================================

unadjusted_model <- lm(
  Wong_stemness ~ CAF_associated_score,
  data = patient_table
)

unadjusted_summary <- summary(
  unadjusted_model
)


# ============================================================
# 10. Histology-adjusted sensitivity model
# ============================================================

adjusted_model <- lm(
  Wong_stemness ~ CAF_associated_score + histology,
  data = patient_table
)

adjusted_summary <- summary(
  adjusted_model
)


# ============================================================
# 11. Print correlation results
# ============================================================

cat("\n============================================\n")
cat("PEARSON CORRELATION\n")
cat("============================================\n\n")

cat(
  "Pearson r = ",
  round(
    unname(pearson_test$estimate),
    6
  ),
  "\n"
)

cat(
  "95% CI = [",
  round(
    pearson_test$conf.int[1],
    6
  ),
  ", ",
  round(
    pearson_test$conf.int[2],
    6
  ),
  "]\n",
  sep = ""
)

cat(
  "p-value = ",
  signif(
    pearson_test$p.value,
    6
  ),
  "\n"
)


cat("\n============================================\n")
cat("SPEARMAN CORRELATION\n")
cat("============================================\n\n")

cat(
  "Spearman rho = ",
  round(
    unname(spearman_test$estimate),
    6
  ),
  "\n"
)

cat(
  "p-value = ",
  signif(
    spearman_test$p.value,
    6
  ),
  "\n"
)


# ============================================================
# 12. Print unadjusted model
# ============================================================

cat("\n============================================\n")
cat("UNADJUSTED LINEAR MODEL\n")
cat("============================================\n\n")

print(unadjusted_summary)


# ============================================================
# 13. Print histology-adjusted model
# ============================================================

cat("\n============================================\n")
cat("HISTOLOGY-ADJUSTED MODEL\n")
cat("============================================\n\n")

print(adjusted_summary)


# ============================================================
# 14. Extract model coefficients
# ============================================================

unadjusted_coef <- coef(
  unadjusted_summary
)

adjusted_coef <- coef(
  adjusted_summary
)


# ============================================================
# 15. Save correlation results
# ============================================================

correlation_results <- data.frame(
  
  analysis = c(
    "Pearson",
    "Spearman"
  ),
  
  estimate = c(
    unname(pearson_test$estimate),
    unname(spearman_test$estimate)
  ),
  
  p_value = c(
    pearson_test$p.value,
    spearman_test$p.value
  ),
  
  stringsAsFactors = FALSE
)

write.csv(
  correlation_results,
  file.path(
    analysis_dir,
    "GSE33532_CAF_Wong_correlations.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. Save model summaries
# ============================================================

capture.output(
  unadjusted_summary,
  file = file.path(
    analysis_dir,
    "GSE33532_CAF_Wong_unadjusted_model.txt"
  )
)

capture.output(
  adjusted_summary,
  file = file.path(
    analysis_dir,
    "GSE33532_CAF_Wong_histology_adjusted_model.txt"
  )
)


# ============================================================
# 17. Scatterplot: CAF vs Wong
# ============================================================

png(
  filename = file.path(
    figure_dir,
    "GSE33532_CAF_vs_Wong_stemness.png"
  ),
  width = 1800,
  height = 1400,
  res = 180
)

plot(
  x,
  y,
  pch = 19,
  xlab = "CAF-associated transcriptional score",
  ylab = "Wong broad stemness score",
  main = "CAF-associated Activity vs Wong Stemness"
)

abline(
  unadjusted_model,
  lty = 2
)

legend(
  "topleft",
  legend = paste(
    "Pearson r =",
    round(
      unname(pearson_test$estimate),
      3
    )
  ),
  bty = "n"
)

dev.off()


# ============================================================
# 18. Scatterplot with histology labels
# ============================================================

histology_symbols <- c(
  "Adenocarcinoma" = 16,
  "Mixed" = 17,
  "Squamous" = 15
)

point_symbols <- histology_symbols[
  as.character(patient_table$histology)
]

png(
  filename = file.path(
    figure_dir,
    "GSE33532_CAF_vs_Wong_stemness_histology.png"
  ),
  width = 1800,
  height = 1400,
  res = 180
)

plot(
  x,
  y,
  pch = point_symbols,
  xlab = "CAF-associated transcriptional score",
  ylab = "Wong broad stemness score",
  main = "CAF-associated Activity vs Wong Stemness by Histology"
)

abline(
  unadjusted_model,
  lty = 2
)

legend(
  "topleft",
  legend = names(histology_symbols),
  pch = histology_symbols,
  bty = "n"
)

dev.off()


# ============================================================
# 19. Patient-level analysis table with residuals
# ============================================================

patient_results <- patient_table

patient_results$Wong_fitted_unadjusted <- fitted(
  unadjusted_model
)

patient_results$Wong_residual_unadjusted <- residuals(
  unadjusted_model
)

patient_results$Wong_fitted_adjusted <- fitted(
  adjusted_model
)

patient_results$Wong_residual_adjusted <- residuals(
  adjusted_model
)

write.csv(
  patient_results,
  file.path(
    analysis_dir,
    "GSE33532_CAF_Wong_patient_level_results.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 20. Save concise summary
# ============================================================

summary_file <- file.path(
  analysis_dir,
  "GSE33532_CAF_Wong_analysis_summary.txt"
)

sink(summary_file)

cat(
  "GSE33532 CAF-associated Activity and Wong Stemness\n"
)

cat(
  "==================================================\n\n"
)

cat(
  "Analysis unit: patient level\n"
)

cat(
  "Number of patients:",
  nrow(patient_table),
  "\n\n"
)

cat(
  "Histology distribution:\n"
)

print(histology_counts)

cat("\n")

cat(
  "Pearson r:",
  round(
    unname(pearson_test$estimate),
    6
  ),
  "\n"
)

cat(
  "Pearson 95% CI:",
  paste(
    round(
      pearson_test$conf.int,
      6
    ),
    collapse = " to "
  ),
  "\n"
)

cat(
  "Pearson p:",
  signif(
    pearson_test$p.value,
    6
  ),
  "\n\n"
)

cat(
  "Spearman rho:",
  round(
    unname(spearman_test$estimate),
    6
  ),
  "\n"
)

cat(
  "Spearman p:",
  signif(
    spearman_test$p.value,
    6
  ),
  "\n\n"
)

cat(
  "Unadjusted CAF coefficient:",
  unadjusted_coef[
    "CAF_associated_score",
    "Estimate"
  ],
  "\n"
)

cat(
  "Unadjusted CAF coefficient p:",
  unadjusted_coef[
    "CAF_associated_score",
    "Pr(>|t|)"
  ],
  "\n"
  
)

cat(
  "Unadjusted R-squared:",
  unadjusted_summary$r.squared,
  "\n\n"
)

cat(
  "Histology-adjusted CAF coefficient:",
  adjusted_coef[
    "CAF_associated_score",
    "Estimate"
  ],
  "\n"
)

cat(
  "Histology-adjusted CAF coefficient p:",
  adjusted_coef[
    "CAF_associated_score",
    "Pr(>|t|)"
  ],
  "\n"
)

cat(
  "Histology-adjusted R-squared:",
  adjusted_summary$r.squared,
  "\n"
)

sink()


# ============================================================
# 21. Final verification
# ============================================================

cat("\n============================================\n")
cat("SCRIPT 22 COMPLETE\n")
cat("============================================\n\n")

cat(
  "Correlation results:\n",
  file.path(
    analysis_dir,
    "GSE33532_CAF_Wong_correlations.csv"
  ),
  "\n\n"
)

cat(
  "Unadjusted model:\n",
  file.path(
    analysis_dir,
    "GSE33532_CAF_Wong_unadjusted_model.txt"
  ),
  "\n\n"
)

cat(
  "Histology-adjusted model:\n",
  file.path(
    analysis_dir,
    "GSE33532_CAF_Wong_histology_adjusted_model.txt"
  ),
  "\n\n"
)

cat(
  "Patient-level results:\n",
  file.path(
    analysis_dir,
    "GSE33532_CAF_Wong_patient_level_results.csv"
  ),
  "\n\n"
)

cat(
  "Summary:\n",
  summary_file,
  "\n\n"
)

cat("Figures:\n")

cat(
  file.path(
    figure_dir,
    "GSE33532_CAF_vs_Wong_stemness.png"
  ),
  "\n"
)

cat(
  file.path(
    figure_dir,
    "GSE33532_CAF_vs_Wong_stemness_histology.png"
  ),
  "\n"
)

