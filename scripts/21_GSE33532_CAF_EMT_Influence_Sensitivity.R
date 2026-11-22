# ============================================================
# Script 21
# GSE33532: Influence and Sensitivity Analysis
#
# Purpose:
# Assess whether the patient-level CAF–EMT association
# is strongly driven by individual patients.
#
# Primary model:
#   Byers76_EMT ~ CAF_associated_score
#
# This script does NOT remove or modify any observations.
# It performs diagnostic and leave-one-out sensitivity analyses.
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

dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


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
  "Byers76_EMT"
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
  anyNA(patient_table$Byers76_EMT)
) {
  stop("Missing values detected in primary variables.")
}


# ============================================================
# 4. Define variables
# ============================================================

x <- patient_table$CAF_associated_score
y <- patient_table$Byers76_EMT

patient_labels <- as.character(patient_table$patient_id)


# ============================================================
# 5. Fit primary unadjusted linear model
# ============================================================

primary_model <- lm(
  Byers76_EMT ~ CAF_associated_score,
  data = patient_table
)

model_summary <- summary(primary_model)

cat("\n============================================\n")
cat("PRIMARY MODEL\n")
cat("============================================\n\n")

print(model_summary)


# ============================================================
# 6. Influence diagnostics
# ============================================================

cooks_d <- cooks.distance(primary_model)

leverage_values <- hatvalues(primary_model)

standardized_residuals <- rstandard(primary_model)

studentized_residuals <- rstudent(primary_model)

dfbetas_values <- dfbetas(primary_model)


# ============================================================
# 7. Create diagnostic table
# ============================================================

diagnostic_table <- data.frame(
  patient_id = patient_table$patient_id,
  histology = patient_table$histology,
  CAF_associated_score = x,
  Byers76_EMT = y,
  Cook_distance = as.numeric(cooks_d),
  leverage = as.numeric(leverage_values),
  standardized_residual = as.numeric(
    standardized_residuals
  ),
  studentized_residual = as.numeric(
    studentized_residuals
  ),
  DFBETA_CAF = as.numeric(
    dfbetas_values[, "CAF_associated_score"]
  ),
  stringsAsFactors = FALSE
)

diagnostic_table <- diagnostic_table[
  order(
    diagnostic_table$Cook_distance,
    decreasing = TRUE
  ),
]


# ============================================================
# 8. Diagnostic thresholds
#
# These are heuristic reference points, not significance tests.
# ============================================================

n <- nrow(patient_table)

p <- length(coef(primary_model))

cook_threshold <- 4 / n

leverage_threshold <- 2 * p / n

standardized_residual_threshold <- 2


cat("\n============================================\n")
cat("DIAGNOSTIC REFERENCE VALUES\n")
cat("============================================\n\n")

cat(
  "Cook's distance reference threshold (4/n): ",
  round(cook_threshold, 4),
  "\n"
)

cat(
  "Leverage reference threshold (2p/n): ",
  round(leverage_threshold, 4),
  "\n"
)

cat(
  "Absolute standardized residual reference: ",
  standardized_residual_threshold,
  "\n\n"
)


# ============================================================
# 9. Flag diagnostic observations
# ============================================================

diagnostic_table$Cook_flag <-
  diagnostic_table$Cook_distance > cook_threshold

diagnostic_table$Leverage_flag <-
  diagnostic_table$leverage > leverage_threshold

diagnostic_table$Residual_flag <-
  abs(
    diagnostic_table$standardized_residual
  ) > standardized_residual_threshold

diagnostic_table$Any_diagnostic_flag <-
  diagnostic_table$Cook_flag |
  diagnostic_table$Leverage_flag |
  diagnostic_table$Residual_flag


# ============================================================
# 10. Print most influential observations
# ============================================================

cat("\n============================================\n")
cat("MOST INFLUENTIAL PATIENTS\n")
cat("============================================\n\n")

print(
  diagnostic_table[
    ,
    c(
      "patient_id",
      "histology",
      "CAF_associated_score",
      "Byers76_EMT",
      "Cook_distance",
      "leverage",
      "standardized_residual",
      "Cook_flag",
      "Leverage_flag",
      "Residual_flag"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 11. Leave-one-out sensitivity analysis
# ============================================================

loo_results <- vector(
  "list",
  length = nrow(patient_table)
)

for (i in seq_len(nrow(patient_table))) {
  
  reduced_table <- patient_table[-i, ]
  
  pearson_test <- cor.test(
    reduced_table$CAF_associated_score,
    reduced_table$Byers76_EMT,
    method = "pearson"
  )
  
  spearman_test <- cor.test(
    reduced_table$CAF_associated_score,
    reduced_table$Byers76_EMT,
    method = "spearman",
    exact = FALSE
  )
  
  loo_model <- lm(
    Byers76_EMT ~ CAF_associated_score,
    data = reduced_table
  )
  
  loo_results[[i]] <- data.frame(
    omitted_patient = patient_table$patient_id[i],
    omitted_histology = patient_table$histology[i],
    
    Pearson_r = unname(
      pearson_test$estimate
    ),
    
    Pearson_p = pearson_test$p.value,
    
    Spearman_rho = unname(
      spearman_test$estimate
    ),
    
    Spearman_p = spearman_test$p.value,
    
    slope = coef(
      loo_model
    )[["CAF_associated_score"]],
    
    stringsAsFactors = FALSE
  )
}

loo_table <- do.call(
  rbind,
  loo_results
)


# ============================================================
# 12. Full-cohort correlations for comparison
# ============================================================

full_pearson <- cor.test(
  x,
  y,
  method = "pearson"
)

full_spearman <- cor.test(
  x,
  y,
  method = "spearman",
  exact = FALSE
)


# ============================================================
# 13. Print leave-one-out results
# ============================================================

cat("\n============================================\n")
cat("FULL COHORT CORRELATIONS\n")
cat("============================================\n\n")

cat(
  "Pearson r = ",
  round(
    unname(full_pearson$estimate),
    4
  ),
  "\n"
)

cat(
  "Pearson p = ",
  signif(
    full_pearson$p.value,
    4
  ),
  "\n\n"
)

cat(
  "Spearman rho = ",
  round(
    unname(full_spearman$estimate),
    4
  ),
  "\n"
)

cat(
  "Spearman p = ",
  signif(
    full_spearman$p.value,
    4
  ),
  "\n\n"
)


cat("\n============================================\n")
cat("LEAVE-ONE-OUT CORRELATION RANGE\n")
cat("============================================\n\n")

cat(
  "Pearson range: ",
  round(
    min(loo_table$Pearson_r),
    4
  ),
  " to ",
  round(
    max(loo_table$Pearson_r),
    4
  ),
  "\n"
)

cat(
  "Spearman range: ",
  round(
    min(loo_table$Spearman_rho),
    4
  ),
  " to ",
  round(
    max(loo_table$Spearman_rho),
    4
  ),
  "\n\n"
)


# ============================================================
# 14. Identify observations with strongest LOO impact
# ============================================================

full_r <- unname(
  full_pearson$estimate
)

full_rho <- unname(
  full_spearman$estimate
)

loo_table$Pearson_change <- (
  loo_table$Pearson_r - full_r
)

loo_table$Spearman_change <- (
  loo_table$Spearman_rho - full_rho
)

loo_table$abs_Pearson_change <- abs(
  loo_table$Pearson_change
)

loo_table$abs_Spearman_change <- abs(
  loo_table$Spearman_change
)

loo_table <- loo_table[
  order(
    loo_table$abs_Pearson_change,
    decreasing = TRUE
  ),
]


cat("\n============================================\n")
cat("LARGEST PEARSON LOO CHANGES\n")
cat("============================================\n\n")

print(
  loo_table[
    ,
    c(
      "omitted_patient",
      "omitted_histology",
      "Pearson_r",
      "Pearson_change"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 15. Cook's distance plot
# ============================================================

png(
  filename = file.path(
    figure_dir,
    "GSE33532_CAF_EMT_Cooks_distance.png"
  ),
  width = 1800,
  height = 1200,
  res = 180
)

plot(
  cooks_d,
  type = "h",
  main = "Cook's Distance: CAF-associated Score → EMT",
  xlab = "Patient index",
  ylab = "Cook's distance"
)

abline(
  h = cook_threshold,
  lty = 2
)

text(
  x = seq_along(cooks_d),
  y = cooks_d,
  labels = patient_labels,
  pos = 3,
  cex = 0.7
)

dev.off()


# ============================================================
# 16. Leave-one-out Pearson correlation plot
# ============================================================

loo_plot_order <- loo_table[
  order(
    as.numeric(
      as.character(
        loo_table$omitted_patient
      )
    )
  ),
]

png(
  filename = file.path(
    figure_dir,
    "GSE33532_CAF_EMT_leave_one_out_Pearson.png"
  ),
  width = 1800,
  height = 1200,
  res = 180
)

plot(
  seq_len(nrow(loo_plot_order)),
  loo_plot_order$Pearson_r,
  type = "b",
  pch = 19,
  xaxt = "n",
  main = "Leave-One-Out Pearson Correlation",
  xlab = "Omitted patient",
  ylab = "Pearson r"
)

axis(
  1,
  at = seq_len(nrow(loo_plot_order)),
  labels = loo_plot_order$omitted_patient,
  las = 2,
  cex.axis = 0.7
)

abline(
  h = full_r,
  lty = 2
)

abline(
  h = 0,
  lty = 3
)

dev.off()


# ============================================================
# 17. Save results
# ============================================================

write.csv(
  diagnostic_table,
  file.path(
    analysis_dir,
    "GSE33532_CAF_EMT_influence_diagnostics.csv"
  ),
  row.names = FALSE
)

write.csv(
  loo_table,
  file.path(
    analysis_dir,
    "GSE33532_CAF_EMT_leave_one_out.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. Save concise text summary
# ============================================================

summary_file <- file.path(
  analysis_dir,
  "GSE33532_CAF_EMT_influence_summary.txt"
)

sink(summary_file)

cat("GSE33532 CAF–EMT Influence and Sensitivity Analysis\n")
cat("====================================================\n\n")

cat("Primary model:\n")
cat("Byers76_EMT ~ CAF_associated_score\n\n")

cat("Number of patients:", n, "\n\n")

cat(
  "Full-cohort Pearson r:",
  round(full_r, 6),
  "\n"
)

cat(
  "Full-cohort Pearson p:",
  signif(full_pearson$p.value, 6),
  "\n\n"
)

cat(
  "Full-cohort Spearman rho:",
  round(full_rho, 6),
  "\n"
)

cat(
  "Full-cohort Spearman p:",
  signif(full_spearman$p.value, 6),
  "\n\n"
)

cat(
  "Cook's distance reference threshold (4/n):",
  round(cook_threshold, 6),
  "\n"
)

cat(
  "Leverage reference threshold (2p/n):",
  round(leverage_threshold, 6),
  "\n\n"
)

cat("Patients exceeding Cook's distance reference:\n")

cook_patients <- diagnostic_table$patient_id[
  diagnostic_table$Cook_flag
]

if (length(cook_patients) == 0) {
  cat("None\n")
} else {
  cat(
    paste(
      cook_patients,
      collapse = ", "
    ),
    "\n"
  )
}

cat("\nPatients exceeding leverage reference:\n")

leverage_patients <- diagnostic_table$patient_id[
  diagnostic_table$Leverage_flag
]

if (length(leverage_patients) == 0) {
  cat("None\n")
} else {
  cat(
    paste(
      leverage_patients,
      collapse = ", "
    ),
    "\n"
  )
}

cat("\nPatients exceeding standardized residual reference:\n")

residual_patients <- diagnostic_table$patient_id[
  diagnostic_table$Residual_flag
]

if (length(residual_patients) == 0) {
  cat("None\n")
} else {
  cat(
    paste(
      residual_patients,
      collapse = ", "
    ),
    "\n"
  )
}

cat("\nLeave-one-out Pearson range:\n")
cat(
  round(
    min(loo_table$Pearson_r),
    6
  ),
  " to ",
  round(
    max(loo_table$Pearson_r),
    6
  ),
  "\n"
)

cat("\nLeave-one-out Spearman range:\n")
cat(
  round(
    min(loo_table$Spearman_rho),
    6
  ),
  " to ",
  round(
    max(loo_table$Spearman_rho),
    6
  ),
  "\n"
)

sink()


# ============================================================
# 19. Final verification
# ============================================================

cat("\n============================================\n")
cat("SCRIPT 21 COMPLETE\n")
cat("============================================\n\n")

cat(
  "Diagnostic table:\n",
  file.path(
    analysis_dir,
    "GSE33532_CAF_EMT_influence_diagnostics.csv"
  ),
  "\n\n"
)

cat(
  "Leave-one-out table:\n",
  file.path(
    analysis_dir,
    "GSE33532_CAF_EMT_leave_one_out.csv"
  ),
  "\n\n"
)

cat(
  "Summary:\n",
  summary_file,
  "\n\n"
)

cat("Figures saved:\n")
cat(
  file.path(
    figure_dir,
    "GSE33532_CAF_EMT_Cooks_distance.png"
  ),
  "\n"
)

cat(
  file.path(
    figure_dir,
    "GSE33532_CAF_EMT_leave_one_out_Pearson.png"
  ),
  "\n"
)