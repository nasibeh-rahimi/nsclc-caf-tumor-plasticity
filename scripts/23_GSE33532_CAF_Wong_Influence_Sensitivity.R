# ============================================================
# Script 23
# GSE33532 CAF–Wong Influence and Leave-One-Out Sensitivity
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 0. Project setup
# ------------------------------------------------------------
source("scripts/00_project_setup.R")


library(ggplot2)
library(dplyr)

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

# ------------------------------------------------------------
# 1. Load patient-level analysis table
# ------------------------------------------------------------

patient_file <- file.path(
  analysis_dir,
  "GSE33532_patient_level_analysis_table.rds"
)

stopifnot(file.exists(patient_file))

patient_data <- readRDS(patient_file)

# ------------------------------------------------------------
# 2. Validate required variables
# ------------------------------------------------------------

required_columns <- c(
  "patient_id",
  "histology",
  "CAF_associated_score",
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
      paste(missing_columns, collapse = ", ")
    )
  )
}

# Remove incomplete observations only if necessary
analysis_data <- patient_data %>%
  filter(
    !is.na(CAF_associated_score),
    !is.na(Wong_stemness)
  )

# Expected patient-level sample size
stopifnot(nrow(analysis_data) == 20)

# Each patient must appear once
stopifnot(
  length(unique(analysis_data$patient_id)) ==
    nrow(analysis_data)
)

# ------------------------------------------------------------
# 3. Define variables
# ------------------------------------------------------------

x <- analysis_data$CAF_associated_score
y <- analysis_data$Wong_stemness

# ------------------------------------------------------------
# 4. Fit primary unadjusted model
# ------------------------------------------------------------

model <- lm(
  Wong_stemness ~ CAF_associated_score,
  data = analysis_data
)

model_summary <- summary(model)

cat("\n============================================================\n")
cat("CAF–Wong primary model\n")
cat("============================================================\n\n")

print(model_summary)

# ------------------------------------------------------------
# 5. Influence diagnostics
# ------------------------------------------------------------

influence <- influence.measures(model)

cooks_distance <- cooks.distance(model)

leverage <- hatvalues(model)

standardized_residuals <- rstandard(model)

studentized_residuals <- rstudent(model)

dfbetas_matrix <- dfbetas(model)

# CAF coefficient DFBETA
dfbeta_caf <- dfbetas_matrix[
  ,
  "CAF_associated_score"
]

# ------------------------------------------------------------
# 6. Heuristic influence thresholds
# ------------------------------------------------------------

n <- nrow(analysis_data)
p <- length(coef(model))

cook_threshold <- 4 / n

leverage_threshold <- 2 * p / n

residual_threshold <- 2

dfbeta_threshold <- 2 / sqrt(n)

# ------------------------------------------------------------
# 7. Build diagnostics table
# ------------------------------------------------------------

diagnostics <- analysis_data %>%
  mutate(
    cooks_distance = cooks_distance,
    leverage = leverage,
    standardized_residual = standardized_residuals,
    studentized_residual = studentized_residuals,
    DFBETA_CAF = dfbeta_caf,
    
    flag_cook = cooks_distance > cook_threshold,
    
    flag_leverage = leverage > leverage_threshold,
    
    flag_residual =
      abs(standardized_residual) > residual_threshold,
    
    flag_DFBETA =
      abs(DFBETA_CAF) > dfbeta_threshold
  ) %>%
  arrange(desc(cooks_distance))

# ------------------------------------------------------------
# 8. Print diagnostic thresholds
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Influence thresholds\n")
cat("============================================================\n\n")

cat(
  "Cook's distance threshold:",
  round(cook_threshold, 4),
  "\n"
)

cat(
  "Leverage threshold:",
  round(leverage_threshold, 4),
  "\n"
)

cat(
  "Standardized residual threshold:",
  residual_threshold,
  "\n"
)

cat(
  "DFBETA threshold:",
  round(dfbeta_threshold, 4),
  "\n"
)

# ------------------------------------------------------------
# 9. Print most influential observations
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Top observations by Cook's distance\n")
cat("============================================================\n\n")

print(
  diagnostics %>%
    select(
      patient_id,
      histology,
      CAF_associated_score,
      Wong_stemness,
      cooks_distance,
      leverage,
      standardized_residual,
      studentized_residual,
      DFBETA_CAF
    ) %>%
    head(10)
)

# ------------------------------------------------------------
# 10. Flag potentially influential observations
# ------------------------------------------------------------

flagged <- diagnostics %>%
  filter(
    flag_cook |
      flag_leverage |
      flag_residual |
      flag_DFBETA
  )

cat("\n============================================================\n")
cat("Potentially influential observations\n")
cat("============================================================\n\n")

if (nrow(flagged) == 0) {
  
  cat("No observations exceeded the predefined heuristic thresholds.\n")
  
} else {
  
  print(
    flagged %>%
      select(
        patient_id,
        histology,
        CAF_associated_score,
        Wong_stemness,
        cooks_distance,
        leverage,
        standardized_residual,
        studentized_residual,
        DFBETA_CAF,
        flag_cook,
        flag_leverage,
        flag_residual,
        flag_DFBETA
      )
  )
}

# ------------------------------------------------------------
# 11. Leave-one-out sensitivity analysis
# ------------------------------------------------------------

loo_results <- vector(
  "list",
  length = nrow(analysis_data)
)

for (i in seq_len(nrow(analysis_data))) {
  
  loo_data <- analysis_data[-i, ]
  
  loo_model <- lm(
    Wong_stemness ~ CAF_associated_score,
    data = loo_data
  )
  
  loo_pearson <- cor.test(
    loo_data$CAF_associated_score,
    loo_data$Wong_stemness,
    method = "pearson"
  )
  
  loo_spearman <- cor.test(
    loo_data$CAF_associated_score,
    loo_data$Wong_stemness,
    method = "spearman",
    exact = FALSE
  )
  
  loo_coef <- coef(summary(loo_model))
  
  loo_results[[i]] <- data.frame(
    
    omitted_patient =
      as.character(
        analysis_data$patient_id[i]
      ),
    
    omitted_histology =
      as.character(
        analysis_data$histology[i]
      ),
    
    n =
      nrow(loo_data),
    
    Pearson_r =
      unname(
        loo_pearson$estimate
      ),
    
    Pearson_p =
      loo_pearson$p.value,
    
    Spearman_rho =
      unname(
        loo_spearman$estimate
      ),
    
    Spearman_p =
      loo_spearman$p.value,
    
    CAF_slope =
      loo_coef[
        "CAF_associated_score",
        "Estimate"
      ],
    
    CAF_slope_p =
      loo_coef[
        "CAF_associated_score",
        "Pr(>|t|)"
      ]
    
  )
}

loo_results <- bind_rows(loo_results)

# ------------------------------------------------------------
# 12. Add change relative to full model
# ------------------------------------------------------------

full_pearson <- cor.test(
  analysis_data$CAF_associated_score,
  analysis_data$Wong_stemness,
  method = "pearson"
)

full_spearman <- cor.test(
  analysis_data$CAF_associated_score,
  analysis_data$Wong_stemness,
  method = "spearman",
  exact = FALSE
)

full_model_coef <- coef(summary(model))

full_slope <-
  full_model_coef[
    "CAF_associated_score",
    "Estimate"
  ]

loo_results <- loo_results %>%
  mutate(
    
    Pearson_change =
      Pearson_r -
      unname(full_pearson$estimate),
    
    Spearman_change =
      Spearman_rho -
      unname(full_spearman$estimate),
    
    CAF_slope_change =
      CAF_slope -
      full_slope
  )

# ------------------------------------------------------------
# 13. Print LOO results
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Leave-one-out Pearson results\n")
cat("============================================================\n\n")

print(
  loo_results %>%
    select(
      omitted_patient,
      omitted_histology,
      Pearson_r,
      Pearson_change
    ) %>%
    arrange(Pearson_r)
)

cat("\n============================================================\n")
cat("LOO Pearson range\n")
cat("============================================================\n\n")

cat(
  "Minimum:",
  round(min(loo_results$Pearson_r), 4),
  "\n"
)

cat(
  "Maximum:",
  round(max(loo_results$Pearson_r), 4),
  "\n"
)

# ------------------------------------------------------------
# 14. Identify largest LOO changes
# ------------------------------------------------------------

largest_pearson_changes <- loo_results %>%
  arrange(desc(abs(Pearson_change))) %>%
  head(5)

cat("\n============================================================\n")
cat("Largest absolute changes in Pearson r\n")
cat("============================================================\n\n")

print(
  largest_pearson_changes %>%
    select(
      omitted_patient,
      omitted_histology,
      Pearson_r,
      Pearson_change
    )
)

# ------------------------------------------------------------
# 15. Save diagnostics table
# ------------------------------------------------------------

diagnostics_file <- file.path(
  analysis_dir,
  "GSE33532_CAF_Wong_influence_diagnostics.csv"
)

write.csv(
  diagnostics,
  diagnostics_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 16. Save leave-one-out results
# ------------------------------------------------------------

loo_file <- file.path(
  analysis_dir,
  "GSE33532_CAF_Wong_leave_one_out.csv"
)

write.csv(
  loo_results,
  loo_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 17. Save influence summary
# ------------------------------------------------------------

summary_file <- file.path(
  analysis_dir,
  "GSE33532_CAF_Wong_influence_summary.txt"
)

summary_lines <- c(
  
  "GSE33532 CAF–Wong Influence and Sensitivity Analysis",
  "=====================================================",
  "",
  paste(
    "Patient-level sample size:",
    nrow(analysis_data)
  ),
  "",
  "Full-model Pearson correlation:",
  paste(
    "r =",
    round(
      unname(full_pearson$estimate),
      4
    ),
    ", p =",
    signif(
      full_pearson$p.value,
      4
    )
  ),
  "",
  "Full-model Spearman correlation:",
  paste(
    "rho =",
    round(
      unname(full_spearman$estimate),
      4
    ),
    ", p =",
    signif(
      full_spearman$p.value,
      4
    )
  ),
  "",
  "Full-model CAF regression slope:",
  paste(
    "beta =",
    round(
      full_slope,
      4
    )
  ),
  "",
  "Influence thresholds:",
  paste(
    "Cook's distance >",
    round(
      cook_threshold,
      4
    )
  ),
  paste(
    "Leverage >",
    round(
      leverage_threshold,
      4
    )
  ),
  paste(
    "|standardized residual| >",
    residual_threshold
  ),
  paste(
    "|DFBETA| >",
    round(
      dfbeta_threshold,
      4
    )
  ),
  "",
  paste(
    "Number of flagged observations:",
    nrow(flagged)
  ),
  "",
  "LOO Pearson range:",
  paste(
    round(
      min(loo_results$Pearson_r),
      4
    ),
    "to",
    round(
      max(loo_results$Pearson_r),
      4
    )
  ),
  "",
  "Largest absolute LOO Pearson change:",
  paste(
    round(
      max(
        abs(
          loo_results$Pearson_change
        )
      ),
      4
    )
  ),
  "",
  "Interpretation:",
  "The purpose of this analysis is sensitivity assessment,",
  "not removal of observations.",
  "Potential influence is evaluated using predefined diagnostic",
  "criteria and leave-one-out analyses.",
  "",
  "The CAF–Wong relationship should not be interpreted as",
  "supported evidence of an association unless the sensitivity",
  "analysis and primary statistical results support that conclusion."
)

writeLines(
  summary_lines,
  summary_file
)

# ------------------------------------------------------------
# 18. Figure: Cook's distance
# ------------------------------------------------------------

cook_plot_data <- diagnostics %>%
  mutate(
    patient_id = factor(
      patient_id,
      levels = patient_id[
        order(cooks_distance)
      ]
    )
  )

p_cook <- ggplot(
  cook_plot_data,
  aes(
    x = patient_id,
    y = cooks_distance
  )
) +
  geom_col() +
  geom_hline(
    yintercept = cook_threshold,
    linetype = "dashed"
  ) +
  labs(
    title = "Cook's Distance: CAF–Wong Model",
    x = "Patient",
    y = "Cook's distance"
  ) +
  theme_bw() +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )

cook_figure <- file.path(
  figure_dir,
  "GSE33532_CAF_Wong_Cooks_distance.png"
)

ggsave(
  cook_figure,
  p_cook,
  width = 9,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 19. Figure: Leave-one-out Pearson correlation
# ------------------------------------------------------------

loo_plot_data <- loo_results %>%
  mutate(
    omitted_patient = factor(
      omitted_patient,
      levels = omitted_patient[
        order(Pearson_r)
      ]
    )
  )

p_loo <- ggplot(
  loo_plot_data,
  aes(
    x = omitted_patient,
    y = Pearson_r
  )
) +
  geom_point(
    size = 3
  ) +
  geom_hline(
    yintercept =
      unname(full_pearson$estimate),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted"
  ) +
  labs(
    title =
      "Leave-One-Out Pearson Correlation: CAF–Wong",
    x = "Omitted patient",
    y = "Pearson r"
  ) +
  theme_bw() +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )

loo_figure <- file.path(
  figure_dir,
  "GSE33532_CAF_Wong_leave_one_out_Pearson.png"
)

ggsave(
  loo_figure,
  p_loo,
  width = 9,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 20. Final verification
# ------------------------------------------------------------

expected_files <- c(
  diagnostics_file,
  loo_file,
  summary_file,
  cook_figure,
  loo_figure
)

cat("\n============================================================\n")
cat("Script 23 verification\n")
cat("============================================================\n\n")

for (f in expected_files) {
  
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

cat("\nScript 23 completed.\n")



