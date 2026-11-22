############################################################
# 20_GSE33532_Primary_CAF_Plasticity_Association.R
#
# Primary question:
#   Is CAF-associated transcriptional activity associated
#   with tumor plasticity/EMT in NSCLC?
#
# Primary variables:
#   CAF_associated_score
#   Byers76_EMT
#
# Unit of analysis:
#   Patient level (n = 20)
#
# Histology:
#   Examined as a potential confounding/context variable
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

figure_dir <- file.path(
  project_dir,
  "figures"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 2. Load patient-level data
############################################################

input_file <- file.path(
  analysis_dir,
  "GSE33532_patient_level_analysis_table.csv"
)

stopifnot(file.exists(input_file))

patient_level <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Patient-level table loaded.\n")
cat("Patients:", nrow(patient_level), "\n")


############################################################
# 3. Validate input
############################################################

required_columns <- c(
  "patient_id",
  "histology",
  "CAF_associated_score",
  "Byers76_EMT"
)

stopifnot(
  all(required_columns %in% colnames(patient_level))
)

stopifnot(
  nrow(patient_level) == 20
)

stopifnot(
  length(unique(patient_level$patient_id)) == 20
)

stopifnot(
  !anyNA(patient_level[, required_columns])
)

patient_level$histology <- factor(
  patient_level$histology,
  levels = c(
    "Adenocarcinoma",
    "Mixed",
    "Squamous"
  )
)

stopifnot(
  !anyNA(patient_level$histology)
)


############################################################
# 4. Descriptive summaries
############################################################

cat("\n============================================\n")
cat("DESCRIPTIVE SUMMARY\n")
cat("============================================\n")

cat("Histology:\n")
print(table(patient_level$histology))

cat("\nCAF-associated score:\n")
print(summary(patient_level$CAF_associated_score))

cat("\nByers76 EMT score:\n")
print(summary(patient_level$Byers76_EMT))


############################################################
# 5. Pearson correlation
############################################################

pearson_test <- cor.test(
  patient_level$CAF_associated_score,
  patient_level$Byers76_EMT,
  method = "pearson"
)

cat("\n============================================\n")
cat("PEARSON CORRELATION\n")
cat("============================================\n")

print(pearson_test)


############################################################
# 6. Spearman correlation
############################################################

spearman_test <- cor.test(
  patient_level$CAF_associated_score,
  patient_level$Byers76_EMT,
  method = "spearman",
  exact = FALSE
)

cat("\n============================================\n")
cat("SPEARMAN CORRELATION\n")
cat("============================================\n")

print(spearman_test)


############################################################
# 7. Linear model without histology adjustment
############################################################

primary_lm <- lm(
  Byers76_EMT ~ CAF_associated_score,
  data = patient_level
)

cat("\n============================================\n")
cat("UNADJUSTED LINEAR MODEL\n")
cat("============================================\n")

print(summary(primary_lm))


############################################################
# 8. Histology-adjusted sensitivity model
############################################################

histology_lm <- lm(
  Byers76_EMT ~ CAF_associated_score + histology,
  data = patient_level
)

cat("\n============================================\n")
cat("HISTOLOGY-ADJUSTED LINEAR MODEL\n")
cat("============================================\n")

print(summary(histology_lm))


############################################################
# 9. Scatterplot: CAF vs EMT
############################################################

png(
  filename = file.path(
    figure_dir,
    "GSE33532_CAF_vs_Byers76_EMT_patient_level.png"
  ),
  width = 1800,
  height = 1400,
  res = 220
)

plot(
  patient_level$CAF_associated_score,
  patient_level$Byers76_EMT,
  pch = 19,
  xlab = "CAF-associated transcriptional program score",
  ylab = "Byers76 EMT score",
  main = "CAF-associated program vs NSCLC EMT"
)

abline(
  primary_lm,
  lwd = 2
)

text(
  patient_level$CAF_associated_score,
  patient_level$Byers76_EMT,
  labels = patient_level$patient_id,
  pos = 3,
  cex = 0.7
)

dev.off()


############################################################
# 10. Histology-aware scatterplot
############################################################

histology_symbols <- c(
  "Adenocarcinoma" = 16,
  "Mixed" = 17,
  "Squamous" = 15
)

png(
  filename = file.path(
    figure_dir,
    "GSE33532_CAF_vs_Byers76_EMT_histology.png"
  ),
  width = 1800,
  height = 1400,
  res = 220
)

plot(
  patient_level$CAF_associated_score,
  patient_level$Byers76_EMT,
  pch = histology_symbols[
    as.character(patient_level$histology)
  ],
  xlab = "CAF-associated transcriptional program score",
  ylab = "Byers76 EMT score",
  main = "CAF-associated program vs NSCLC EMT",
  cex = 1.2
)

abline(
  primary_lm,
  lwd = 2
)

legend(
  "topleft",
  legend = names(histology_symbols),
  pch = histology_symbols,
  title = "Histology",
  bty = "n"
)

dev.off()


############################################################
# 11. Save statistical results
############################################################

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
  )
)

write.csv(
  correlation_results,
  file.path(
    analysis_dir,
    "GSE33532_CAF_vs_Byers76_EMT_correlations.csv"
  ),
  row.names = FALSE
)

capture.output(
  summary(primary_lm),
  file = file.path(
    analysis_dir,
    "GSE33532_CAF_vs_Byers76_EMT_linear_model.txt"
  )
)

capture.output(
  summary(histology_lm),
  file = file.path(
    analysis_dir,
    "GSE33532_CAF_vs_Byers76_EMT_histology_adjusted_model.txt"
  )
)


############################################################
# 12. Final confirmation
############################################################

cat("\n============================================\n")
cat("SCRIPT 20 COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

cat(
  "Patient-level observations:",
  nrow(patient_level),
  "\n"
)

cat(
  "Primary analysis: CAF-associated score vs Byers76 EMT\n"
)

cat(
  "Figures saved to:\n",
  figure_dir,
  "\n"
)

cat(
  "Statistical results saved to:\n",
  analysis_dir,
  "\n"
)