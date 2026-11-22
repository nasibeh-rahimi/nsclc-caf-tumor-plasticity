
# ============================================================
# Script 27: Final Project Figures
# Purpose: Generate three presentation-ready figures
#          from previously generated analysis outputs.
#
# This script does not rerun statistical analyses.
# Scripts 00-26 remain unchanged.
# ============================================================


# 1. Load project setup ---------------------------------------

if (!file.exists(file.path("scripts", "00_project_setup.R"))) {
  stop(
    "Please set the working directory to the project root: ",
    "D:/NSCLC_Plasticity2"
  )
}

source(file.path("scripts", "00_project_setup.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})


# 2. Define input and output directories ----------------------

analysis_dir <- file.path(result_dir, "analysis")

if (!dir.exists(analysis_dir)) {
  stop("Analysis directory not found: ", analysis_dir)
}

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 3. Figure 1: CAF-associated program in discovery cohort ------

caf_file <- file.path(
  analysis_dir,
  "GSE149955_CAF_signature_discovery_scores.csv"
)

if (!file.exists(caf_file)) {
  stop("Missing input file: ", caf_file)
}

caf_data <- read.csv(
  caf_file,
  stringsAsFactors = FALSE
)

required_caf_columns <- c(
  "CAF_associated_score",
  "group"
)

if (!all(required_caf_columns %in% names(caf_data))) {
  stop("Required columns are missing from the CAF discovery table.")
}

if (anyNA(caf_data$CAF_associated_score)) {
  stop("Missing values detected in CAF-associated scores.")
}

caf_data$group <- factor(
  caf_data$group,
  levels = c("NF", "CAF")
)

# Make jittered point positions reproducible
set.seed(123)


p1 <- ggplot(
  caf_data,
  aes(x = group, y = CAF_associated_score)
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.10,
    size = 2.2,
    alpha = 0.75
  ) +
  labs(
    title = "CAF-associated transcriptional program",
    subtitle = "GSE149955 discovery cohort: 25 matched NF–CAF pairs",
    x = NULL,
    y = "CAF-associated score"
  ) +
  theme_classic(base_size = 13)

ggsave(
  filename = file.path(
    figure_dir,
    "Figure1_CAF_program_in_discovery_cohort.png"
  ),
  plot = p1,
  width = 7,
  height = 5.5,
  dpi = 300
)


# 4. Figure 2: Primary CAF–EMT association ---------------------

patient_file <- file.path(
  analysis_dir,
  "GSE33532_patient_level_analysis_table.csv"
)

if (!file.exists(patient_file)) {
  stop("Missing input file: ", patient_file)
}

patient_data <- read.csv(
  patient_file,
  stringsAsFactors = FALSE
)

required_patient_columns <- c(
  "patient_id",
  "CAF_associated_score",
  "Byers76_EMT",
  "Wong_stemness"
)

if (!all(required_patient_columns %in% names(patient_data))) {
  stop("Required columns are missing from the patient-level table.")
}

if (nrow(patient_data) != 20) {
  stop(
    "Expected 20 independent patients; found ",
    nrow(patient_data),
    "."
  )
}

if (anyDuplicated(patient_data$patient_id)) {
  stop("Duplicate patient IDs detected.")
}

score_columns <- c(
  "CAF_associated_score",
  "Byers76_EMT",
  "Wong_stemness"
)

if (anyNA(patient_data[, score_columns])) {
  stop("Missing values detected in patient-level scores.")
}

# Primary Pearson correlation
ct_emt <- cor.test(
  patient_data$CAF_associated_score,
  patient_data$Byers76_EMT,
  method = "pearson"
)

p2 <- ggplot(
  patient_data,
  aes(
    x = CAF_associated_score,
    y = Byers76_EMT
  )
) +
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf(
      "Pearson r = %.3f\np = %.3f",
      unname(ct_emt$estimate),
      ct_emt$p.value
    ),
    hjust = 1.05,
    vjust = 1.2,
    size = 4.2
  ) +
  labs(
    title = "CAF-associated program versus EMT-associated state",
    subtitle = "GSE33532 patient-level primary analysis (n = 20)",
    x = "CAF-associated transcriptional score",
    y = "Byers76 EMT-associated score"
  ) +
  theme_classic(base_size = 13)

ggsave(
  filename = file.path(
    figure_dir,
    "Figure2_CAF_vs_Byers76_EMT_primary.png"
  ),
  plot = p2,
  width = 7,
  height = 5.5,
  dpi = 300
)


# 5. Figure 3: Secondary CAF–Wong association ------------------

ct_wong <- cor.test(
  patient_data$CAF_associated_score,
  patient_data$Wong_stemness,
  method = "pearson"
)

p3 <- ggplot(
  patient_data,
  aes(
    x = CAF_associated_score,
    y = Wong_stemness
  )
) +
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf(
      "Pearson r = %.3f\np = %.3f",
      unname(ct_wong$estimate),
      ct_wong$p.value
    ),
    hjust = 1.05,
    vjust = 1.2,
    size = 4.2
  ) +
  labs(
    title = "CAF-associated program versus broad stemness",
    subtitle = "GSE33532 patient-level secondary analysis (n = 20)",
    x = "CAF-associated transcriptional score",
    y = "Wong stemness-associated score"
  ) +
  theme_classic(base_size = 13)

ggsave(
  filename = file.path(
    figure_dir,
    "Figure3_CAF_vs_Wong_stemness_secondary.png"
  ),
  plot = p3,
  width = 7,
  height = 5.5,
  dpi = 300
)


# 6. Verify generated figures ---------------------------------

final_figures <- c(
  "Figure1_CAF_program_in_discovery_cohort.png",
  "Figure2_CAF_vs_Byers76_EMT_primary.png",
  "Figure3_CAF_vs_Wong_stemness_secondary.png"
)

cat("\nFinal figure generation completed.\n\n")

for (f in final_figures) {
  output_path <- file.path(figure_dir, f)
  
  cat(
    ifelse(file.exists(output_path), "OK: ", "MISSING: "),
    output_path,
    "\n"
  )
}
