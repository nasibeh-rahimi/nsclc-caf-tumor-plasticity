############################################################
# 00_project_setup.R
# NSCLC Plasticity Framework
#
# Purpose:
# Central project configuration script.
#
# Responsibilities:
#   1. Define project paths
#   2. Create required directories
#   3. Load required R packages
#   4. Set reproducible analysis options
#
# All downstream analysis scripts should source this file.
############################################################
library(readxl)
library(tidyr)
library(janitor)

############################################################
# 1. PROJECT ROOT
############################################################
# Identify the project root from the current working directory.
# The project root must contain the "scripts" directory.

current_dir <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

if (basename(current_dir) == "scripts") {
  current_dir <- dirname(current_dir)
}

if (!dir.exists(file.path(current_dir, "scripts"))) {
  stop(
    paste(
      "Project root not found.",
      "Please set the working directory to the project folder",
      "before running the analysis scripts."
    )
  )
}

project_dir <- current_dir


############################################################
# 2. PROJECT DIRECTORIES
############################################################

# Main directories
data_dir      <- file.path(project_dir, "data")
processed_dir <- file.path(data_dir, "processed")

script_dir <- file.path(
  project_dir,
  "scripts"
)

result_dir <- file.path(
  project_dir,
  "results"
)

figure_dir <- file.path(
  project_dir,
  "figures"
)

table_dir <- file.path(
  project_dir,
  "tables"
)

# Compatibility aliases for downstream scripts
results_dir <- result_dir
figures_dir <- figure_dir

# Result subdirectories
expression_dir <- file.path(
  result_dir,
  "expression"
)

plasticity_dir <- file.path(
  result_dir,
  "plasticity"
)

correlation_dir <- file.path(
  result_dir,
  "correlation"
)

enrichment_dir <- file.path(
  result_dir,
  "enrichment"
)

validation_dir <- file.path(
  result_dir,
  "validation"
)

log_dir <- file.path(
  result_dir,
  "logs"
)


############################################################
# 3. CREATE DIRECTORIES
############################################################

project_dirs <- c(
  
  data_dir,
  processed_dir,
  
  script_dir,
  
  result_dir,
  expression_dir,
  plasticity_dir,
  correlation_dir,
  enrichment_dir,
  validation_dir,
  log_dir,
  
  figure_dir,
  table_dir
)

for (dir in project_dirs) {
  
  if (!dir.exists(dir)) {
    
    dir.create(
      dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
  }
}


############################################################
# 4. LOAD REQUIRED PACKAGES
############################################################

required_packages <- c(
  
  # Microarray / Bioconductor
  "affy",
  "Biobase",
  "limma",
  "AnnotationDbi",
  "hgu133plus2.db",
  "GSVA",
  
  # Data manipulation
  "dplyr",
  "readr",
  "readxl",
  "tidyr",
  "purrr",
  "janitor",
  
  # Gene signatures
  "msigdbr",
  
  # Visualization
  "ggplot2",
  "pheatmap",
  
  # Matrix / numerical analysis
  "matrixStats"
)


############################################################
# 5. CHECK AND LOAD PACKAGES
############################################################

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    stop(
      paste0(
        "Required package '",
        pkg,
        "' is not installed. ",
        "Please install it before running the project."
      )
    )
    
  }
  
  suppressPackageStartupMessages(
    library(
      pkg,
      character.only = TRUE
    )
  )
}


############################################################
# 6. GENERAL R OPTIONS
############################################################

options(
  stringsAsFactors = FALSE
)


############################################################
# 7. PROJECT INFORMATION
############################################################

cat("\n")
cat("============================================================\n")
cat(" NSCLC Plasticity Framework\n")
cat(" Project setup loaded successfully\n")
cat("============================================================\n")

cat("\nProject directory:\n")
cat(project_dir, "\n")

cat("\nScripts:\n")
cat(script_dir, "\n")

cat("\nProcessed data:\n")
cat(processed_dir, "\n")

cat("\nResults:\n")
cat(result_dir, "\n")

cat("\nFigures:\n")
cat(figure_dir, "\n")

cat("\nTables:\n")
cat(table_dir, "\n")

cat("\nProject ready.\n")
cat("============================================================\n")