############################################################
## 07_GSE149955_metadata.R
## Build sample metadata table from GEO Series Matrix
############################################################

source("scripts/00_project_setup.R")

library(stringr)
library(dplyr)

#=========================================
# Read Series Matrix
#=========================================

series_file <- "data/metadata/GSE149955_series_matrix.txt"

lines <- readLines(series_file)

cat("Series Matrix loaded\n")

#=========================================
# Helper function
#=========================================

extract_field <- function(field){
  
  x <- lines[grep(paste0("^", field), lines)]
  
  x <- strsplit(x, "\t")[[1]]
  
  x <- x[-1]
  
  x <- gsub('"', "", x)
  
  return(x)
  
}

#=========================================
# Extract metadata
#=========================================

gsm <- extract_field("!Sample_geo_accession")

title <- extract_field("!Sample_title")

source_name <- extract_field("!Sample_source_name_ch1")

#=========================================
# Build metadata table
#=========================================

sample_info <- data.frame(
  
  sample = gsm,
  
  title = title,
  
  source = source_name,
  
  stringsAsFactors = FALSE
  
)

#=========================================
# Group
#=========================================

sample_info$group <- ifelse(
  
  grepl("Cancer associated fibroblast",
        sample_info$title,
        ignore.case = TRUE),
  
  "CAF",
  
  "NF"
  
)

#=========================================
# Patient ID
#=========================================

sample_info$patient <- str_extract(
  
  sample_info$title,
  
  "[0-9]{4,}(_[12])?"
  
)

#=========================================
# Reorder columns
#=========================================
sample_info <- sample_info[
  ,
  c("sample",
    "patient",
    "group",
    "title",
    "source")
]

#=========================================
# Check metadata
#=========================================

print(head(sample_info))

cat("\n")

cat("Samples :", nrow(sample_info), "\n")

cat("Patients:", length(unique(sample_info$patient)), "\n")

print(table(sample_info$group))

#=========================================
# Save
#=========================================

write.csv(
  
  sample_info,
  
  file.path(
    
    table_dir,
    
    "GSE149955_sample_info.csv"
    
  ),
  
  row.names = FALSE
  
)

cat("\nMetadata saved successfully.\n")


############################################################
# Metadata validation
############################################################

# Every patient should have exactly two samples
patient_counts <- table(sample_info$patient)

stopifnot(all(patient_counts == 2))

# Every patient should have one NF and one CAF
group_check <- table(sample_info$patient,
                     sample_info$group)

stopifnot(all(group_check == 1))

cat("✓ Metadata validation passed\n")
