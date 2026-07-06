#!/usr/bin/env Rscript

# ==============================================================================
# Median cytotoxicity score per tetramer-specific T cell population
# ==============================================================================
#
# Purpose:
#   This script calculates the median cytotoxicity score for each antigen-specific
#   T cell population within each sample.
#
#   The output table contains one row per:
#       sample x cohort x antigen
#
#   This donor-level summary table is intended for downstream plotting and
#   statistical testing.
#
# Statistical analysis:
#   No statistical testing is performed in this R script.
#
#   Statistical comparisons for the corresponding figure were performed in
#   BioRender.com using a two-way ANOVA with Bonferroni multiple-comparisons
#   testing.
#
# Input:
#   A Seurat RDS object containing baseline EBV tetramer-enriched T cells.
#
# Required metadata columns:
#   - sample
#   - cohort
#   - antigen
#   - cytotoxicity_score
#
# Output:
#   A CSV file containing:
#   - sample
#   - cohort
#   - antigen
#   - n_cells
#   - median_cytotoxicity_score
#   - mean_cytotoxicity_score
#
# ==============================================================================


# ----------------------------- #
# 0. Load packages
# ----------------------------- #

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  
  # Input Seurat RDS.
  input_rds = "path/to/Publication_RDS/Baseline_EBV_integrated_seurat.rds",
  
  # Output directory.
  out_dir = "path/to/output/Stats",
  
  # Output files.
  output_csv = "median_cytotoxicity_score_per_sample_antigen.csv",
  session_info_file = "sessionInfo_median_cytotoxicity_score_per_sample_antigen.txt",
  
  # Metadata columns.
  sample_col = "sample",
  cohort_col = "cohort",
  antigen_col = "antigen",
  cyto_col = "cytotoxicity_score",
  
  # Optional cohort order.
  cohort_order = c("Control", "MS"),
  
  # Statistical note.
  stats_note = paste(
    "Statistical comparisons for the corresponding figure were performed in",
    "BioRender.com using a two-way ANOVA with Bonferroni multiple-comparisons",
    "testing."
  )
)


# ----------------------------- #
# 2. Helper functions
# ----------------------------- #

create_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}


check_required_columns <- function(df, required_cols, object_name = "dataframe") {
  
  missing_cols <- setdiff(required_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    
    possible_score_cols <- grep(
      "cyto|score|module",
      colnames(df),
      ignore.case = TRUE,
      value = TRUE
    )
    
    possible_antigen_cols <- grep(
      "antigen|tetramer|specific",
      colnames(df),
      ignore.case = TRUE,
      value = TRUE
    )
    
    stop(
      "Missing required columns in ",
      object_name,
      ": ",
      paste(missing_cols, collapse = ", "),
      "\n\nPossible cytotoxicity/module-score columns:\n",
      paste(possible_score_cols, collapse = "\n"),
      "\n\nPossible antigen/tetramer columns:\n",
      paste(possible_antigen_cols, collapse = "\n")
    )
  }
  
  invisible(TRUE)
}


clean_cohort <- function(x) {
  
  x_clean <- str_to_lower(str_trim(as.character(x)))
  
  case_when(
    x_clean %in% c("ms", "pwms", "multiple sclerosis", "rrms", "ppms", "spms") ~ "MS",
    x_clean %in% c(
      "control",
      "controls",
      "healthy control",
      "healthy controls",
      "hc",
      "nms",
      "non-ms",
      "non_ms"
    ) ~ "Control",
    TRUE ~ as.character(x)
  )
}


clean_antigen <- function(x) {
  str_squish(as.character(x))
}


# ----------------------------- #
# 3. Prepare output paths
# ----------------------------- #

create_dir(config$out_dir)

output_csv_path <- file.path(config$out_dir, config$output_csv)
session_info_path <- file.path(config$out_dir, config$session_info_file)


# ----------------------------- #
# 4. Read Seurat RDS
# ----------------------------- #

if (!file.exists(config$input_rds)) {
  stop("Input RDS file does not exist: ", config$input_rds)
}

message("Reading Seurat RDS:")
message(config$input_rds)

merged <- readRDS(config$input_rds)

if (!inherits(merged, "Seurat")) {
  stop("Input RDS is not a Seurat object.")
}

meta <- merged@meta.data %>%
  rownames_to_column("cell_id")

message("Cells in Seurat object: ", nrow(meta))
message("Metadata columns available: ", ncol(meta))


# ----------------------------- #
# 5. Check required metadata columns
# ----------------------------- #

required_cols <- c(
  config$sample_col,
  config$cohort_col,
  config$antigen_col,
  config$cyto_col
)

check_required_columns(
  df = meta,
  required_cols = required_cols,
  object_name = "Seurat metadata"
)


# ----------------------------- #
# 6. Clean metadata
# ----------------------------- #

cyto_df <- meta %>%
  transmute(
    sample = as.character(.data[[config$sample_col]]),
    cohort = clean_cohort(.data[[config$cohort_col]]),
    antigen = clean_antigen(.data[[config$antigen_col]]),
    cytotoxicity_score = suppressWarnings(as.numeric(.data[[config$cyto_col]]))
  ) %>%
  filter(
    !is.na(sample),
    sample != "",
    !is.na(cohort),
    cohort != "",
    !is.na(antigen),
    antigen != "",
    !is.na(cytotoxicity_score)
  ) %>%
  mutate(
    cohort = factor(
      cohort,
      levels = config$cohort_order
    )
  )

message("Cells retained after filtering: ", nrow(cyto_df))

if (nrow(cyto_df) == 0) {
  stop("No usable cells remained after filtering.")
}


# ----------------------------- #
# 7. Calculate median cytotoxicity score
# ----------------------------- #
#
# One row is produced per sample x cohort x antigen.
#
# This is the table used for BioRender plotting/statistical testing.
#
# No statistical testing is performed here.
#
# BioRender statistical method:
#   Two-way ANOVA with Bonferroni multiple-comparisons testing.
#
# ----------------------------- #

median_cyto_by_antigen <- cyto_df %>%
  group_by(sample, cohort, antigen) %>%
  summarise(
    n_cells = n(),
    median_cytotoxicity_score = median(cytotoxicity_score, na.rm = TRUE),
    mean_cytotoxicity_score = mean(cytotoxicity_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(cohort, sample, antigen)

print(median_cyto_by_antigen)


# ----------------------------- #
# 8. Save outputs
# ----------------------------- #

write_csv(
  median_cyto_by_antigen,
  output_csv_path
)

writeLines(
  c(
    "Median cytotoxicity score per sample and antigen",
    "================================================",
    "",
    "Input RDS:",
    config$input_rds,
    "",
    "Output CSV:",
    output_csv_path,
    "",
    "Metadata columns used:",
    paste0("sample column: ", config$sample_col),
    paste0("cohort column: ", config$cohort_col),
    paste0("antigen column: ", config$antigen_col),
    paste0("cytotoxicity score column: ", config$cyto_col),
    "",
    "Analysis performed in this script:",
    "Median cytotoxicity score was calculated per sample, cohort and antigen.",
    "",
    "Statistical analysis note:",
    config$stats_note,
    "",
    "Important:",
    "No statistical testing was performed in this R script.",
    "",
    "Cells in Seurat object:",
    as.character(nrow(meta)),
    "",
    "Cells retained after filtering:",
    as.character(nrow(cyto_df)),
    "",
    "Rows exported:",
    as.character(nrow(median_cyto_by_antigen)),
    "",
    "Exported data preview:",
    capture.output(print(median_cyto_by_antigen)),
    "",
    "Session information:",
    capture.output(sessionInfo())
  ),
  session_info_path
)

message("\nMedian cytotoxicity score export complete.")
message("Saved CSV: ", output_csv_path)
message("Saved session info: ", session_info_path)
message("")
message("No statistical testing was performed in this R script.")
message(config$stats_note)