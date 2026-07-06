#!/usr/bin/env Rscript

# ==============================================================================
# Dot-level frequency tables by HLA A2/B7/DR15 combination
# ==============================================================================
#
# Purpose:
#   This script creates source-data tables for latent-specific CD8+ T-cell
#   frequency by HLA-A2, HLA-B7 and HLA-DR15 combination.
#
#   It exports:
#     1. One row per dot for the pooled HLA-combination graph.
#     2. One row per dot for the cohort-split MS vs Control graph.
#     3. Counts by HLA condition.
#     4. Counts by HLA condition and cohort.
#     5. Descriptive summaries by HLA condition.
#     6. Descriptive summaries by HLA condition and cohort.
# Input: 
#     - Uses correlation_variables also available in source_data.xlsx sheet FigS1A-C
# HLA grouping:
#   The following non-triple-negative groups are retained:
#     - A2+B7+DR15+
#     - A2+B7+DR15-
#     - A2+B7-DR15+
#     - A2+B7-DR15-
#     - A2-B7+DR15+
#     - A2-B7+DR15-
#     - A2-B7-DR15+
#
#   Triple-negative participants are excluded:
#     - A2-B7-DR15-
#
# BioRender note:
#   The exported dot-level source-data tables from this script were imported into
#   BioRender.com for graphing and statistical testing.
#
# Notes:
#   - Paths are intentionally generic for GitHub.
#   - Update the config block before running.
#   - Raw frequency is exported as T_cell_freq.
#   - A log10-transformed plotting value is also exported as T_cell_freq_log10.
#   - No graphing or statistical testing is performed in this script.
#
# ==============================================================================


# ----------------------------- #
# 0. Load packages
# ----------------------------- #

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(openxlsx)
  library(readr)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  
  # Input Excel file
  input_file = "path/to/correlation_variables.xlsx",
  
  # Input sheet
  input_sheet = "Full",
  
  # Output directory
  output_dir = "outputs/HLA_A2_B7_DR15_dot_level_tables",
  
  # Column names in the input file
  frequency_col = "T_cell_freq",
  cohort_col = "cohort",
  hla_a2_col = "HLA-A2",
  hla_b7_col = "HLA-B7",
  hla_dr15_col = "HLA-DR15",
  
  # Candidate sample/donor column names.
  # The first one found will be used.
  sample_col_candidates = c("sample", "id", "ID", "donor_id", "participant_id"),
  
  # Pseudocount used only to create the plotting/log-transformed value.
  # Raw frequency remains unchanged.
  pseudocount = 0.001,
  
  # Cohort labels
  control_label = "Control",
  ms_label = "MS",
  
  # If TRUE, unresolved non-missing HLA strings are treated as positive.
  # Use FALSE for stricter handling.
  treat_unresolved_hla_as_positive = FALSE,
  
  # Output files
  pooled_dot_csv = "dot_level_frequencies_by_HLA_condition_pooled.csv",
  cohort_dot_csv = "dot_level_frequencies_by_HLA_condition_and_cohort.csv",
  counts_hla_csv = "counts_by_HLA_condition.csv",
  counts_hla_cohort_csv = "counts_by_HLA_condition_and_cohort.csv",
  summary_hla_csv = "summary_frequency_by_HLA_condition.csv",
  summary_hla_cohort_csv = "summary_frequency_by_HLA_condition_and_cohort.csv",
  excluded_rows_csv = "excluded_rows_missing_or_triple_negative.csv",
  workbook_xlsx = "HLA_A2_B7_DR15_dot_level_frequency_tables.xlsx",
  session_info_file = "sessionInfo_HLA_A2_B7_DR15_dot_level_tables.txt"
)


# ----------------------------- #
# 2. Helper functions
# ----------------------------- #

create_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
}


read_input_excel <- function(input_file, input_sheet) {
  
  if (!file.exists(input_file)) {
    stop("Input file does not exist: ", input_file)
  }
  
  readxl::read_excel(input_file, sheet = input_sheet)
}


check_required_columns <- function(df, required_cols) {
  
  missing_cols <- setdiff(required_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  invisible(TRUE)
}


find_first_existing_col <- function(df, candidates) {
  
  hit <- candidates[candidates %in% colnames(df)]
  
  if (length(hit) == 0) {
    return(NA_character_)
  }
  
  hit[1]
}


to_numeric_clean <- function(x) {
  
  z <- trimws(as.character(x))
  z[z %in% c("", "NA", "N/A", "na", "n/a", "NaN", "nan", "null", "NULL")] <- NA
  z <- gsub(",", "", z)
  z <- gsub("<", "", z)
  z <- gsub(">", "", z)
  
  suppressWarnings(as.numeric(z))
}


to_binary_hla <- function(x, treat_unresolved_as_positive = FALSE) {
  
  if (is.logical(x)) {
    return(ifelse(is.na(x), NA_real_, as.numeric(x)))
  }
  
  if (is.numeric(x)) {
    return(ifelse(is.na(x), NA_real_, ifelse(x > 0, 1, 0)))
  }
  
  z <- tolower(trimws(as.character(x)))
  z[z %in% c("", "na", "n/a", "nan", "null", "none")] <- NA
  
  out <- rep(NA_real_, length(z))
  
  out[z %in% c(
    "1", "yes", "y", "positive", "pos",
    "carrier", "present", "true", "+"
  )] <- 1
  
  out[z %in% c(
    "0", "no", "n", "negative", "neg",
    "non-carrier", "absent", "false", "-"
  )] <- 0
  
  unresolved <- is.na(out) & !is.na(z)
  
  if (any(unresolved) && isTRUE(treat_unresolved_as_positive)) {
    out[unresolved] <- 1
  }
  
  out
}


normalise_cohort <- function(x, control_label = "Control", ms_label = "MS") {
  
  z <- tolower(trimws(as.character(x)))
  z[z %in% c("", "na", "n/a", "nan", "null", "none")] <- NA
  
  out <- rep(NA_character_, length(z))
  
  out[z %in% c("1", "ms", "pwms", "rrms", "multiple sclerosis", "case", "patient")] <- ms_label
  
  out[z %in% c(
    "0", "control", "controls", "ctrl", "healthy control",
    "healthy controls", "hc", "nms", "non-ms", "nonms"
  )] <- control_label
  
  factor(out, levels = c(control_label, ms_label))
}


pm <- function(x) {
  
  case_when(
    x == 1 ~ "+",
    x == 0 ~ "-",
    TRUE ~ NA_character_
  )
}


format_range <- function(x1, x2, digits = 6) {
  
  if (is.na(x1) || is.na(x2)) {
    return(NA_character_)
  }
  
  paste0(round(x1, digits), " - ", round(x2, digits))
}


summarise_frequency <- function(df, group_cols) {
  
  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n = n(),
      n_zero = sum(T_cell_freq == 0, na.rm = TRUE),
      mean = mean(T_cell_freq, na.rm = TRUE),
      sd = sd(T_cell_freq, na.rm = TRUE),
      sem = sd / sqrt(n),
      median = median(T_cell_freq, na.rm = TRUE),
      q1 = quantile(T_cell_freq, 0.25, na.rm = TRUE),
      q3 = quantile(T_cell_freq, 0.75, na.rm = TRUE),
      min_value = min(T_cell_freq, na.rm = TRUE),
      max_value = max(T_cell_freq, na.rm = TRUE),
      IQR = format_range(q1, q3),
      min_max = format_range(min_value, max_value),
      .groups = "drop"
    ) %>%
    arrange(across(all_of(group_cols)))
}


# ----------------------------- #
# 3. Load data and check columns
# ----------------------------- #

create_dir(config$output_dir)

dt <- read_input_excel(
  input_file = config$input_file,
  input_sheet = config$input_sheet
)

dt <- as.data.frame(dt)
colnames(dt) <- trimws(colnames(dt))

required_cols <- c(
  config$frequency_col,
  config$cohort_col,
  config$hla_a2_col,
  config$hla_b7_col,
  config$hla_dr15_col
)

check_required_columns(dt, required_cols)

sample_col <- find_first_existing_col(
  df = dt,
  candidates = config$sample_col_candidates
)

if (is.na(sample_col)) {
  stop(
    "Could not find sample column. Tried: ",
    paste(config$sample_col_candidates, collapse = ", ")
  )
}


# ----------------------------- #
# 4. Create HLA-combination groups
# ----------------------------- #

hla_condition_levels <- c(
  "A2+B7+DR15+",
  "A2+B7+DR15-",
  "A2+B7-DR15+",
  "A2+B7-DR15-",
  "A2-B7+DR15+",
  "A2-B7+DR15-",
  "A2-B7-DR15+"
)

prepared_df <- dt %>%
  mutate(
    sample = as.character(.data[[sample_col]]),
    
    T_cell_freq = to_numeric_clean(.data[[config$frequency_col]]),
    
    HLA_A2 = to_binary_hla(
      .data[[config$hla_a2_col]],
      treat_unresolved_as_positive = config$treat_unresolved_hla_as_positive
    ),
    
    HLA_B7 = to_binary_hla(
      .data[[config$hla_b7_col]],
      treat_unresolved_as_positive = config$treat_unresolved_hla_as_positive
    ),
    
    HLA_DR15 = to_binary_hla(
      .data[[config$hla_dr15_col]],
      treat_unresolved_as_positive = config$treat_unresolved_hla_as_positive
    ),
    
    cohort_label = normalise_cohort(
      .data[[config$cohort_col]],
      control_label = config$control_label,
      ms_label = config$ms_label
    ),
    
    HLA_condition_raw = paste0(
      "A2", pm(HLA_A2),
      "B7", pm(HLA_B7),
      "DR15", pm(HLA_DR15)
    ),
    
    HLA_condition_raw = ifelse(
      is.na(HLA_A2) | is.na(HLA_B7) | is.na(HLA_DR15),
      NA_character_,
      HLA_condition_raw
    ),
    
    exclusion_reason = case_when(
      is.na(sample) | sample == "" ~ "Missing sample ID",
      is.na(cohort_label) ~ "Missing or unrecognised cohort",
      is.na(T_cell_freq) ~ "Missing T_cell_freq",
      !is.na(T_cell_freq) & T_cell_freq < 0 ~ "Negative T_cell_freq",
      is.na(HLA_A2) | is.na(HLA_B7) | is.na(HLA_DR15) ~ "Missing or unrecognised HLA status",
      HLA_condition_raw == "A2-B7-DR15-" ~ "Triple-negative HLA combination excluded",
      TRUE ~ NA_character_
    ),
    
    HLA_condition = ifelse(
      HLA_condition_raw == "A2-B7-DR15-",
      NA_character_,
      HLA_condition_raw
    ),
    
    HLA_condition = factor(
      HLA_condition,
      levels = hla_condition_levels
    ),
    
    # Plotting value only. Raw T_cell_freq is retained and exported.
    T_cell_freq_log10 = log10(T_cell_freq + config$pseudocount)
  )

dot_df <- prepared_df %>%
  filter(is.na(exclusion_reason)) %>%
  mutate(dot_id = row_number())

excluded_rows <- prepared_df %>%
  filter(!is.na(exclusion_reason)) %>%
  select(
    sample,
    exclusion_reason,
    cohort_label,
    HLA_A2,
    HLA_B7,
    HLA_DR15,
    HLA_condition_raw,
    T_cell_freq
  )


# ----------------------------- #
# 5. Dot-level source-data tables
# ----------------------------- #

dot_table_pooled <- dot_df %>%
  select(
    dot_id,
    sample,
    HLA_condition,
    HLA_A2,
    HLA_B7,
    HLA_DR15,
    T_cell_freq,
    T_cell_freq_log10
  ) %>%
  arrange(
    HLA_condition,
    T_cell_freq
  )

dot_table_cohort_split <- dot_df %>%
  select(
    dot_id,
    sample,
    cohort_label,
    HLA_condition,
    HLA_A2,
    HLA_B7,
    HLA_DR15,
    T_cell_freq,
    T_cell_freq_log10
  ) %>%
  arrange(
    HLA_condition,
    cohort_label,
    T_cell_freq
  )


# ----------------------------- #
# 6. Counts and summaries
# ----------------------------- #

counts_by_HLA <- dot_df %>%
  count(
    HLA_condition,
    name = "n_dots"
  ) %>%
  arrange(HLA_condition)

counts_by_HLA_and_cohort <- dot_df %>%
  count(
    HLA_condition,
    cohort_label,
    name = "n_dots"
  ) %>%
  arrange(HLA_condition, cohort_label)

summary_by_HLA <- summarise_frequency(
  df = dot_df,
  group_cols = c("HLA_condition")
)

summary_by_HLA_and_cohort <- summarise_frequency(
  df = dot_df,
  group_cols = c("HLA_condition", "cohort_label")
)


# ----------------------------- #
# 7. Export CSV files
# ----------------------------- #

pooled_dot_path <- file.path(config$output_dir, config$pooled_dot_csv)
cohort_dot_path <- file.path(config$output_dir, config$cohort_dot_csv)
counts_hla_path <- file.path(config$output_dir, config$counts_hla_csv)
counts_hla_cohort_path <- file.path(config$output_dir, config$counts_hla_cohort_csv)
summary_hla_path <- file.path(config$output_dir, config$summary_hla_csv)
summary_hla_cohort_path <- file.path(config$output_dir, config$summary_hla_cohort_csv)
excluded_rows_path <- file.path(config$output_dir, config$excluded_rows_csv)

write_csv(dot_table_pooled, pooled_dot_path)
write_csv(dot_table_cohort_split, cohort_dot_path)
write_csv(counts_by_HLA, counts_hla_path)
write_csv(counts_by_HLA_and_cohort, counts_hla_cohort_path)
write_csv(summary_by_HLA, summary_hla_path)
write_csv(summary_by_HLA_and_cohort, summary_hla_cohort_path)
write_csv(excluded_rows, excluded_rows_path)


# ----------------------------- #
# 8. Export one Excel workbook with all tables
# ----------------------------- #

workbook_path <- file.path(config$output_dir, config$workbook_xlsx)

wb <- createWorkbook()

addWorksheet(wb, "Dot_table_pooled")
writeData(wb, "Dot_table_pooled", dot_table_pooled)

addWorksheet(wb, "Dot_table_cohort_split")
writeData(wb, "Dot_table_cohort_split", dot_table_cohort_split)

addWorksheet(wb, "Counts_by_HLA")
writeData(wb, "Counts_by_HLA", counts_by_HLA)

addWorksheet(wb, "Counts_by_HLA_cohort")
writeData(wb, "Counts_by_HLA_cohort", counts_by_HLA_and_cohort)

addWorksheet(wb, "Summary_by_HLA")
writeData(wb, "Summary_by_HLA", summary_by_HLA)

addWorksheet(wb, "Summary_by_HLA_cohort")
writeData(wb, "Summary_by_HLA_cohort", summary_by_HLA_and_cohort)

addWorksheet(wb, "Excluded_rows")
writeData(wb, "Excluded_rows", excluded_rows)

addWorksheet(wb, "Notes")
writeData(
  wb,
  "Notes",
  data.frame(
    item = c(
      "Input file",
      "Input sheet",
      "Sample column used",
      "Frequency column",
      "HLA columns",
      "Excluded group",
      "Pseudocount",
      "Pooled dot table",
      "Cohort split dot table",
      "BioRender note"
    ),
    value = c(
      config$input_file,
      config$input_sheet,
      sample_col,
      config$frequency_col,
      paste(
        config$hla_a2_col,
        config$hla_b7_col,
        config$hla_dr15_col,
        sep = "; "
      ),
      "A2-B7-DR15- triple-negative individuals excluded",
      config$pseudocount,
      "One row per dot, grouped only by HLA_condition",
      "One row per dot, grouped by HLA_condition and cohort_label",
      "The exported dot-level source-data tables were imported into BioRender.com for graphing and statistical testing."
    )
  )
)

saveWorkbook(
  wb,
  workbook_path,
  overwrite = TRUE
)


# ----------------------------- #
# 9. Save session information
# ----------------------------- #

session_info_path <- file.path(config$output_dir, config$session_info_file)

writeLines(
  c(
    "Configuration:",
    capture.output(str(config)),
    "",
    "Sample column used:",
    sample_col,
    "",
    "Number of input rows:",
    as.character(nrow(dt)),
    "",
    "Number of retained dot-level rows:",
    as.character(nrow(dot_df)),
    "",
    "Number of excluded rows:",
    as.character(nrow(excluded_rows)),
    "",
    "Counts by HLA condition:",
    capture.output(print(counts_by_HLA)),
    "",
    "Counts by HLA condition and cohort:",
    capture.output(print(counts_by_HLA_and_cohort)),
    "",
    "BioRender note:",
    "The exported dot-level source-data tables were imported into BioRender.com for graphing and statistical testing.",
    "",
    "Output files:",
    pooled_dot_path,
    cohort_dot_path,
    counts_hla_path,
    counts_hla_cohort_path,
    summary_hla_path,
    summary_hla_cohort_path,
    excluded_rows_path,
    workbook_path,
    "",
    "Session information:",
    capture.output(sessionInfo())
  ),
  session_info_path
)


# ----------------------------- #
# 10. Print checks
# ----------------------------- #

message("Done.")
message("Output folder: ", config$output_dir)
message("Excel workbook: ", workbook_path)
message("Session info: ", session_info_path)

message("\nCounts by HLA condition:")
print(counts_by_HLA)

message("\nCounts by HLA condition and cohort:")
print(counts_by_HLA_and_cohort)

message("\nFirst rows of pooled dot table:")
print(head(dot_table_pooled, 20))

message("\nFirst rows of cohort-split dot table:")
print(head(dot_table_cohort_split, 20))
