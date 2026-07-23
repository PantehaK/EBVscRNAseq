#!/usr/bin/env Rscript

# ==============================================================================
# Figures S1B-C: HLA A2/B7/DR15 frequency analysis
# ==============================================================================
#
# Description
#   Generates dot-level tables, descriptive summaries and non-parametric
#   statistical tests for latent-specific CD8+ T-cell frequency across
#   HLA-A2, HLA-B7 and HLA-DR15 carrier combinations.
#
# Statistical tests
#   - Kruskal-Wallis omnibus tests:
#       * pooled
#       * Control only
#       * MS only
#   - Pairwise two-sided Wilcoxon rank-sum tests between HLA combinations:
#       * pooled
#       * Control only
#       * MS only
#   - Two-sided Wilcoxon rank-sum tests comparing MS versus Control within
#     each HLA combination.
#
# Multiple testing
#   Raw P-values and Benjamini-Hochberg FDR-adjusted P-values are exported.
#
# Usage
#   From the repository root:
#
#   Rscript scripts/figure_s1bc_hla_analysis.R
#
#   Or specify the input, output directory and Excel sheet:
#
#   Rscript scripts/figure_s1bc_hla_analysis.R \
#     data/source_data_Fig.S1B-C.xlsx \
#     results/figure_s1bc_hla \
#     Full
#
# Positional arguments
#   1. Input Excel file
#      Default: data/source_data_Fig.S1B-C.xlsx
#
#   2. Output directory
#      Default: results/figure_s1bc_hla
#
#   3. Excel sheet name or number
#      Default: Full
#
# Required input columns
#   T_cell_freq
#   cohort
#   HLA-A2
#   HLA-B7
#   HLA-DR15
#
# Sample identifier column
#   One of: sample, id, ID
#
# Outputs
#   - dot-level CSV tables
#   - count and summary CSV tables
#   - statistical-result CSV tables
#   - one Excel workbook containing all tables
#   - session_info.txt
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Dependencies
# ------------------------------------------------------------------------------

required_packages <- c(
  "readxl",
  "dplyr",
  "tidyr",
  "tibble",
  "stringr",
  "openxlsx",
  "purrr"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(openxlsx)
  library(purrr)
})


# ------------------------------------------------------------------------------
# 2. Command-line arguments and paths
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

input_file <- if (length(args) >= 1) {
  args[[1]]
} else {
  file.path("data", "source_data_Fig.S1B-C.xlsx")
}

out_dir <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path("results", "figure_s1bc_hla")
}

input_sheet <- if (length(args) >= 3) {
  if (grepl("^[0-9]+$", args[[3]])) {
    as.integer(args[[3]])
  } else {
    args[[3]]
  }
} else {
  "Full"
}

if (!file.exists(input_file)) {
  stop(
    "Input file does not exist: ",
    input_file,
    "\nRun the script from the repository root or provide the input path explicitly."
  )
}

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

pseudocount <- 0.001

session_info_file <- file.path(
  out_dir,
  "session_info.txt"
)


# ------------------------------------------------------------------------------
# 3. Load data
# ------------------------------------------------------------------------------

message("Reading input file: ", input_file)
message("Using Excel sheet: ", input_sheet)

dt <- readxl::read_excel(
  path = input_file,
  sheet = input_sheet
)

# ------------------------------------------------------------------------------
# 4. Helper functions
# ------------------------------------------------------------------------------

to_binary_hla <- function(x) {
  
  if (is.logical(x)) {
    return(ifelse(is.na(x), NA_real_, as.numeric(x)))
  }
  
  if (is.numeric(x)) {
    return(ifelse(is.na(x), NA_real_, ifelse(x > 0, 1, 0)))
  }
  
  z <- tolower(trimws(as.character(x)))
  z[z %in% c("", "na", "n/a", "nan", "null")] <- NA
  
  out <- rep(NA_real_, length(z))
  
  out[z %in% c(
    "1", "yes", "y", "positive", "pos",
    "carrier", "present", "true", "+"
  )] <- 1
  
  out[z %in% c(
    "0", "no", "n", "negative", "neg",
    "non-carrier", "absent", "false", "-"
  )] <- 0
  
  # If anything unresolved but non-missing remains, assume present.
  # Remove this if your HLA columns are definitely clean 0/1.
  unresolved <- is.na(out) & !is.na(z)
  out[unresolved] <- 1
  
  out
}

pm <- function(x) {
  case_when(
    x == 1 ~ "+",
    x == 0 ~ "-",
    TRUE ~ NA_character_
  )
}

# ============================================================
# Check required columns
# ============================================================

required_cols <- c(
  "T_cell_freq",
  "cohort",
  "HLA-A2",
  "HLA-B7",
  "HLA-DR15"
)

missing_cols <- setdiff(required_cols, colnames(dt))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# Sample column can be sample, id, or ID
sample_col <- case_when(
  "sample" %in% colnames(dt) ~ "sample",
  "id" %in% colnames(dt) ~ "id",
  "ID" %in% colnames(dt) ~ "ID",
  TRUE ~ NA_character_
)

if (is.na(sample_col)) {
  stop("Could not find sample column. Expected one of: sample, id, ID")
}

# ============================================================
# Create HLA-combination groups
# ============================================================

hla_condition_levels <- c(
  "A2+B7+DR15+",
  "A2+B7+DR15-",
  "A2+B7-DR15+",
  "A2+B7-DR15-",
  "A2-B7+DR15+",
  "A2-B7+DR15-",
  "A2-B7-DR15+"
)

dot_df <- dt %>%
  mutate(
    sample = as.character(.data[[sample_col]]),
    
    T_cell_freq = suppressWarnings(
      as.numeric(
        gsub(",", "", trimws(as.character(T_cell_freq)))
      )
    ),
    
    HLA_A2 = to_binary_hla(`HLA-A2`),
    HLA_B7 = to_binary_hla(`HLA-B7`),
    HLA_DR15 = to_binary_hla(`HLA-DR15`),
    
    cohort_label = case_when(
      cohort == 1 ~ "MS",
      cohort == 0 ~ "Control",
      as.character(cohort) %in% c("MS", "ms") ~ "MS",
      as.character(cohort) %in% c("Control", "Controls", "control", "controls") ~ "Control",
      TRUE ~ NA_character_
    ),
    cohort_label = factor(cohort_label, levels = c("Control", "MS")),
    
    HLA_condition = paste0(
      "A2", pm(HLA_A2),
      "B7", pm(HLA_B7),
      "DR15", pm(HLA_DR15)
    ),
    
    HLA_condition = ifelse(
      is.na(HLA_A2) | is.na(HLA_B7) | is.na(HLA_DR15),
      NA_character_,
      HLA_condition
    ),
    
    # Exclude triple negatives
    HLA_condition = ifelse(
      HLA_condition == "A2-B7-DR15-",
      NA_character_,
      HLA_condition
    ),
    
    HLA_condition = factor(
      HLA_condition,
      levels = hla_condition_levels
    ),
    
    # This is the plotted value if using log10 axis
    T_cell_freq_log10 = log10(T_cell_freq + pseudocount)
  ) %>%
  filter(
    !is.na(sample),
    sample != "",
    !is.na(cohort_label),
    !is.na(HLA_condition),
    !is.na(T_cell_freq),
    T_cell_freq >= 0
  ) %>%
  mutate(
    dot_id = row_number()
  )

# ============================================================
# Table 1: pooled HLA-combination dot table
# One row = one dot in the pooled graph
# ============================================================

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

# ============================================================
# Table 2: cohort-split dot table
# One row = one dot in the MS vs Control split graph
# ============================================================

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

# ============================================================
# Extra useful summaries
# ============================================================

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

summary_by_HLA <- dot_df %>%
  group_by(HLA_condition) %>%
  summarise(
    n = n(),
    n_zero = sum(T_cell_freq == 0, na.rm = TRUE),
    mean = mean(T_cell_freq, na.rm = TRUE),
    sd = sd(T_cell_freq, na.rm = TRUE),
    median = median(T_cell_freq, na.rm = TRUE),
    IQR = paste0(
      round(quantile(T_cell_freq, 0.25, na.rm = TRUE), 6),
      " - ",
      round(quantile(T_cell_freq, 0.75, na.rm = TRUE), 6)
    ),
    min_max = paste0(
      round(min(T_cell_freq, na.rm = TRUE), 6),
      " - ",
      round(max(T_cell_freq, na.rm = TRUE), 6)
    ),
    .groups = "drop"
  ) %>%
  arrange(HLA_condition)

summary_by_HLA_and_cohort <- dot_df %>%
  group_by(HLA_condition, cohort_label) %>%
  summarise(
    n = n(),
    n_zero = sum(T_cell_freq == 0, na.rm = TRUE),
    mean = mean(T_cell_freq, na.rm = TRUE),
    sd = sd(T_cell_freq, na.rm = TRUE),
    median = median(T_cell_freq, na.rm = TRUE),
    IQR = paste0(
      round(quantile(T_cell_freq, 0.25, na.rm = TRUE), 6),
      " - ",
      round(quantile(T_cell_freq, 0.75, na.rm = TRUE), 6)
    ),
    min_max = paste0(
      round(min(T_cell_freq, na.rm = TRUE), 6),
      " - ",
      round(max(T_cell_freq, na.rm = TRUE), 6)
    ),
    .groups = "drop"
  ) %>%
  arrange(HLA_condition, cohort_label)


# ============================================================
# Statistical tests
# ============================================================
#
# Tests use the untransformed T_cell_freq values.
#
# 1. Kruskal-Wallis omnibus tests:
#    - pooled across all participants
#    - separately within Control and MS
#
# 2. Pairwise Wilcoxon rank-sum tests between HLA combinations:
#    - pooled across all participants
#    - separately within Control and MS
#
# 3. Wilcoxon rank-sum tests comparing MS vs Control:
#    - separately within each HLA combination
#
# Both raw P-values and Benjamini-Hochberg FDR-adjusted
# P-values are exported.
# ============================================================

format_test_p <- function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

safe_kruskal <- function(data, analysis_set, cohort_subset) {

  tmp <- data %>%
    filter(
      !is.na(HLA_condition),
      !is.na(T_cell_freq)
    ) %>%
    mutate(
      HLA_condition = droplevels(HLA_condition)
    )

  n_groups <- n_distinct(tmp$HLA_condition)

  if (nrow(tmp) < 3 || n_groups < 2) {
    return(
      tibble(
        analysis_set = analysis_set,
        cohort_subset = cohort_subset,
        test = "Kruskal-Wallis test",
        n_total = nrow(tmp),
        n_groups = n_groups,
        statistic = NA_real_,
        df = NA_real_,
        raw_p = NA_real_,
        note = "Not tested: fewer than two represented HLA groups"
      )
    )
  }

  test_result <- tryCatch(
    kruskal.test(
      T_cell_freq ~ HLA_condition,
      data = tmp
    ),
    error = function(e) NULL
  )

  if (is.null(test_result)) {
    return(
      tibble(
        analysis_set = analysis_set,
        cohort_subset = cohort_subset,
        test = "Kruskal-Wallis test",
        n_total = nrow(tmp),
        n_groups = n_groups,
        statistic = NA_real_,
        df = NA_real_,
        raw_p = NA_real_,
        note = "Test failed"
      )
    )
  }

  tibble(
    analysis_set = analysis_set,
    cohort_subset = cohort_subset,
    test = "Kruskal-Wallis test",
    n_total = nrow(tmp),
    n_groups = n_groups,
    statistic = unname(test_result$statistic),
    df = unname(test_result$parameter),
    raw_p = test_result$p.value,
    note = NA_character_
  )
}

kruskal_results <- bind_rows(
  safe_kruskal(
    data = dot_df,
    analysis_set = "Pooled HLA comparison",
    cohort_subset = "All"
  ),
  safe_kruskal(
    data = dot_df %>% filter(cohort_label == "Control"),
    analysis_set = "HLA comparison within cohort",
    cohort_subset = "Control"
  ),
  safe_kruskal(
    data = dot_df %>% filter(cohort_label == "MS"),
    analysis_set = "HLA comparison within cohort",
    cohort_subset = "MS"
  )
) %>%
  mutate(
    FDR_BH = p.adjust(raw_p, method = "BH"),
    raw_p_label = format_test_p(raw_p),
    FDR_BH_label = format_test_p(FDR_BH)
  ) %>%
  select(
    analysis_set,
    cohort_subset,
    test,
    n_total,
    n_groups,
    statistic,
    df,
    raw_p,
    raw_p_label,
    FDR_BH,
    FDR_BH_label,
    note
  )

run_pairwise_HLA_wilcox <- function(
    data,
    analysis_set,
    cohort_subset
) {

  tmp <- data %>%
    filter(
      !is.na(HLA_condition),
      !is.na(T_cell_freq)
    ) %>%
    mutate(
      HLA_condition = droplevels(HLA_condition)
    )

  represented_groups <- levels(tmp$HLA_condition)

  if (length(represented_groups) < 2) {
    return(
      tibble(
        analysis_set = character(),
        cohort_subset = character(),
        group_1 = character(),
        group_2 = character(),
        n_group_1 = integer(),
        n_group_2 = integer(),
        median_group_1 = double(),
        median_group_2 = double(),
        median_difference_group_1_minus_group_2 = double(),
        statistic_W = double(),
        raw_p = double(),
        note = character()
      )
    )
  }

  group_pairs <- combn(
    represented_groups,
    2,
    simplify = FALSE
  )

  map_dfr(
    group_pairs,
    function(pair) {

      pair_df <- tmp %>%
        filter(HLA_condition %in% pair) %>%
        mutate(
          HLA_condition = factor(
            as.character(HLA_condition),
            levels = pair
          )
        )

      group_1_values <- pair_df %>%
        filter(HLA_condition == pair[1]) %>%
        pull(T_cell_freq)

      group_2_values <- pair_df %>%
        filter(HLA_condition == pair[2]) %>%
        pull(T_cell_freq)

      n_group_1 <- length(group_1_values)
      n_group_2 <- length(group_2_values)

      median_group_1 <- ifelse(
        n_group_1 > 0,
        median(group_1_values, na.rm = TRUE),
        NA_real_
      )

      median_group_2 <- ifelse(
        n_group_2 > 0,
        median(group_2_values, na.rm = TRUE),
        NA_real_
      )

      if (n_group_1 == 0 || n_group_2 == 0) {
        return(
          tibble(
            analysis_set = analysis_set,
            cohort_subset = cohort_subset,
            group_1 = pair[1],
            group_2 = pair[2],
            n_group_1 = n_group_1,
            n_group_2 = n_group_2,
            median_group_1 = median_group_1,
            median_group_2 = median_group_2,
            median_difference_group_1_minus_group_2 =
              median_group_1 - median_group_2,
            statistic_W = NA_real_,
            raw_p = NA_real_,
            note = "Not tested: one group has no observations"
          )
        )
      }

      test_result <- tryCatch(
        suppressWarnings(
          wilcox.test(
            T_cell_freq ~ HLA_condition,
            data = pair_df,
            exact = FALSE
          )
        ),
        error = function(e) NULL
      )

      tibble(
        analysis_set = analysis_set,
        cohort_subset = cohort_subset,
        group_1 = pair[1],
        group_2 = pair[2],
        n_group_1 = n_group_1,
        n_group_2 = n_group_2,
        median_group_1 = median_group_1,
        median_group_2 = median_group_2,
        median_difference_group_1_minus_group_2 =
          median_group_1 - median_group_2,
        statistic_W = ifelse(
          is.null(test_result),
          NA_real_,
          unname(test_result$statistic)
        ),
        raw_p = ifelse(
          is.null(test_result),
          NA_real_,
          test_result$p.value
        ),
        note = ifelse(
          is.null(test_result),
          "Test failed",
          NA_character_
        )
      )
    }
  )
}

pairwise_HLA_results <- bind_rows(
  run_pairwise_HLA_wilcox(
    data = dot_df,
    analysis_set = "Pooled HLA comparison",
    cohort_subset = "All"
  ),
  run_pairwise_HLA_wilcox(
    data = dot_df %>% filter(cohort_label == "Control"),
    analysis_set = "HLA comparison within cohort",
    cohort_subset = "Control"
  ),
  run_pairwise_HLA_wilcox(
    data = dot_df %>% filter(cohort_label == "MS"),
    analysis_set = "HLA comparison within cohort",
    cohort_subset = "MS"
  )
) %>%
  group_by(analysis_set, cohort_subset) %>%
  mutate(
    # Primary FDR correction: within each logical family of comparisons
    FDR_BH_within_set = p.adjust(raw_p, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    # Also supplied for maximum transparency across every pairwise HLA test
    FDR_BH_all_pairwise = p.adjust(raw_p, method = "BH"),
    raw_p_label = format_test_p(raw_p),
    FDR_BH_within_set_label = format_test_p(FDR_BH_within_set),
    FDR_BH_all_pairwise_label = format_test_p(FDR_BH_all_pairwise)
  ) %>%
  select(
    analysis_set,
    cohort_subset,
    group_1,
    group_2,
    n_group_1,
    n_group_2,
    median_group_1,
    median_group_2,
    median_difference_group_1_minus_group_2,
    statistic_W,
    raw_p,
    raw_p_label,
    FDR_BH_within_set,
    FDR_BH_within_set_label,
    FDR_BH_all_pairwise,
    FDR_BH_all_pairwise_label,
    note
  )

run_cohort_test_within_HLA <- function(data) {

  represented_groups <- levels(droplevels(data$HLA_condition))

  map_dfr(
    represented_groups,
    function(hla_group) {

      tmp <- data %>%
        filter(
          HLA_condition == hla_group,
          !is.na(cohort_label),
          !is.na(T_cell_freq)
        ) %>%
        mutate(
          cohort_label = droplevels(cohort_label)
        )

      control_values <- tmp %>%
        filter(cohort_label == "Control") %>%
        pull(T_cell_freq)

      ms_values <- tmp %>%
        filter(cohort_label == "MS") %>%
        pull(T_cell_freq)

      n_control <- length(control_values)
      n_MS <- length(ms_values)

      median_control <- ifelse(
        n_control > 0,
        median(control_values, na.rm = TRUE),
        NA_real_
      )

      median_MS <- ifelse(
        n_MS > 0,
        median(ms_values, na.rm = TRUE),
        NA_real_
      )

      if (n_control == 0 || n_MS == 0) {
        return(
          tibble(
            HLA_condition = hla_group,
            n_Control = n_control,
            n_MS = n_MS,
            median_Control = median_control,
            median_MS = median_MS,
            median_difference_MS_minus_Control =
              median_MS - median_control,
            statistic_W = NA_real_,
            raw_p = NA_real_,
            note = "Not tested: only one cohort represented"
          )
        )
      }

      test_result <- tryCatch(
        suppressWarnings(
          wilcox.test(
            T_cell_freq ~ cohort_label,
            data = tmp,
            exact = FALSE
          )
        ),
        error = function(e) NULL
      )

      tibble(
        HLA_condition = hla_group,
        n_Control = n_control,
        n_MS = n_MS,
        median_Control = median_control,
        median_MS = median_MS,
        median_difference_MS_minus_Control =
          median_MS - median_control,
        statistic_W = ifelse(
          is.null(test_result),
          NA_real_,
          unname(test_result$statistic)
        ),
        raw_p = ifelse(
          is.null(test_result),
          NA_real_,
          test_result$p.value
        ),
        note = ifelse(
          is.null(test_result),
          "Test failed",
          NA_character_
        )
      )
    }
  )
}

cohort_within_HLA_results <- run_cohort_test_within_HLA(dot_df) %>%
  mutate(
    FDR_BH = p.adjust(raw_p, method = "BH"),
    raw_p_label = format_test_p(raw_p),
    FDR_BH_label = format_test_p(FDR_BH)
  ) %>%
  select(
    HLA_condition,
    n_Control,
    n_MS,
    median_Control,
    median_MS,
    median_difference_MS_minus_Control,
    statistic_W,
    raw_p,
    raw_p_label,
    FDR_BH,
    FDR_BH_label,
    note
  )


# ============================================================
# Export CSV files
# ============================================================

write.csv(
  dot_table_pooled,
  file = file.path(out_dir, "dot_level_frequencies_by_HLA_condition_pooled.csv"),
  row.names = FALSE
)

write.csv(
  dot_table_cohort_split,
  file = file.path(out_dir, "dot_level_frequencies_by_HLA_condition_and_cohort.csv"),
  row.names = FALSE
)

write.csv(
  counts_by_HLA,
  file = file.path(out_dir, "counts_by_HLA_condition.csv"),
  row.names = FALSE
)

write.csv(
  counts_by_HLA_and_cohort,
  file = file.path(out_dir, "counts_by_HLA_condition_and_cohort.csv"),
  row.names = FALSE
)

write.csv(
  summary_by_HLA,
  file = file.path(out_dir, "summary_frequency_by_HLA_condition.csv"),
  row.names = FALSE
)

write.csv(
  summary_by_HLA_and_cohort,
  file = file.path(out_dir, "summary_frequency_by_HLA_condition_and_cohort.csv"),
  row.names = FALSE
)


write.csv(
  kruskal_results,
  file = file.path(out_dir, "stats_Kruskal_Wallis_omnibus.csv"),
  row.names = FALSE
)

write.csv(
  pairwise_HLA_results,
  file = file.path(out_dir, "stats_pairwise_HLA_Wilcoxon.csv"),
  row.names = FALSE
)

write.csv(
  cohort_within_HLA_results,
  file = file.path(out_dir, "stats_MS_vs_Control_within_HLA.csv"),
  row.names = FALSE
)

# ============================================================
# Export one Excel workbook with all tables
# ============================================================

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


addWorksheet(wb, "Stats_Kruskal_Wallis")
writeData(wb, "Stats_Kruskal_Wallis", kruskal_results)

addWorksheet(wb, "Stats_pairwise_HLA")
writeData(wb, "Stats_pairwise_HLA", pairwise_HLA_results)

addWorksheet(wb, "Stats_cohort_within_HLA")
writeData(wb, "Stats_cohort_within_HLA", cohort_within_HLA_results)

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
      "Excluded group",
      "Pseudocount",
      "Pooled dot table",
      "Cohort split dot table",
      "Statistical scale",
      "Omnibus tests",
      "Pairwise HLA tests",
      "Cohort tests",
      "Primary pairwise FDR",
      "Additional global pairwise FDR"
    ),
    value = c(
      input_file,
      input_sheet,
      sample_col,
      "T_cell_freq",
      "A2-B7-DR15- triple-negative individuals excluded",
      pseudocount,
      "One row per dot, grouped only by HLA_condition",
      "One row per dot, grouped by HLA_condition and cohort_label",
      "All tests use raw, untransformed T_cell_freq values",
      "Kruskal-Wallis tests compare HLA combinations pooled and within each cohort",
      "Two-sided Wilcoxon rank-sum tests compare every represented HLA pair",
      "Two-sided Wilcoxon rank-sum tests compare MS vs Control within each HLA combination",
      "FDR_BH_within_set adjusts pairwise tests separately within pooled, Control-only and MS-only sets",
      "FDR_BH_all_pairwise adjusts across every pairwise HLA test in the script"
    )
  )
)

xlsx_out <- file.path(out_dir, "HLA_A2_B7_DR15_dot_level_frequency_tables.xlsx")

saveWorkbook(
  wb,
  xlsx_out,
  overwrite = TRUE
)

# ============================================================
# Print checks
# ============================================================

capture.output(
  sessionInfo(),
  file = session_info_file
)

message("Done.")
message("Output folder: ", out_dir)
message("Excel workbook: ", xlsx_out)
message("Session information: ", session_info_file)

message("\nCounts by HLA condition:")
print(counts_by_HLA)

message("\nCounts by HLA condition and cohort:")
print(counts_by_HLA_and_cohort)

message("\nFirst rows of pooled dot table:")
print(head(dot_table_pooled, 20))

message("\nFirst rows of cohort-split dot table:")
print(head(dot_table_cohort_split, 20))

message("\nKruskal-Wallis omnibus tests:")
print(kruskal_results)

message("\nPairwise HLA Wilcoxon tests:")
print(pairwise_HLA_results)

message("\nMS vs Control within each HLA combination:")
print(cohort_within_HLA_results)

