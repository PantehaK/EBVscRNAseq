#!/usr/bin/env Rscript

# ==============================================================================
# D50 of EBV tetramer-specific TRB clonotypes from baseline RDS
# ==============================================================================
#
# Purpose:
#   This script calculates D50 for baseline EBV tetramer-enriched TRB clonotypes.
#
#   Clone definition:
#     TRB_cdr3 + TRB_v_gene + TRB_j_gene
#
#   Filtering:
#     - Only cells from batches beginning with GEMEBV are included.
#     - Only tetramers listed in config$tet_order are included.
#     - Only cells with complete TRB_cdr3, TRB_v_gene, and TRB_j_gene are used.
#
#   D50 calculation:
#     For each id x tetramer specificity:
#       1. Count cells per TRB clonotype.
#       2. Optionally keep only clonotypes above a minimum clone-size cutoff.
#       3. Sort clonotypes from largest to smallest.
#       4. D50_count = number of top clonotypes needed to account for >=50%
#          of cells in the D50-eligible clonotype pool.
#       5. D50_fraction = D50_count / number of D50-eligible clonotypes.
#       6. D50_percent = 100 * D50_fraction.
#
#   Important:
#     - Primary D50 is calculated at id x tetramer specificity.
#     - Primary Wilcoxon + BH statistics are run on the id x tetramer D50 table.
#     - Antigen-level D50 is secondary and calculated after collapsing clone sizes
#       by id x antigen.
#
# Notes:
#   - Paths are intentionally generic for GitHub.
#   - Update the config block before running.
#   - The uploaded analysis used min_clone_size_for_D50 = 0, which means all
#     clonotypes with clone_size > 0 are included in D50. Set this to 3 if you
#     want to restrict D50 to clonotypes with at least 3 cells.
#
# ==============================================================================


# ----------------------------- #
# 0. Load packages
# ----------------------------- #

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readr)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  
  # Input RDS
  input_rds = "path/to/Publication_RDS/Baseline_EBV_integrated_seurat.rds",
  
  # Output directory
  out_dir = "path/to/Publication_data/EBV_baseline/Clonal_expansion",
  
  # Metadata columns
  id_col = "id",
  cohort_col = "cohort",
  batch_col = "batch",
  tetramer_col = "tetramer",
  antigen_col = "antigen",
  
  # TRB columns used to define clone identity
  trb_cdr3_col = "TRB_cdr3",
  trb_v_col = "TRB_v_gene",
  trb_j_col = "TRB_j_gene",
  
  # Batch filter
  gemebv_batch_prefix = "GEMEBV",
  
  # Tetramers to include
  tet_order = c(
    "GLCT", "RAKF",
    "RPPI", "FLRG", "QAKW",
    "AVFD*", "LLDF", "YLQQ",
    "CLGG", "FLYA",
    "RPQK*", "YNLR*"
  ),
  
  # D50 settings
  # Set to 0 to include all clonotypes with clone_size > 0.
  # Set to 3 to calculate D50 using only clonotypes with at least 3 cells.
  min_clone_size_for_D50 = 0,
  
  # D50 metric to use for MS-vs-Control statistics.
  # Options:
  #   "D50_count"
  #   "D50_fraction"
  #   "D50_percent"
  d50_metric_for_stats = "D50_percent",
  
  # Cohort labels
  control_label = "Control",
  ms_label = "MS",
  
  # Statistics settings
  min_n_per_cohort_for_test = 2,
  p_adjust_method = "BH",
  fold_change_offset = 1e-9,
  
  # Output files
  clone_size_csv = "GEMEBV_TRB_clone_sizes_by_id_tetramer.csv",
  clone_size_d50_eligible_csv = "GEMEBV_TRB_clone_sizes_D50_eligible_by_id_tetramer.csv",
  id_tetramer_d50_csv = "GEMEBV_TRB_D50_by_id_tetramer.csv",
  id_antigen_d50_csv = "GEMEBV_TRB_D50_by_id_antigen.csv",
  cohort_tetramer_summary_csv = "GEMEBV_TRB_D50_summary_by_cohort_tetramer.csv",
  
  tetramer_stats_input_csv = "GEMEBV_TRB_D50_input_by_tetramer.csv",
  tetramer_stats_csv = "GEMEBV_TRB_D50_stats_by_tetramer_BH.csv",
  
  antigen_stats_input_csv = "GEMEBV_TRB_D50_input_by_antigen.csv",
  antigen_stats_csv = "GEMEBV_TRB_D50_stats_by_antigen_BH.csv",
  
  session_info_file = "sessionInfo_GEMEBV_TRB_D50.txt"
)


# ----------------------------- #
# 2. Helper functions
# ----------------------------- #

create_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}


load_seurat_object <- function(path) {
  
  if (!file.exists(path)) {
    stop("Input RDS file does not exist: ", path)
  }
  
  obj <- readRDS(path)
  
  if (!inherits(obj, "Seurat")) {
    stop("Input RDS does not contain a Seurat object: ", path)
  }
  
  obj
}


check_required_columns <- function(df, required_cols, object_name = "dataframe") {
  
  missing_cols <- setdiff(required_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ",
      object_name,
      ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  invisible(TRUE)
}


clean_chr <- function(x) {
  str_trim(as.character(x))
}


clean_cohort <- function(x, control_label = "Control", ms_label = "MS") {
  
  x_raw <- as.character(x)
  x_clean <- str_to_lower(str_trim(x_raw))
  
  case_when(
    x_clean %in% c(
      "control", "controls", "healthy control", "healthy controls",
      "hc", "nms", "non-ms", "non_ms"
    ) ~ control_label,
    
    x_clean %in% c(
      "ms", "pwms", "multiple sclerosis", "rrms", "ppms", "spms"
    ) ~ ms_label,
    
    TRUE ~ x_raw
  )
}


clean_tetramer <- function(x) {
  
  # Handles either "GLCT" or "A02_GLCT" style tetramer names.
  x %>%
    clean_chr() %>%
    str_remove("^[^_]+_")
}


calc_D50 <- function(clone_sizes_vec) {
  
  clone_sizes_vec <- clone_sizes_vec[!is.na(clone_sizes_vec)]
  clone_sizes_vec <- clone_sizes_vec[clone_sizes_vec > 0]
  
  if (length(clone_sizes_vec) == 0) {
    return(tibble(
      eligible_cells_for_D50 = 0,
      eligible_clonotypes_for_D50 = 0,
      D50_count = NA_real_,
      D50_fraction = NA_real_,
      D50_percent = NA_real_
    ))
  }
  
  clone_sizes_sorted <- sort(clone_sizes_vec, decreasing = TRUE)
  cumulative_fraction <- cumsum(clone_sizes_sorted) / sum(clone_sizes_sorted)
  
  D50_count <- which(cumulative_fraction >= 0.5)[1]
  n_clonotypes <- length(clone_sizes_sorted)
  
  tibble(
    eligible_cells_for_D50 = sum(clone_sizes_sorted),
    eligible_clonotypes_for_D50 = n_clonotypes,
    D50_count = D50_count,
    D50_fraction = D50_count / n_clonotypes,
    D50_percent = 100 * D50_fraction
  )
}


safe_wilcox_p <- function(
    values,
    groups,
    control_label = "Control",
    ms_label = "MS",
    min_n_per_cohort = 2
) {
  
  test_df <- tibble(
    value = suppressWarnings(as.numeric(values)),
    group = as.character(groups)
  ) %>%
    filter(
      !is.na(value),
      group %in% c(control_label, ms_label)
    )
  
  n_control <- sum(test_df$group == control_label)
  n_ms <- sum(test_df$group == ms_label)
  
  if (n_control < min_n_per_cohort || n_ms < min_n_per_cohort) {
    return(NA_real_)
  }
  
  if (length(unique(test_df$value)) < 2) {
    return(1)
  }
  
  suppressWarnings(
    wilcox.test(
      value ~ group,
      data = test_df,
      exact = FALSE
    )$p.value
  )
}


format_p <- function(p) {
  
  case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}


run_mw_stats <- function(
    df,
    group_cols,
    metric_col = "D50",
    control_label = "Control",
    ms_label = "MS",
    min_n_per_cohort = 2,
    p_adjust_method = "BH",
    fold_change_offset = 1e-9
) {
  
  check_required_columns(
    df,
    required_cols = c(group_cols, "Cohort", metric_col),
    object_name = "stats input"
  )
  
  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_Control = sum(Cohort == control_label),
      n_MS = sum(Cohort == ms_label),
      
      p_value = safe_wilcox_p(
        values = .data[[metric_col]],
        groups = Cohort,
        control_label = control_label,
        ms_label = ms_label,
        min_n_per_cohort = min_n_per_cohort
      ),
      
      median_Control = median(.data[[metric_col]][Cohort == control_label], na.rm = TRUE),
      median_MS = median(.data[[metric_col]][Cohort == ms_label], na.rm = TRUE),
      
      diff_median = median_MS - median_Control,
      fold_change = median_MS / (median_Control + fold_change_offset),
      
      .groups = "drop"
    ) %>%
    mutate(
      p_adj = p.adjust(p_value, method = p_adjust_method),
      p_value_label = format_p(p_value),
      p_adj_label = format_p(p_adj)
    ) %>%
    arrange(p_adj, p_value)
}


# ----------------------------- #
# 3. Load RDS and extract metadata
# ----------------------------- #

create_dir(config$out_dir)

message("Loading baseline EBV object:")
message(config$input_rds)

merged <- load_seurat_object(config$input_rds)

meta <- merged@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell_barcode")

required_cols <- c(
  config$id_col,
  config$cohort_col,
  config$batch_col,
  config$tetramer_col,
  config$antigen_col,
  config$trb_cdr3_col,
  config$trb_v_col,
  config$trb_j_col
)

check_required_columns(
  df = meta,
  required_cols = required_cols,
  object_name = "baseline EBV metadata"
)


# ----------------------------- #
# 4. Clean metadata and build TRB-only clone key
# ----------------------------- #

df <- meta %>%
  mutate(
    id = clean_chr(.data[[config$id_col]]),
    cohort = clean_cohort(
      .data[[config$cohort_col]],
      control_label = config$control_label,
      ms_label = config$ms_label
    ),
    batch = clean_chr(.data[[config$batch_col]]),
    antigen = clean_chr(.data[[config$antigen_col]]),
    tetramer_clean = clean_tetramer(.data[[config$tetramer_col]]),
    
    TRB_cdr3 = clean_chr(.data[[config$trb_cdr3_col]]),
    TRB_v_gene = clean_chr(.data[[config$trb_v_col]]),
    TRB_j_gene = clean_chr(.data[[config$trb_j_col]])
  ) %>%
  filter(
    str_starts(batch, config$gemebv_batch_prefix),
    tetramer_clean %in% config$tet_order,
    
    !is.na(id), id != "",
    !is.na(cohort), cohort != "",
    !is.na(antigen), antigen != "",
    
    !is.na(TRB_cdr3), TRB_cdr3 != "",
    !is.na(TRB_v_gene), TRB_v_gene != "",
    !is.na(TRB_j_gene), TRB_j_gene != ""
  ) %>%
  mutate(
    clone_key = paste(
      TRB_cdr3,
      TRB_v_gene,
      TRB_j_gene,
      sep = "|"
    ),
    clone_definition = "TRB_cdr3|TRB_v_gene|TRB_j_gene"
  )

message("Cells retained after GEMEBV/tetramer/TRB filtering: ", nrow(df))

if (nrow(df) == 0) {
  stop("No cells remained after filtering. Check batch, tetramer, id, and TRB metadata columns.")
}

message("Tetramer counts after filtering:")
print(table(df$tetramer_clean, useNA = "ifany"))


# ----------------------------- #
# 5. Clone sizes within id x tetramer specificity
# ----------------------------- #

clone_sizes <- df %>%
  count(
    id,
    cohort,
    tetramer_clean,
    antigen,
    clone_key,
    clone_definition,
    name = "clone_size"
  )

# Keep total cells/clonotypes before applying the optional D50 clone-size cutoff.
id_tet_totals_all <- clone_sizes %>%
  group_by(id, cohort, tetramer_clean, antigen) %>%
  summarise(
    total_tetramer_specific_cells_all_clones = sum(clone_size),
    total_clonotypes_all = n(),
    .groups = "drop"
  )

clone_sizes_D50 <- clone_sizes %>%
  filter(clone_size >= config$min_clone_size_for_D50)


# ----------------------------- #
# 6. D50 per id x tetramer
# ----------------------------- #

D50_summary <- clone_sizes_D50 %>%
  group_by(id, cohort, tetramer_clean, antigen) %>%
  summarise(
    calc_D50(clone_size),
    D50_definition = paste0(
      "TRB clone; D50 calculated using clonotypes with clone_size >= ",
      config$min_clone_size_for_D50,
      " within id x tetramer"
    ),
    .groups = "drop"
  ) %>%
  right_join(
    id_tet_totals_all,
    by = c("id", "cohort", "tetramer_clean", "antigen")
  ) %>%
  mutate(
    eligible_cells_for_D50 = replace_na(eligible_cells_for_D50, 0),
    eligible_clonotypes_for_D50 = replace_na(eligible_clonotypes_for_D50, 0),
    
    D50_definition = if_else(
      is.na(D50_definition),
      paste0(
        "D50 not calculated: no TRB clonotypes with clone_size >= ",
        config$min_clone_size_for_D50,
        " within id x tetramer"
      ),
      D50_definition
    ),
    
    tetramer_clean = factor(tetramer_clean, levels = config$tet_order)
  ) %>%
  arrange(id, tetramer_clean)


# ----------------------------- #
# 7. D50 per id x antigen
# ----------------------------- #
#
# This is a secondary table. It collapses clone sizes across tetramers within
# antigen and then calculates D50 within id x antigen.
#
# ----------------------------- #

clone_sizes_antigen <- df %>%
  count(
    id,
    cohort,
    antigen,
    clone_key,
    name = "clone_size"
  )

id_antigen_totals_all <- clone_sizes_antigen %>%
  group_by(id, cohort, antigen) %>%
  summarise(
    total_antigen_specific_cells_all_clones = sum(clone_size),
    total_clonotypes_all = n(),
    .groups = "drop"
  )

clone_sizes_antigen_D50 <- clone_sizes_antigen %>%
  filter(clone_size >= config$min_clone_size_for_D50)

D50_summary_by_antigen <- clone_sizes_antigen_D50 %>%
  group_by(id, cohort, antigen) %>%
  summarise(
    calc_D50(clone_size),
    D50_definition = paste0(
      "TRB clone; D50 calculated using clonotypes with clone_size >= ",
      config$min_clone_size_for_D50,
      " within id x antigen"
    ),
    .groups = "drop"
  ) %>%
  right_join(
    id_antigen_totals_all,
    by = c("id", "cohort", "antigen")
  ) %>%
  mutate(
    eligible_cells_for_D50 = replace_na(eligible_cells_for_D50, 0),
    eligible_clonotypes_for_D50 = replace_na(eligible_clonotypes_for_D50, 0),
    
    D50_definition = if_else(
      is.na(D50_definition),
      paste0(
        "D50 not calculated: no TRB clonotypes with clone_size >= ",
        config$min_clone_size_for_D50,
        " within id x antigen"
      ),
      D50_definition
    )
  ) %>%
  arrange(id, antigen)


# ----------------------------- #
# 8. Quick cohort/tetramer summary
# ----------------------------- #

quick_summary <- D50_summary %>%
  group_by(cohort, tetramer_clean, antigen) %>%
  summarise(
    n_ids = n(),
    n_ids_with_D50 = sum(!is.na(D50_percent)),
    median_D50_count = median(D50_count, na.rm = TRUE),
    median_D50_fraction = median(D50_fraction, na.rm = TRUE),
    median_D50_percent = median(D50_percent, na.rm = TRUE),
    median_eligible_clonotypes_for_D50 =
      median(eligible_clonotypes_for_D50, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(tetramer_clean, cohort)


# ----------------------------- #
# 9. Wilcoxon + BH tests
# ----------------------------- #

if (!config$d50_metric_for_stats %in% colnames(D50_summary)) {
  stop(
    "config$d50_metric_for_stats is not in D50_summary: ",
    config$d50_metric_for_stats
  )
}

# Primary stats table:
# This is intentionally derived from D50_summary, which is id x tetramer.
df_fc_tetramer <- D50_summary %>%
  transmute(
    ID = id,
    Cohort = cohort,
    Tetramer = as.character(tetramer_clean),
    Antigen = antigen,
    
    D50 = .data[[config$d50_metric_for_stats]],
    D50_count,
    D50_fraction,
    D50_percent,
    
    eligible_cells_for_D50,
    eligible_clonotypes_for_D50,
    total_tetramer_specific_cells_all_clones,
    total_clonotypes_all,
    D50_definition
  ) %>%
  filter(
    !is.na(Tetramer),
    !is.na(Cohort),
    !is.na(D50),
    Cohort %in% c(config$control_label, config$ms_label)
  )

mw_results_tetramer <- run_mw_stats(
  df = df_fc_tetramer,
  group_cols = c("Tetramer", "Antigen"),
  metric_col = "D50",
  control_label = config$control_label,
  ms_label = config$ms_label,
  min_n_per_cohort = config$min_n_per_cohort_for_test,
  p_adjust_method = config$p_adjust_method,
  fold_change_offset = config$fold_change_offset
)

# Secondary stats table:
# This is derived from D50_summary_by_antigen.
df_fc_antigen <- D50_summary_by_antigen %>%
  transmute(
    ID = id,
    Cohort = cohort,
    Antigen = antigen,
    
    D50 = .data[[config$d50_metric_for_stats]],
    D50_count,
    D50_fraction,
    D50_percent,
    
    eligible_cells_for_D50,
    eligible_clonotypes_for_D50,
    total_antigen_specific_cells_all_clones,
    total_clonotypes_all,
    D50_definition
  ) %>%
  filter(
    !is.na(Antigen),
    !is.na(Cohort),
    !is.na(D50),
    Cohort %in% c(config$control_label, config$ms_label)
  )

mw_results_antigen <- run_mw_stats(
  df = df_fc_antigen,
  group_cols = "Antigen",
  metric_col = "D50",
  control_label = config$control_label,
  ms_label = config$ms_label,
  min_n_per_cohort = config$min_n_per_cohort_for_test,
  p_adjust_method = config$p_adjust_method,
  fold_change_offset = config$fold_change_offset
)


# ----------------------------- #
# 10. Save outputs
# ----------------------------- #

clone_size_path <- file.path(config$out_dir, config$clone_size_csv)
clone_size_d50_eligible_path <- file.path(config$out_dir, config$clone_size_d50_eligible_csv)
id_tetramer_d50_path <- file.path(config$out_dir, config$id_tetramer_d50_csv)
id_antigen_d50_path <- file.path(config$out_dir, config$id_antigen_d50_csv)
quick_summary_path <- file.path(config$out_dir, config$cohort_tetramer_summary_csv)

tetramer_stats_input_path <- file.path(config$out_dir, config$tetramer_stats_input_csv)
tetramer_stats_path <- file.path(config$out_dir, config$tetramer_stats_csv)

antigen_stats_input_path <- file.path(config$out_dir, config$antigen_stats_input_csv)
antigen_stats_path <- file.path(config$out_dir, config$antigen_stats_csv)

session_info_path <- file.path(config$out_dir, config$session_info_file)

write_csv(clone_sizes, clone_size_path)
write_csv(clone_sizes_D50, clone_size_d50_eligible_path)
write_csv(D50_summary, id_tetramer_d50_path)
write_csv(D50_summary_by_antigen, id_antigen_d50_path)
write_csv(quick_summary, quick_summary_path)

write_csv(df_fc_tetramer, tetramer_stats_input_path)
write_csv(mw_results_tetramer, tetramer_stats_path)

write_csv(df_fc_antigen, antigen_stats_input_path)
write_csv(mw_results_antigen, antigen_stats_path)

writeLines(
  c(
    "Configuration:",
    capture.output(str(config)),
    "",
    "Cells retained after GEMEBV/tetramer/TRB filtering:",
    as.character(nrow(df)),
    "",
    "Clone definition:",
    unique(df$clone_definition),
    "",
    "D50 metric used for statistics:",
    config$d50_metric_for_stats,
    "",
    "Primary stats dataset:",
    "df_fc_tetramer, derived from D50_summary at id x tetramer level.",
    "",
    "Secondary stats dataset:",
    "df_fc_antigen, derived from D50_summary_by_antigen at id x antigen level.",
    "",
    "Output files:",
    clone_size_path,
    clone_size_d50_eligible_path,
    id_tetramer_d50_path,
    id_antigen_d50_path,
    quick_summary_path,
    tetramer_stats_input_path,
    tetramer_stats_path,
    antigen_stats_input_path,
    antigen_stats_path,
    "",
    "Tetramer-level Wilcoxon + BH results:",
    capture.output(print(mw_results_tetramer)),
    "",
    "Antigen-level Wilcoxon + BH results:",
    capture.output(print(mw_results_antigen)),
    "",
    "Session information:",
    capture.output(sessionInfo())
  ),
  session_info_path
)


# ----------------------------- #
# 11. Print quick outputs
# ----------------------------- #

message("\nD50 workflow complete.")
message("Clone-size output: ", clone_size_path)
message("D50-eligible clone-size output: ", clone_size_d50_eligible_path)
message("id x tetramer D50 output: ", id_tetramer_d50_path)
message("id x antigen D50 output: ", id_antigen_d50_path)
message("Tetramer D50 stats output: ", tetramer_stats_path)
message("Antigen D50 stats output: ", antigen_stats_path)
message("Session info: ", session_info_path)

message("\nD50 metric used for statistics:")
message(config$d50_metric_for_stats)

message("\nQuick cohort/tetramer summary:")
print(quick_summary, n = Inf)

message("\nTetramer-level Wilcoxon + BH results:")
print(mw_results_tetramer, n = Inf)

message("\nAntigen-level Wilcoxon + BH results:")
print(mw_results_antigen, n = Inf)
