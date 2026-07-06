#!/usr/bin/env Rscript

# ==============================================================================
# Expanded EBV tetramer-specific TRB clonotype proportions from baseline RDS
# ==============================================================================
#
# Purpose:
#   This script calculates the proportion of expanded TRB-defined clonotypes in
#   baseline EBV tetramer-enriched cells.
#
#   Clone definition:
#     TRB_cdr3 + TRB_v_gene + TRB_j_gene
#
#   Expansion definition:
#     Within each id x tetramer specificity, a clone is called expanded when:
#
#       clone_size / total tetramer-specific cells > 0.10
#
#   Important:
#     - Expansion is calculated within id x tetramer specificity.
#     - Only cells from batches beginning with GEMEBV are included.
#     - No minimum clone-size cutoff is applied.
#
#   Outputs:
#     1. Clone-level expansion table.
#     2. id x tetramer expansion summary.
#     3. id x antigen expansion summary.
#     4. Wilcoxon + BH/FDR statistics by tetramer.
#     5. Wilcoxon + BH/FDR statistics by antigen.
#     6. Session information.
#
# Notes:
#   - Paths are intentionally generic for GitHub.
#   - Update the config block before running.
#   - The primary statistics table is the tetramer-level table because the
#     expansion metric is defined at id x tetramer specificity.
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
  
  # Expansion threshold
  # Expanded = clone proportion within id x tetramer > this value
  expanded_clone_prop_cutoff = 0.10,
  
  # Cohort labels
  control_label = "Control",
  ms_label = "MS",
  
  # Statistics settings
  expansion_metric_for_stats = "proportion_expanded_clonotypes",
  min_n_per_cohort_for_test = 2,
  p_adjust_method = "BH",
  fold_change_offset = 1e-9,
  
  # Output files
  clone_level_csv = "GEMEBV_TRB_clone_level_expansion_by_id_tetramer.csv",
  id_tetramer_csv = "GEMEBV_TRB_expanded_clone_proportions_by_id_tetramer.csv",
  id_antigen_csv = "GEMEBV_TRB_expanded_clone_proportions_by_id_antigen.csv",
  cohort_tetramer_summary_csv = "GEMEBV_TRB_expanded_clone_summary_by_cohort_tetramer.csv",
  
  tetramer_stats_input_csv = "GEMEBV_TRB_expanded_clone_proportion_input_by_tetramer.csv",
  tetramer_stats_csv = "GEMEBV_TRB_expanded_clone_proportion_stats_by_tetramer_BH.csv",
  
  antigen_stats_input_csv = "GEMEBV_TRB_expanded_clone_proportion_input_by_antigen.csv",
  antigen_stats_csv = "GEMEBV_TRB_expanded_clone_proportion_stats_by_antigen_BH.csv",
  
  session_info_file = "sessionInfo_GEMEBV_TRB_expanded_clone_proportions.txt"
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
    x_clean %in% c("control", "controls", "healthy control", "healthy controls", "hc", "nms", "non-ms", "non_ms") ~ control_label,
    x_clean %in% c("ms", "pwms", "multiple sclerosis", "rrms", "ppms", "spms") ~ ms_label,
    TRUE ~ x_raw
  )
}


clean_tetramer <- function(x) {
  
  # Handles either "GLCT" or "A02_GLCT" style tetramer names.
  x %>%
    clean_chr() %>%
    str_remove("^[^_]+_")
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
    metric_col = "Expansion",
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

# Denominator = all tetramer-specific cells for that id x tetramer.
id_tet_totals <- clone_sizes %>%
  group_by(id, cohort, tetramer_clean, antigen) %>%
  summarise(
    total_tetramer_specific_cells = sum(clone_size),
    total_clonotypes = n(),
    .groups = "drop"
  )

# Expansion is calculated within id x tetramer specificity.
clone_expansion <- clone_sizes %>%
  left_join(
    id_tet_totals,
    by = c("id", "cohort", "tetramer_clean", "antigen")
  ) %>%
  mutate(
    clone_proportion_within_tetramer =
      clone_size / total_tetramer_specific_cells,
    
    expanded_clone =
      clone_proportion_within_tetramer > config$expanded_clone_prop_cutoff,
    
    expansion_definition = paste0(
      "TRB clone; clone proportion within id x tetramer > ",
      config$expanded_clone_prop_cutoff,
      "; no minimum clone-size cutoff"
    )
  ) %>%
  arrange(
    id,
    tetramer_clean,
    desc(clone_size)
  )


# ----------------------------- #
# 6. id x tetramer summary
# ----------------------------- #

expanded_summary <- clone_expansion %>%
  group_by(id, cohort, tetramer_clean, antigen) %>%
  summarise(
    total_tetramer_specific_cells = first(total_tetramer_specific_cells),
    total_clonotypes = first(total_clonotypes),
    
    expanded_clonotypes = sum(expanded_clone, na.rm = TRUE),
    expanded_cells = sum(clone_size[expanded_clone], na.rm = TRUE),
    
    proportion_expanded_clonotypes =
      expanded_clonotypes / total_clonotypes,
    
    proportion_expanded_cells =
      expanded_cells / total_tetramer_specific_cells,
    
    clone_definition = first(clone_definition),
    expansion_definition = first(expansion_definition),
    
    .groups = "drop"
  ) %>%
  mutate(
    tetramer_clean = factor(tetramer_clean, levels = config$tet_order)
  ) %>%
  arrange(id, tetramer_clean)


# ----------------------------- #
# 7. id x antigen summary
# ----------------------------- #
#
# The expansion call is still made at id x tetramer first.
# This table collapses those clone-level calls across tetramers by antigen.
#
# ----------------------------- #

expanded_summary_by_antigen <- clone_expansion %>%
  group_by(id, cohort, antigen) %>%
  summarise(
    total_tetramer_specific_cells = sum(clone_size),
    total_clonotypes = n(),
    
    expanded_clonotypes = sum(expanded_clone, na.rm = TRUE),
    expanded_cells = sum(clone_size[expanded_clone], na.rm = TRUE),
    
    proportion_expanded_clonotypes =
      expanded_clonotypes / total_clonotypes,
    
    proportion_expanded_cells =
      expanded_cells / total_tetramer_specific_cells,
    
    clone_definition = first(clone_definition),
    expansion_definition = first(expansion_definition),
    
    .groups = "drop"
  ) %>%
  arrange(id, antigen)


# ----------------------------- #
# 8. Quick cohort/tetramer summary
# ----------------------------- #

quick_summary <- expanded_summary %>%
  group_by(cohort, tetramer_clean, antigen) %>%
  summarise(
    n_ids = n(),
    median_prop_expanded_clonotypes =
      median(proportion_expanded_clonotypes, na.rm = TRUE),
    median_prop_expanded_cells =
      median(proportion_expanded_cells, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(tetramer_clean, cohort)


# ----------------------------- #
# 9. Wilcoxon + BH tests
# ----------------------------- #

if (!config$expansion_metric_for_stats %in% colnames(expanded_summary)) {
  stop(
    "config$expansion_metric_for_stats is not in expanded_summary: ",
    config$expansion_metric_for_stats
  )
}

# Primary stats dataset:
# This is intentionally built from expanded_summary, which is id x tetramer.
# That preserves the correct tetramer-specific denominator.
df_fc_tetramer <- expanded_summary %>%
  transmute(
    ID = id,
    Cohort = cohort,
    Tetramer = as.character(tetramer_clean),
    Antigen = antigen,
    Expansion = .data[[config$expansion_metric_for_stats]],
    total_tetramer_specific_cells,
    total_clonotypes,
    expanded_clonotypes,
    expanded_cells,
    proportion_expanded_cells,
    clone_definition,
    expansion_definition
  ) %>%
  filter(
    !is.na(Tetramer),
    !is.na(Cohort),
    !is.na(Expansion),
    Cohort %in% c(config$control_label, config$ms_label)
  )

mw_results_tetramer <- run_mw_stats(
  df = df_fc_tetramer,
  group_cols = c("Tetramer", "Antigen"),
  metric_col = "Expansion",
  control_label = config$control_label,
  ms_label = config$ms_label,
  min_n_per_cohort = config$min_n_per_cohort_for_test,
  p_adjust_method = config$p_adjust_method,
  fold_change_offset = config$fold_change_offset
)

# Secondary stats dataset:
# This is built from expanded_summary_by_antigen after expansion has already
# been called at id x tetramer.
df_fc_antigen <- expanded_summary_by_antigen %>%
  transmute(
    ID = id,
    Cohort = cohort,
    Antigen = antigen,
    Expansion = .data[[config$expansion_metric_for_stats]],
    total_tetramer_specific_cells,
    total_clonotypes,
    expanded_clonotypes,
    expanded_cells,
    proportion_expanded_cells,
    clone_definition,
    expansion_definition
  ) %>%
  filter(
    !is.na(Antigen),
    !is.na(Cohort),
    !is.na(Expansion),
    Cohort %in% c(config$control_label, config$ms_label)
  )

mw_results_antigen <- run_mw_stats(
  df = df_fc_antigen,
  group_cols = "Antigen",
  metric_col = "Expansion",
  control_label = config$control_label,
  ms_label = config$ms_label,
  min_n_per_cohort = config$min_n_per_cohort_for_test,
  p_adjust_method = config$p_adjust_method,
  fold_change_offset = config$fold_change_offset
)


# ----------------------------- #
# 10. Save outputs
# ----------------------------- #

clone_level_path <- file.path(config$out_dir, config$clone_level_csv)
id_tetramer_path <- file.path(config$out_dir, config$id_tetramer_csv)
id_antigen_path <- file.path(config$out_dir, config$id_antigen_csv)
quick_summary_path <- file.path(config$out_dir, config$cohort_tetramer_summary_csv)

tetramer_stats_input_path <- file.path(config$out_dir, config$tetramer_stats_input_csv)
tetramer_stats_path <- file.path(config$out_dir, config$tetramer_stats_csv)

antigen_stats_input_path <- file.path(config$out_dir, config$antigen_stats_input_csv)
antigen_stats_path <- file.path(config$out_dir, config$antigen_stats_csv)

session_info_path <- file.path(config$out_dir, config$session_info_file)

write_csv(clone_expansion, clone_level_path)
write_csv(expanded_summary, id_tetramer_path)
write_csv(expanded_summary_by_antigen, id_antigen_path)
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
    "Expansion definition:",
    unique(clone_expansion$expansion_definition),
    "",
    "Primary stats dataset:",
    "df_fc_tetramer, derived from expanded_summary at id x tetramer level.",
    "",
    "Secondary stats dataset:",
    "df_fc_antigen, derived from expanded_summary_by_antigen after expansion was called at id x tetramer level.",
    "",
    "Output files:",
    clone_level_path,
    id_tetramer_path,
    id_antigen_path,
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

message("\nExpanded TRB clone proportion workflow complete.")
message("Clone-level output: ", clone_level_path)
message("id x tetramer output: ", id_tetramer_path)
message("id x antigen output: ", id_antigen_path)
message("Tetramer stats output: ", tetramer_stats_path)
message("Antigen stats output: ", antigen_stats_path)
message("Session info: ", session_info_path)

message("\nExpansion definition:")
message(unique(clone_expansion$expansion_definition))

message("\nQuick cohort/tetramer summary:")
print(quick_summary, n = Inf)

message("\nTetramer-level Wilcoxon + BH results:")
print(mw_results_tetramer, n = Inf)

message("\nAntigen-level Wilcoxon + BH results:")
print(mw_results_antigen, n = Inf)
