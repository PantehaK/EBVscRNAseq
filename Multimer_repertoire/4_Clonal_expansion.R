#!/usr/bin/env Rscript

# ==============================================================================
# Expanded clonotype analysis for baseline EBV TCR data
# ==============================================================================
#
# Purpose:
#   This script calculates the proportion of expanded TRB-defined clonotypes in
#   either:
#     1. virus-negative/background T cells, or
#     2. multimer-specific T cells grouped by tetramer.
#
# Notes:
#   - Clonotypes are defined using TRB CDR3 amino acid + TRBV + TRBJ.
#   - Expanded clonotypes are defined as clone_size > 2 by default.
#   - expanded_fraction = expanded clonotypes / total clonotypes.
#   - expanded_cell_fraction = cells in expanded clonotypes / total cells.
#
# ==============================================================================


# ----------------------------- #
# 0. Load packages
# ----------------------------- #

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(readxl)
  library(stringr)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  dataset_id = "baseline_ebv_expanded_clonotypes",
  
  input_tcr_file = "path/to/input/Baseline_EBV_TCR_module_score.csv",
  
  output_dir = "path/to/output/Publication_data/EBV_baseline/Expanded_clonotypes",
  
  clone_sizes_csv = "path/to/output/Publication_data/EBV_baseline/Expanded_clonotypes/TRB_clone_sizes.csv",
  expanded_summary_csv = "path/to/output/Publication_data/EBV_baseline/Expanded_clonotypes/expanded_clonotype_summary.csv",
  cohort_summary_csv = "path/to/output/Publication_data/EBV_baseline/Expanded_clonotypes/expanded_clonotype_cohort_summary.csv",
  stats_csv = "path/to/output/Publication_data/EBV_baseline/Expanded_clonotypes/expanded_clonotype_MS_vs_Control_stats.csv",
  session_info_file = "path/to/output/Publication_data/EBV_baseline/Expanded_clonotypes/sessionInfo_expanded_clonotypes.txt",
  
  # Choose one:
  #   "virus_negative"      = matches your current script
  #   "multimer_specific"   = groups by tetramer
  analysis_mode = "virus_negative",
  
  # Used only when analysis_mode = "virus_negative"
  virus_col = "virus",
  virus_keep = "negative",
  
  # Used only when analysis_mode = "multimer_specific"
  tetramer_col = "tetramer",
  tetramers_keep = c(
    "GLCT", "CLGG", "YLQQ", "FLRG", "FLYA", "QAKW",
    "YNLR*", "RPPI", "RPQK*", "LLDF", "RAKF"
  ),
  
  # Metadata columns.
  sample_id_col = "id",
  cohort_col = "cohort",
  
  # TRB clonotype definition.
  trb_cdr3_col = "TRB_cdr3",
  trb_v_col = "TRB_v_gene",
  trb_j_col = "TRB_j_gene",
  
  # Expansion threshold.
  # clone_size > expanded_clone_min_cells is considered expanded.
  expanded_clone_min_cells = 2,
  
  # Optional statistics.
  run_wilcoxon = TRUE,
  cohort_levels = c("Control", "MS"),
  
  verbose = TRUE
)


# ----------------------------- #
# 2. Helper functions
# ----------------------------- #

create_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}


create_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}


setup_dirs <- function(cfg) {
  create_dir(cfg$output_dir)
  create_parent_dir(cfg$clone_sizes_csv)
  create_parent_dir(cfg$expanded_summary_csv)
  create_parent_dir(cfg$cohort_summary_csv)
  create_parent_dir(cfg$stats_csv)
  create_parent_dir(cfg$session_info_file)
}


read_tcr_table <- function(path) {
  if (!file.exists(path)) {
    stop("Input TCR file does not exist: ", path)
  }
  
  ext <- tools::file_ext(path) |> tolower()
  
  if (ext %in% c("xlsx", "xls")) {
    read_xlsx(path)
  } else if (ext %in% c("csv", "txt")) {
    read_csv(path, show_col_types = FALSE)
  } else {
    stop("Unsupported input file type: ", ext)
  }
}


check_required_columns <- function(df, cols) {
  missing_cols <- setdiff(cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "Input table is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
}


filter_valid_tcrs <- function(df, cfg) {
  required_cols <- c(
    cfg$sample_id_col,
    cfg$cohort_col,
    cfg$trb_cdr3_col,
    cfg$trb_v_col,
    cfg$trb_j_col
  )
  
  if (cfg$analysis_mode == "virus_negative") {
    required_cols <- c(required_cols, cfg$virus_col)
  }
  
  if (cfg$analysis_mode == "multimer_specific") {
    required_cols <- c(required_cols, cfg$tetramer_col)
  }
  
  check_required_columns(df, required_cols)
  
  df_filtered <- df %>%
    filter(
      !is.na(.data[[cfg$sample_id_col]]),
      !is.na(.data[[cfg$cohort_col]]),
      !is.na(.data[[cfg$trb_cdr3_col]]),
      .data[[cfg$trb_cdr3_col]] != "",
      !is.na(.data[[cfg$trb_v_col]]),
      .data[[cfg$trb_v_col]] != "",
      !is.na(.data[[cfg$trb_j_col]]),
      .data[[cfg$trb_j_col]] != ""
    )
  
  if (cfg$analysis_mode == "virus_negative") {
    df_filtered <- df_filtered %>%
      filter(.data[[cfg$virus_col]] == cfg$virus_keep) %>%
      mutate(analysis_group = cfg$virus_keep)
  }
  
  if (cfg$analysis_mode == "multimer_specific") {
    df_filtered <- df_filtered %>%
      filter(.data[[cfg$tetramer_col]] %in% cfg$tetramers_keep) %>%
      mutate(analysis_group = .data[[cfg$tetramer_col]])
  }
  
  df_filtered %>%
    mutate(
      clone_key_trb = paste(
        .data[[cfg$trb_cdr3_col]],
        .data[[cfg$trb_v_col]],
        .data[[cfg$trb_j_col]],
        sep = "|"
      )
    )
}


calculate_clone_sizes <- function(df_valid, cfg) {
  df_valid %>%
    count(
      id = .data[[cfg$sample_id_col]],
      cohort = .data[[cfg$cohort_col]],
      analysis_group,
      clone_key_trb,
      name = "clone_size"
    )
}


summarise_expanded_clonotypes <- function(clone_sizes, cfg) {
  clone_sizes %>%
    mutate(
      is_expanded = clone_size > cfg$expanded_clone_min_cells
    ) %>%
    group_by(id, cohort, analysis_group) %>%
    summarise(
      total_cells = sum(clone_size),
      total_clonotypes = n(),
      expanded_clonotypes = sum(is_expanded),
      expanded_cells = sum(clone_size[is_expanded]),
      expanded_fraction = expanded_clonotypes / total_clonotypes,
      expanded_cell_fraction = expanded_cells / total_cells,
      .groups = "drop"
    )
}


summarise_by_cohort <- function(expanded_df) {
  expanded_df %>%
    group_by(analysis_group, cohort) %>%
    summarise(
      n_samples = n_distinct(id),
      median_total_cells = median(total_cells, na.rm = TRUE),
      median_total_clonotypes = median(total_clonotypes, na.rm = TRUE),
      median_expanded_clonotypes = median(expanded_clonotypes, na.rm = TRUE),
      median_expanded_fraction = median(expanded_fraction, na.rm = TRUE),
      mean_expanded_fraction = mean(expanded_fraction, na.rm = TRUE),
      sd_expanded_fraction = sd(expanded_fraction, na.rm = TRUE),
      median_expanded_cell_fraction = median(expanded_cell_fraction, na.rm = TRUE),
      mean_expanded_cell_fraction = mean(expanded_cell_fraction, na.rm = TRUE),
      sd_expanded_cell_fraction = sd(expanded_cell_fraction, na.rm = TRUE),
      .groups = "drop"
    )
}


run_wilcoxon_stats <- function(expanded_df, cfg) {
  if (!isTRUE(cfg$run_wilcoxon)) {
    return(tibble())
  }
  
  dat <- expanded_df %>%
    filter(cohort %in% cfg$cohort_levels) %>%
    mutate(cohort = factor(cohort, levels = cfg$cohort_levels))
  
  metrics <- c("expanded_fraction", "expanded_cell_fraction")
  
  stats_df <- dat %>%
    group_by(analysis_group) %>%
    group_modify(~ {
      group_dat <- .x
      
      bind_rows(lapply(metrics, function(metric) {
        n_control <- sum(group_dat$cohort == cfg$cohort_levels[1])
        n_ms <- sum(group_dat$cohort == cfg$cohort_levels[2])
        
        median_control <- median(
          group_dat[[metric]][group_dat$cohort == cfg$cohort_levels[1]],
          na.rm = TRUE
        )
        
        median_ms <- median(
          group_dat[[metric]][group_dat$cohort == cfg$cohort_levels[2]],
          na.rm = TRUE
        )
        
        if (n_control < 2 || n_ms < 2) {
          return(tibble(
            metric = metric,
            n_control = n_control,
            n_ms = n_ms,
            wilcox_p = NA_real_,
            median_control = median_control,
            median_ms = median_ms,
            median_difference_ms_minus_control = median_ms - median_control
          ))
        }
        
        test <- wilcox.test(
          group_dat[[metric]] ~ group_dat$cohort,
          exact = FALSE
        )
        
        tibble(
          metric = metric,
          n_control = n_control,
          n_ms = n_ms,
          wilcox_p = test$p.value,
          median_control = median_control,
          median_ms = median_ms,
          median_difference_ms_minus_control = median_ms - median_control
        )
      }))
    }) %>%
    ungroup() %>%
    mutate(
      p_adj_BH = p.adjust(wilcox_p, method = "BH")
    )
  
  stats_df
}


# ----------------------------- #
# 3. Run pipeline
# ----------------------------- #

setup_dirs(config)

message("Loading baseline EBV TCR table...")
df <- read_tcr_table(config$input_tcr_file)

message("Filtering valid TCRs...")
df_valid <- filter_valid_tcrs(df, config)

message("Calculating TRB clone sizes...")
clone_sizes <- calculate_clone_sizes(df_valid, config)

write.csv(
  clone_sizes,
  config$clone_sizes_csv,
  row.names = FALSE
)

message("Summarising expanded clonotypes...")
expanded_df <- summarise_expanded_clonotypes(clone_sizes, config)

write.csv(
  expanded_df,
  config$expanded_summary_csv,
  row.names = FALSE
)

message("Summarising expanded clonotypes by cohort...")
cohort_summary <- summarise_by_cohort(expanded_df)

write.csv(
  cohort_summary,
  config$cohort_summary_csv,
  row.names = FALSE
)

message("Running optional MS vs Control statistics...")
stats_df <- run_wilcoxon_stats(expanded_df, config)

write.csv(
  stats_df,
  config$stats_csv,
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  config$session_info_file
)

message("\nExpanded clonotype analysis complete.")
message("Analysis mode: ", config$analysis_mode)
message("Saved clone sizes to: ", config$clone_sizes_csv)
message("Saved expanded clonotype summary to: ", config$expanded_summary_csv)
message("Saved cohort summary to: ", config$cohort_summary_csv)
message("Saved statistics to: ", config$stats_csv)
message("Saved session info to: ", config$session_info_file)