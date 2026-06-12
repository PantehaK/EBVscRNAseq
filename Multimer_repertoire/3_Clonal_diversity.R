#!/usr/bin/env Rscript

# ==============================================================================
# D50 clonal diversity analysis for baseline EBV-specific T cells
# ==============================================================================
#
# Purpose:
#   This script calculates D50 clonal diversity for baseline EBV-specific CD8+
#   T cells using TRB-defined clonotypes.
#
#   D50 is calculated as the minimum number of clonotypes required to account for
#   50% of tetramer-specific cells within each sample and tetramer.
#
#   It saves:
#     1. clone-size table per sample/tetramer/clonotype,
#     2. D50 table per sample/tetramer,
#     3. cohort-level D50 summary,
#     4. optional Wilcoxon statistics comparing MS vs Control,
#     5. session information.
#
# Expected input:
#   A baseline EBV TCR metadata table, for example:
#     Baseline_EBV_TCR_module_score.csv
#
# Notes:
#   - Clonotypes are defined using TRB CDR3 amino acid + TRBV + TRBJ.
#   - Lower D50 means stronger clonal focusing.
#   - Higher D50 means broader clonal diversity.
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
  library(purrr)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  dataset_id = "baseline_ebv",
  
  input_tcr_file = "path/to/input/Baseline_EBV_TCR_module_score.csv",
  
  output_dir = "path/to/output/Publication_data/EBV_baseline/D50_clonal_diversity",
  
  clone_sizes_csv = "path/to/output/Publication_data/EBV_baseline/D50_clonal_diversity/Baseline_EBV_TRB_clone_sizes.csv",
  d50_csv = "path/to/output/Publication_data/EBV_baseline/D50_clonal_diversity/D50_baseline_EBV.csv",
  d50_summary_csv = "path/to/output/Publication_data/EBV_baseline/D50_clonal_diversity/D50_baseline_EBV_summary_by_cohort.csv",
  d50_stats_csv = "path/to/output/Publication_data/EBV_baseline/D50_clonal_diversity/D50_baseline_EBV_MS_vs_Control_stats.csv",
  session_info_file = "path/to/output/Publication_data/EBV_baseline/D50_clonal_diversity/sessionInfo_D50_baseline_EBV.txt",
  
  tetramers_keep = c(
    "GLCT", "CLGG", "YLQQ", "FLRG", "FLYA", "QAKW",
    "YNLR*", "RPPI", "RPQK*", "LLDF", "RAKF"
  ),
  
  # Required metadata columns.
  sample_id_col = "id",
  cohort_col = "cohort",
  tetramer_col = "tetramer",
  antigen_col = "antigen",
  
  # TRB clonotype definition.
  trb_cdr3_col = "TRB_cdr3",
  trb_v_col = "TRB_v_gene",
  trb_j_col = "TRB_j_gene",
  
  # Optional minimum total tetramer-specific cells per sample/tetramer.
  # Set to 10 if you want to enforce a minimum cell count for diversity estimates.
  min_total_tetramer_cells = 1,
  
  # Cohort comparison.
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
  create_parent_dir(cfg$d50_csv)
  create_parent_dir(cfg$d50_summary_csv)
  create_parent_dir(cfg$d50_stats_csv)
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


make_trb_clone_key <- function(df, cfg) {
  df |>
    mutate(
      clone_key_trb = paste(
        .data[[cfg$trb_cdr3_col]],
        .data[[cfg$trb_v_col]],
        .data[[cfg$trb_j_col]],
        sep = "|"
      )
    )
}


filter_tcr_table <- function(df, cfg) {
  required_cols <- c(
    cfg$sample_id_col,
    cfg$cohort_col,
    cfg$tetramer_col,
    cfg$antigen_col,
    cfg$trb_cdr3_col,
    cfg$trb_v_col,
    cfg$trb_j_col
  )
  
  check_required_columns(df, required_cols)
  
  df |>
    filter(
      .data[[cfg$tetramer_col]] %in% cfg$tetramers_keep,
      !is.na(.data[[cfg$sample_id_col]]),
      !is.na(.data[[cfg$cohort_col]]),
      !is.na(.data[[cfg$antigen_col]]),
      !is.na(.data[[cfg$trb_cdr3_col]]),
      .data[[cfg$trb_cdr3_col]] != "",
      !is.na(.data[[cfg$trb_v_col]]),
      .data[[cfg$trb_v_col]] != "",
      !is.na(.data[[cfg$trb_j_col]]),
      .data[[cfg$trb_j_col]] != ""
    ) |>
    make_trb_clone_key(cfg)
}


calculate_clone_sizes <- function(df, cfg) {
  df |>
    count(
      id = .data[[cfg$sample_id_col]],
      cohort = .data[[cfg$cohort_col]],
      tetramer = .data[[cfg$tetramer_col]],
      antigen = .data[[cfg$antigen_col]],
      clone_key_trb,
      name = "clone_size"
    )
}


calc_D50_count <- function(clone_sizes) {
  total_cells <- sum(clone_sizes, na.rm = TRUE)
  
  if (total_cells == 0 || length(clone_sizes) == 0) {
    return(NA_integer_)
  }
  
  cumulative_fraction <- cumsum(sort(clone_sizes, decreasing = TRUE)) / total_cells
  
  which(cumulative_fraction >= 0.5)[1]
}


calculate_D50 <- function(clone_sizes, cfg) {
  d50_df <- clone_sizes |>
    group_by(id, cohort, tetramer, antigen) |>
    summarise(
      total_tetramer_cells = sum(clone_size),
      n_clonotypes = n(),
      D50_count = calc_D50_count(clone_size),
      D50_fraction = D50_count / n_clonotypes,
      D50_percent = 100 * D50_fraction,
      .groups = "drop"
    ) |>
    filter(total_tetramer_cells >= cfg$min_total_tetramer_cells)
  
  d50_df
}


summarise_D50_by_cohort <- function(d50_df) {
  d50_df |>
    group_by(tetramer, antigen, cohort) |>
    summarise(
      n_samples = n_distinct(id),
      median_total_tetramer_cells = median(total_tetramer_cells, na.rm = TRUE),
      median_n_clonotypes = median(n_clonotypes, na.rm = TRUE),
      median_D50_count = median(D50_count, na.rm = TRUE),
      median_D50_percent = median(D50_percent, na.rm = TRUE),
      mean_D50_percent = mean(D50_percent, na.rm = TRUE),
      sd_D50_percent = sd(D50_percent, na.rm = TRUE),
      .groups = "drop"
    )
}


run_wilcoxon_by_tetramer <- function(d50_df, cfg) {
  if (!isTRUE(cfg$run_wilcoxon)) {
    return(tibble())
  }
  
  d50_df <- d50_df |>
    filter(cohort %in% cfg$cohort_levels) |>
    mutate(cohort = factor(cohort, levels = cfg$cohort_levels))
  
  d50_df |>
    group_by(tetramer, antigen) |>
    group_modify(~ {
      dat <- .x
      
      n_control <- sum(dat$cohort == cfg$cohort_levels[1])
      n_ms <- sum(dat$cohort == cfg$cohort_levels[2])
      
      if (n_control < 2 || n_ms < 2) {
        return(tibble(
          n_control = n_control,
          n_ms = n_ms,
          wilcox_p = NA_real_,
          median_control = median(dat$D50_percent[dat$cohort == cfg$cohort_levels[1]], na.rm = TRUE),
          median_ms = median(dat$D50_percent[dat$cohort == cfg$cohort_levels[2]], na.rm = TRUE),
          median_difference_ms_minus_control = NA_real_
        ))
      }
      
      test <- wilcox.test(
        D50_percent ~ cohort,
        data = dat,
        exact = FALSE
      )
      
      median_control <- median(
        dat$D50_percent[dat$cohort == cfg$cohort_levels[1]],
        na.rm = TRUE
      )
      
      median_ms <- median(
        dat$D50_percent[dat$cohort == cfg$cohort_levels[2]],
        na.rm = TRUE
      )
      
      tibble(
        n_control = n_control,
        n_ms = n_ms,
        wilcox_p = test$p.value,
        median_control = median_control,
        median_ms = median_ms,
        median_difference_ms_minus_control = median_ms - median_control
      )
    }) |>
    ungroup() |>
    mutate(
      p_adj_BH = p.adjust(wilcox_p, method = "BH")
    )
}


# ----------------------------- #
# 3. Run pipeline
# ----------------------------- #

setup_dirs(config)

message("Loading baseline EBV TCR table...")
df <- read_tcr_table(config$input_tcr_file)

message("Filtering TCR table and defining TRB clonotypes...")
df_filt <- filter_tcr_table(df, config)

message("Calculating clone sizes...")
clone_sizes <- calculate_clone_sizes(df_filt, config)

write.csv(
  clone_sizes,
  config$clone_sizes_csv,
  row.names = FALSE
)

message("Calculating D50...")
D50_df <- calculate_D50(clone_sizes, config)

write.csv(
  D50_df,
  config$d50_csv,
  row.names = FALSE
)

message("Summarising D50 by cohort...")
D50_summary <- summarise_D50_by_cohort(D50_df)

write.csv(
  D50_summary,
  config$d50_summary_csv,
  row.names = FALSE
)

message("Running optional MS vs Control Wilcoxon tests...")
D50_stats <- run_wilcoxon_by_tetramer(D50_df, config)

write.csv(
  D50_stats,
  config$d50_stats_csv,
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  config$session_info_file
)

message("\nD50 clonal diversity analysis complete.")
message("Saved clone sizes to: ", config$clone_sizes_csv)
message("Saved D50 table to: ", config$d50_csv)
message("Saved cohort summary to: ", config$d50_summary_csv)
message("Saved statistics to: ", config$d50_stats_csv)
message("Saved session info to: ", config$session_info_file)