#!/usr/bin/env Rscript

# ==============================================================================
# Mixed association matrix and cohort sensitivity analyses
# ==============================================================================
#
# Purpose:
#   This script performs mixed association testing between cohort, EBV antibody
#   titres and latent-specific CD8+ T-cell frequency, then generates a lower-
#   triangle square-tile association matrix.
#
# Analyses:
#   1. Overall mixed association matrix in fully matched samples.
#   2. Cohort-only sensitivity analyses for T_cell_freq associations:
#        - Overall
#        - MS only
#        - Control only
#
# Input:
#   correlation_variables.xlsx supplied in source_data sheet Fig1SA-C
#     - sample
#     - cohort
#     - T_cell_freq
#     - EBNA1_IgG
#     - VCA_IgG
#     - EA_IgG
#
# Tests:
#   - Continuous-continuous:
#       Spearman correlation.
#       Effect = Spearman rho.
#
#   - Binary-continuous:
#       Wilcoxon rank-sum test.
#       Effect = rank-biserial correlation.
#
# Cohort sensitivity rules:
#   - IgG vs T_cell_freq = Spearman correlation.
#   - cohort vs T_cell_freq = Wilcoxon rank-sum test.
#   - cohort is not tested inside MS-only or Control-only strata because it is fixed.
#
# Matrix plot:
#   - Lower triangle only.
#   - Square tiles are coloured by association effect size.
#   - Diagonal self-comparisons are fixed at rho/effect = 1.
#   - Numeric values are shown inside tiles.
#   - Asterisks show raw P-values:
#       *   P < 0.05
#       **  P < 0.01
#       *** P < 0.001
#
# Notes:
#   - Paths are intentionally generic for GitHub.
#   - Update the config block before running.
#   - Exact raw and FDR-adjusted P-values are exported.
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
  library(purrr)
  library(ggplot2)
  library(openxlsx)
  library(scales)
  library(stringr)
  library(readr)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  
  # Input Excel file
  input_file = "path/to/correlation_variables.xlsx",
  
  # Sheet name or index.
  # Use NULL to read the first sheet.
  input_sheet = NULL,
  
  # Output directory.
  output_dir = "outputs/mixed_association_matrix",
  
  # Minimum complete samples required.
  min_complete_samples = 6,
  
  # Column candidates. The first matching column is used.
  sample_col_candidates = c("sample", "Sample", "id", "ID", "donor_id"),
  cohort_col_candidates = c("cohort", "diagnosis", "group"),
  tcell_col_candidates = c(
    "T_cell_freq",
    "T cell freq",
    "T_cell_frequency",
    "T cell frequency",
    "Latent CD8 T cell frequency",
    "Latent CD8+ T cell frequency"
  ),
  ebna1_col_candidates = c(
    "EBNA1_IgG",
    "EBNA1 IgG",
    "EBNA1_titre",
    "EBNA1 titre"
  ),
  vca_col_candidates = c(
    "VCA_IgG",
    "VCA IgG",
    "VCA_titre",
    "VCA titre"
  ),
  ea_col_candidates = c(
    "EA_IgG",
    "EA IgG",
    "EA_titre",
    "EA titre"
  ),
  
  # Output filenames
  matrix_png = "overall_matrix_Tcell_cohort_square_effect.png",
  matrix_pdf = "overall_matrix_Tcell_cohort_square_effect.pdf",
  workbook_xlsx = "overall_matrix_plus_cohort_sensitivity_results.xlsx",
  matched_samples_csv = "matched_samples_used_for_matrix.csv",
  association_results_csv = "overall_matrix_results.csv",
  sensitivity_results_csv = "cohort_sensitivity_results.csv",
  sensitivity_counts_csv = "cohort_sensitivity_counts.csv",
  plot_source_csv = "matrix_plot_source_data.csv",
  session_info_file = "sessionInfo_overall_matrix_cohort_sensitivity.txt",
  
  # Plot size
  plot_width = 8.4,
  plot_height = 7.4,
  plot_dpi = 600
)


# ----------------------------- #
# 2. Helper functions
# ----------------------------- #

create_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}


read_input_excel <- function(input_file, input_sheet = NULL) {
  
  if (!file.exists(input_file)) {
    stop("Input file does not exist: ", input_file)
  }
  
  if (is.null(input_sheet)) {
    readxl::read_xlsx(input_file)
  } else {
    readxl::read_xlsx(input_file, sheet = input_sheet)
  }
}


find_col_any <- function(df, possible_names) {
  
  hits <- names(df)[tolower(names(df)) %in% tolower(possible_names)]
  
  if (length(hits) == 0) {
    stop(
      "Missing required column. Tried: ",
      paste(possible_names, collapse = ", ")
    )
  }
  
  hits[1]
}


to_continuous <- function(x) {
  
  z <- trimws(as.character(x))
  z[z %in% c("", "NA", "N/A", "na", "n/a", "NaN", "nan", "null", "NULL")] <- NA
  z <- gsub(",", "", z)
  z <- gsub("<", "", z)
  z <- gsub(">", "", z)
  
  suppressWarnings(as.numeric(z))
}


to_binary_cohort <- function(x) {
  
  if (is.numeric(x)) {
    return(
      ifelse(
        is.na(x),
        NA_real_,
        ifelse(x == 1, 1, ifelse(x == 0, 0, NA_real_))
      )
    )
  }
  
  z <- tolower(trimws(as.character(x)))
  z[z %in% c("", "na", "n/a", "nan", "null", "NA")] <- NA
  
  out <- rep(NA_real_, length(z))
  
  out[z %in% c(
    "1",
    "ms",
    "pwms",
    "rrms",
    "multiple sclerosis",
    "case",
    "patient"
  )] <- 1
  
  out[z %in% c(
    "0",
    "control",
    "controls",
    "ctrl",
    "healthy control",
    "healthy controls",
    "hc",
    "nms",
    "non-ms",
    "nonms"
  )] <- 0
  
  out
}


format_p <- function(p) {
  
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}


p_to_stars <- function(p) {
  
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}


rank_biserial <- function(binary_x, continuous_y) {
  
  keep <- complete.cases(binary_x, continuous_y)
  x <- binary_x[keep]
  y <- continuous_y[keep]
  
  if (length(unique(x)) < 2) return(NA_real_)
  
  n1 <- sum(x == 1)
  n0 <- sum(x == 0)
  
  if (n1 == 0 || n0 == 0) return(NA_real_)
  
  ranks <- rank(y, ties.method = "average")
  R1 <- sum(ranks[x == 1])
  U1 <- R1 - n1 * (n1 + 1) / 2
  
  (2 * U1 / (n1 * n0)) - 1
}


safe_spearman <- function(data, x, y) {
  
  tmp <- data %>%
    filter(
      !is.na(.data[[x]]),
      !is.na(.data[[y]])
    )
  
  n <- nrow(tmp)
  
  if (
    n < 3 ||
    length(unique(tmp[[x]])) < 2 ||
    length(unique(tmp[[y]])) < 2
  ) {
    return(data.frame(
      test = "Spearman correlation",
      n = n,
      effect = NA_real_,
      raw_p = NA_real_,
      note = "Not tested: too few samples or no variation"
    ))
  }
  
  test_res <- suppressWarnings(
    cor.test(
      tmp[[x]],
      tmp[[y]],
      method = "spearman",
      exact = FALSE
    )
  )
  
  data.frame(
    test = "Spearman correlation",
    n = n,
    effect = unname(test_res$estimate),
    raw_p = test_res$p.value,
    note = NA_character_
  )
}


safe_wilcox <- function(data, binary, continuous) {
  
  tmp <- data %>%
    filter(
      !is.na(.data[[binary]]),
      !is.na(.data[[continuous]])
    )
  
  n <- nrow(tmp)
  
  if (
    n < 3 ||
    length(unique(tmp[[binary]])) < 2 ||
    length(unique(tmp[[continuous]])) < 2
  ) {
    return(data.frame(
      test = "Wilcoxon rank-sum test",
      n = n,
      effect = NA_real_,
      raw_p = NA_real_,
      note = "Not tested: too few samples or no variation"
    ))
  }
  
  test_res <- suppressWarnings(
    wilcox.test(
      tmp[[continuous]] ~ factor(tmp[[binary]]),
      exact = FALSE
    )
  )
  
  data.frame(
    test = "Wilcoxon rank-sum test",
    n = n,
    effect = rank_biserial(tmp[[binary]], tmp[[continuous]]),
    raw_p = test_res$p.value,
    note = NA_character_
  )
}


# ----------------------------- #
# 3. Paths
# ----------------------------- #

create_dir(config$output_dir)

out_plot_png <- file.path(config$output_dir, config$matrix_png)
out_plot_pdf <- file.path(config$output_dir, config$matrix_pdf)
out_excel <- file.path(config$output_dir, config$workbook_xlsx)
out_matched_csv <- file.path(config$output_dir, config$matched_samples_csv)
out_association_csv <- file.path(config$output_dir, config$association_results_csv)
out_sensitivity_csv <- file.path(config$output_dir, config$sensitivity_results_csv)
out_sensitivity_counts_csv <- file.path(config$output_dir, config$sensitivity_counts_csv)
out_plot_source_csv <- file.path(config$output_dir, config$plot_source_csv)
out_session_info <- file.path(config$output_dir, config$session_info_file)


# ----------------------------- #
# 4. Read and standardise data
# ----------------------------- #

df_raw <- read_input_excel(
  input_file = config$input_file,
  input_sheet = config$input_sheet
) %>%
  as.data.frame()

names(df_raw) <- trimws(names(df_raw))

sample_col <- find_col_any(df_raw, config$sample_col_candidates)
cohort_col <- find_col_any(df_raw, config$cohort_col_candidates)
tcell_col <- find_col_any(df_raw, config$tcell_col_candidates)
ebna1_col <- find_col_any(df_raw, config$ebna1_col_candidates)
vca_col <- find_col_any(df_raw, config$vca_col_candidates)
ea_col <- find_col_any(df_raw, config$ea_col_candidates)

df_raw <- df_raw %>%
  transmute(
    sample = .data[[sample_col]],
    cohort = .data[[cohort_col]],
    T_cell_freq = .data[[tcell_col]],
    EBNA1_IgG = .data[[ebna1_col]],
    VCA_IgG = .data[[vca_col]],
    EA_IgG = .data[[ea_col]]
  )


# ----------------------------- #
# 5. Variable definitions
# ----------------------------- #

vars <- c(
  "cohort",
  "EBNA1_IgG",
  "VCA_IgG",
  "EA_IgG",
  "T_cell_freq"
)

binary_vars <- c("cohort")

continuous_vars <- c(
  "EBNA1_IgG",
  "VCA_IgG",
  "EA_IgG",
  "T_cell_freq"
)

var_labels <- c(
  "cohort" = "Cohort",
  "EBNA1_IgG" = "EBNA1\nIgG",
  "VCA_IgG" = "VCA\nIgG",
  "EA_IgG" = "EA\nIgG",
  "T_cell_freq" = "Latent CD8+ T cell\nfrequency"
)

var_labels <- var_labels[vars]

get_type <- function(v) {
  
  if (v %in% binary_vars) return("binary")
  if (v %in% continuous_vars) return("continuous")
  
  NA_character_
}


# ----------------------------- #
# 6. Clean data and keep fully matched samples
# ----------------------------- #

df <- df_raw %>%
  mutate(
    sample = as.character(sample),
    cohort = to_binary_cohort(cohort),
    EBNA1_IgG = to_continuous(EBNA1_IgG),
    VCA_IgG = to_continuous(VCA_IgG),
    EA_IgG = to_continuous(EA_IgG),
    T_cell_freq = to_continuous(T_cell_freq),
    cohort_label = case_when(
      cohort == 1 ~ "MS",
      cohort == 0 ~ "Control",
      TRUE ~ NA_character_
    ),
    cohort_label = factor(
      cohort_label,
      levels = c("Control", "MS")
    )
  )

df_matched <- df %>%
  filter(if_all(all_of(vars), ~ !is.na(.x))) %>%
  filter(cohort %in% c(0, 1))

message("Original samples: ", nrow(df))
message("Fully matched samples used: ", nrow(df_matched))

if (nrow(df_matched) < config$min_complete_samples) {
  stop("Too few fully matched samples remain.")
}


# ----------------------------- #
# 7. Counts and descriptive summaries
# ----------------------------- #

counts_cohort <- df_matched %>%
  count(cohort_label, name = "n")

cohort_T_cell_summary <- df_matched %>%
  group_by(cohort_label) %>%
  summarise(
    n = n(),
    median_T_cell_freq = median(T_cell_freq, na.rm = TRUE),
    IQR_T_cell_freq = IQR(T_cell_freq, na.rm = TRUE),
    mean_T_cell_freq = mean(T_cell_freq, na.rm = TRUE),
    sd_T_cell_freq = sd(T_cell_freq, na.rm = TRUE),
    .groups = "drop"
  )

print(counts_cohort)
print(cohort_T_cell_summary)


# ==============================================================================
# PART A: Overall mixed association matrix
# ==============================================================================

pair_association <- function(var1, var2, data) {
  
  x <- data[[var1]]
  y <- data[[var2]]
  
  type1 <- get_type(var1)
  type2 <- get_type(var2)
  
  keep <- complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  
  n <- length(x)
  
  if (n < 3) {
    return(data.frame(
      var1 = var1,
      var2 = var2,
      type1 = type1,
      type2 = type2,
      test = "Not tested",
      effect = NA_real_,
      raw_p = NA_real_,
      n = n,
      direction = "Insufficient complete observations"
    ))
  }
  
  # Continuous-continuous associations.
  if (type1 == "continuous" && type2 == "continuous") {
    
    if (length(unique(x)) < 2 || length(unique(y)) < 2) {
      return(data.frame(
        var1 = var1,
        var2 = var2,
        type1 = type1,
        type2 = type2,
        test = "Not tested",
        effect = NA_real_,
        raw_p = NA_real_,
        n = n,
        direction = "No variation"
      ))
    }
    
    test_res <- suppressWarnings(
      cor.test(x, y, method = "spearman", exact = FALSE)
    )
    
    return(data.frame(
      var1 = var1,
      var2 = var2,
      type1 = type1,
      type2 = type2,
      test = "Spearman correlation",
      effect = unname(test_res$estimate),
      raw_p = test_res$p.value,
      n = n,
      direction = paste0(
        "Positive = ",
        var1,
        " increases as ",
        var2,
        " increases"
      )
    ))
  }
  
  # Binary-continuous associations.
  if (type1 == "binary" && type2 == "continuous") {
    
    if (length(unique(x)) < 2 || length(unique(y)) < 2) {
      return(data.frame(
        var1 = var1,
        var2 = var2,
        type1 = type1,
        type2 = type2,
        test = "Not tested",
        effect = NA_real_,
        raw_p = NA_real_,
        n = n,
        direction = "No variation"
      ))
    }
    
    test_res <- suppressWarnings(
      wilcox.test(y ~ factor(x), exact = FALSE)
    )
    
    effect <- rank_biserial(x, y)
    
    return(data.frame(
      var1 = var1,
      var2 = var2,
      type1 = type1,
      type2 = type2,
      test = "Wilcoxon rank-sum test",
      effect = effect,
      raw_p = test_res$p.value,
      n = n,
      direction = paste0(
        "Positive = ",
        var2,
        " higher when ",
        var1,
        " = 1"
      )
    ))
  }
  
  if (type1 == "continuous" && type2 == "binary") {
    
    if (length(unique(x)) < 2 || length(unique(y)) < 2) {
      return(data.frame(
        var1 = var1,
        var2 = var2,
        type1 = type1,
        type2 = type2,
        test = "Not tested",
        effect = NA_real_,
        raw_p = NA_real_,
        n = n,
        direction = "No variation"
      ))
    }
    
    test_res <- suppressWarnings(
      wilcox.test(x ~ factor(y), exact = FALSE)
    )
    
    effect <- rank_biserial(y, x)
    
    return(data.frame(
      var1 = var1,
      var2 = var2,
      type1 = type1,
      type2 = type2,
      test = "Wilcoxon rank-sum test",
      effect = effect,
      raw_p = test_res$p.value,
      n = n,
      direction = paste0(
        "Positive = ",
        var1,
        " higher when ",
        var2,
        " = 1"
      )
    ))
  }
  
  data.frame(
    var1 = var1,
    var2 = var2,
    type1 = type1,
    type2 = type2,
    test = "Not tested",
    effect = NA_real_,
    raw_p = NA_real_,
    n = n,
    direction = "Unsupported variable type combination"
  )
}

pair_grid <- expand.grid(
  i = seq_along(vars),
  j = seq_along(vars)
) %>%
  filter(i > j) %>%
  mutate(
    var1 = vars[i],
    var2 = vars[j]
  )

association_results <- purrr::map2_dfr(
  pair_grid$var1,
  pair_grid$var2,
  ~ pair_association(.x, .y, df_matched)
) %>%
  mutate(
    BH_p = p.adjust(raw_p, method = "BH"),
    raw_p_label = format_p(raw_p),
    BH_p_label = format_p(BH_p),
    sig_raw = p_to_stars(raw_p),
    sig_BH = p_to_stars(BH_p),
    abs_effect = abs(effect),
    effect_label = ifelse(is.na(effect), "", sprintf("%.2f", effect))
  )

cohort_T_cell_test <- association_results %>%
  filter(
    (var1 == "T_cell_freq" & var2 == "cohort") |
      (var1 == "cohort" & var2 == "T_cell_freq")
  )

print(cohort_T_cell_test)


# ==============================================================================
# PART B: Cohort sensitivity analyses
# ==============================================================================

run_sensitivity_set <- function(data, stratum_name, test_cohort = FALSE) {
  
  out <- bind_rows(
    safe_spearman(data, "EBNA1_IgG", "T_cell_freq") %>%
      mutate(
        stratum = stratum_name,
        exposure = "EBNA1_IgG",
        outcome = "T_cell_freq",
        interpretation = "Spearman rho: positive means higher T cell frequency with higher EBNA1 IgG"
      ),
    
    safe_spearman(data, "VCA_IgG", "T_cell_freq") %>%
      mutate(
        stratum = stratum_name,
        exposure = "VCA_IgG",
        outcome = "T_cell_freq",
        interpretation = "Spearman rho: positive means higher T cell frequency with higher VCA IgG"
      ),
    
    safe_spearman(data, "EA_IgG", "T_cell_freq") %>%
      mutate(
        stratum = stratum_name,
        exposure = "EA_IgG",
        outcome = "T_cell_freq",
        interpretation = "Spearman rho: positive means higher T cell frequency with higher EA IgG"
      )
  )
  
  if (test_cohort) {
    out <- bind_rows(
      out,
      safe_wilcox(data, "cohort", "T_cell_freq") %>%
        mutate(
          stratum = stratum_name,
          exposure = "cohort",
          outcome = "T_cell_freq",
          interpretation = "Rank-biserial effect: positive means higher T cell frequency in MS"
        )
    )
  }
  
  out
}

sensitivity_results <- bind_rows(
  run_sensitivity_set(
    data = df_matched,
    stratum_name = "Overall",
    test_cohort = TRUE
  ),
  
  run_sensitivity_set(
    data = df_matched %>% filter(cohort == 1),
    stratum_name = "MS only",
    test_cohort = FALSE
  ),
  
  run_sensitivity_set(
    data = df_matched %>% filter(cohort == 0),
    stratum_name = "Control only",
    test_cohort = FALSE
  )
) %>%
  group_by(stratum) %>%
  mutate(
    FDR_within_stratum = p.adjust(raw_p, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    FDR_all_sensitivity = p.adjust(raw_p, method = "BH"),
    raw_p_label = format_p(raw_p),
    FDR_within_stratum_label = format_p(FDR_within_stratum),
    FDR_all_sensitivity_label = format_p(FDR_all_sensitivity),
    sig_raw = p_to_stars(raw_p),
    effect_label = ifelse(is.na(effect), "NA", sprintf("%.2f", effect))
  ) %>%
  select(
    stratum,
    exposure,
    outcome,
    test,
    n,
    effect,
    effect_label,
    raw_p,
    raw_p_label,
    sig_raw,
    FDR_within_stratum,
    FDR_within_stratum_label,
    FDR_all_sensitivity,
    FDR_all_sensitivity_label,
    interpretation,
    note
  )

print(sensitivity_results)


summarise_stratum <- function(data, stratum_name) {
  
  data %>%
    summarise(
      stratum = stratum_name,
      n = n(),
      n_MS = sum(cohort == 1, na.rm = TRUE),
      n_Control = sum(cohort == 0, na.rm = TRUE),
      median_T_cell_freq = ifelse(
        n() > 0,
        median(T_cell_freq, na.rm = TRUE),
        NA_real_
      ),
      IQR_T_cell_freq = ifelse(
        n() > 0,
        IQR(T_cell_freq, na.rm = TRUE),
        NA_real_
      ),
      mean_T_cell_freq = ifelse(
        n() > 0,
        mean(T_cell_freq, na.rm = TRUE),
        NA_real_
      ),
      sd_T_cell_freq = ifelse(
        n() > 1,
        sd(T_cell_freq, na.rm = TRUE),
        NA_real_
      )
    )
}

sensitivity_counts <- bind_rows(
  summarise_stratum(df_matched, "Overall"),
  summarise_stratum(df_matched %>% filter(cohort == 1), "MS only"),
  summarise_stratum(df_matched %>% filter(cohort == 0), "Control only")
)

print(sensitivity_counts)


# ==============================================================================
# PART C: Lower-triangle square association matrix plot
# ==============================================================================

diag_results <- data.frame(
  var1 = vars,
  var2 = vars,
  type1 = sapply(vars, get_type),
  type2 = sapply(vars, get_type),
  test = "Self-correlation",
  effect = 1,
  raw_p = NA_real_,
  n = nrow(df_matched),
  direction = "Self-correlation"
) %>%
  mutate(
    BH_p = NA_real_,
    raw_p_label = "",
    BH_p_label = "",
    sig_raw = "",
    sig_BH = "",
    abs_effect = 1,
    effect_label = "1.00"
  )

association_plot_results <- association_results %>%
  mutate(
    rho = effect,
    abs_rho = abs(effect),
    rho_label = ifelse(is.na(rho), "", sprintf("%.2f", rho))
  )

plot_lower <- bind_rows(
  association_plot_results %>%
    transmute(
      x = var2,
      y = var1,
      rho,
      abs_rho,
      raw_p,
      sig_raw,
      test,
      rho_label
    ),
  diag_results %>%
    transmute(
      x = var2,
      y = var1,
      rho = effect,
      abs_rho = abs(effect),
      raw_p,
      sig_raw,
      test,
      rho_label = effect_label
    )
) %>%
  mutate(
    x = factor(x, levels = vars),
    y = factor(y, levels = rev(vars))
  )

custom_cols <- colorRampPalette(c(
  "#011f4b",
  "#03396c",
  "#005b96",
  "#6497b1",
  "#b3cde0",
  "#ffefea",
  "#fbd9d3",
  "#ffb09c",
  "#fe5757",
  "#cb2424",
  "#900000"
))(100)

p <- ggplot() +
  
  geom_tile(
    data = plot_lower,
    aes(
      x = x,
      y = y,
      fill = rho
    ),
    colour = "grey45",
    linewidth = 0.55
  ) +
  
  geom_text(
    data = plot_lower,
    aes(
      x = x,
      y = y,
      label = rho_label
    ),
    colour = "black",
    size = 4.3,
    fontface = "bold"
  ) +
  
  geom_text(
    data = plot_lower %>%
      filter(sig_raw != "", test != "Self-correlation"),
    aes(
      x = x,
      y = y,
      label = sig_raw
    ),
    colour = "black",
    size = 4.6,
    fontface = "bold",
    nudge_y = 0.28
  ) +
  
  scale_fill_gradientn(
    colours = custom_cols,
    limits = c(-1, 1),
    oob = scales::squish,
    breaks = c(-1, -0.5, 0, 0.5, 1),
    name = expression(rho)
  ) +
  
  scale_x_discrete(
    position = "top",
    drop = FALSE,
    labels = var_labels
  ) +
  
  scale_y_discrete(
    drop = FALSE,
    labels = var_labels
  ) +
  
  coord_fixed() +
  
  labs(
    x = NULL,
    y = NULL,
    title = "Matched-only mixed association matrix",
    subtitle = paste0(
      "Fully matched samples only: n = ",
      nrow(df_matched),
      ". Lower triangle shows association effect as coloured squares. ",
      "Numbers are Spearman rho for continuous-continuous tests and rank-biserial effect for binary-continuous tests. ",
      "Asterisks indicate raw P-values: * < 0.05, ** < 0.01, *** < 0.001."
    )
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.text.x = element_text(
      angle = 45,
      hjust = 0,
      vjust = 0,
      size = 15,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 15,
      colour = "black"
    ),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 10.5),
    plot.margin = margin(10, 14, 10, 10)
  )

print(p)

ggsave(
  filename = out_plot_png,
  plot = p,
  width = config$plot_width,
  height = config$plot_height,
  dpi = config$plot_dpi,
  bg = "white"
)

ggsave(
  filename = out_plot_pdf,
  plot = p,
  width = config$plot_width,
  height = config$plot_height,
  bg = "white"
)


# ----------------------------- #
# 8. Export CSV files
# ----------------------------- #

write_csv(df_matched, out_matched_csv)
write_csv(association_results, out_association_csv)
write_csv(sensitivity_results, out_sensitivity_csv)
write_csv(sensitivity_counts, out_sensitivity_counts_csv)
write_csv(plot_lower, out_plot_source_csv)


# ==============================================================================
# PART D: Save Excel workbook
# ==============================================================================

wb <- createWorkbook()

addWorksheet(wb, "Matched_samples")
writeData(wb, "Matched_samples", df_matched)

addWorksheet(wb, "Overall_matrix_results")
writeData(wb, "Overall_matrix_results", association_results)

addWorksheet(wb, "Matrix_plot_source")
writeData(wb, "Matrix_plot_source", plot_lower)

addWorksheet(wb, "Cohort_T_cell_freq_test")
writeData(wb, "Cohort_T_cell_freq_test", cohort_T_cell_test)

addWorksheet(wb, "Cohort_T_cell_freq_summary")
writeData(wb, "Cohort_T_cell_freq_summary", cohort_T_cell_summary)

addWorksheet(wb, "Counts_cohort")
writeData(wb, "Counts_cohort", counts_cohort)

addWorksheet(wb, "Sensitivity_results")
writeData(wb, "Sensitivity_results", sensitivity_results)

addWorksheet(wb, "Sensitivity_counts")
writeData(wb, "Sensitivity_counts", sensitivity_counts)

addWorksheet(wb, "Notes")
writeData(
  wb,
  "Notes",
  data.frame(
    item = c(
      "Input file",
      "Input sheet",
      "Sample column used",
      "Cohort column used",
      "T-cell frequency column used",
      "Overall matrix",
      "Cohort-only sensitivity analysis",
      "Overall sensitivity",
      "Cohort-stratified sensitivity",
      "Continuous-continuous tests",
      "Binary-continuous tests",
      "Cohort coding",
      "FDR correction",
      "Asterisks",
      "Matrix visualisation"
    ),
    note = c(
      config$input_file,
      ifelse(is.null(config$input_sheet), "First sheet", as.character(config$input_sheet)),
      sample_col,
      cohort_col,
      tcell_col,
      "The matrix shows overall association patterns in fully matched samples.",
      "Sensitivity analyses repeat associations with latent CD8+ T cell frequency across overall, MS-only, and Control-only strata.",
      "Overall sensitivity includes EBNA1/VCA/EA IgG correlations with T_cell_freq and cohort comparison of T_cell_freq.",
      "In MS-only and Control-only strata, only EBNA1/VCA/EA IgG correlations with T_cell_freq are performed because cohort is fixed.",
      "Continuous-continuous associations use Spearman correlation. Effect = Spearman rho.",
      "Binary-continuous associations use Wilcoxon rank-sum test. Effect = rank-biserial correlation.",
      "Cohort is coded 1 = MS, 0 = Control.",
      "FDR_within_stratum adjusts P-values within each sensitivity stratum. FDR_all_sensitivity adjusts across all sensitivity tests.",
      "Asterisks on the matrix indicate raw P-values: * < 0.05, ** < 0.01, *** < 0.001. Exact raw and FDR P-values are saved in Excel.",
      "The lower triangle uses square tiles coloured by effect size. Numeric values are shown inside each tile."
    )
  )
)

saveWorkbook(wb, out_excel, overwrite = TRUE)


# ----------------------------- #
# 9. Session information
# ----------------------------- #

writeLines(
  c(
    "Configuration:",
    capture.output(str(config)),
    "",
    "Detected columns:",
    paste("sample:", sample_col),
    paste("cohort:", cohort_col),
    paste("T_cell_freq:", tcell_col),
    paste("EBNA1_IgG:", ebna1_col),
    paste("VCA_IgG:", vca_col),
    paste("EA_IgG:", ea_col),
    "",
    "Original samples:",
    as.character(nrow(df)),
    "",
    "Fully matched samples:",
    as.character(nrow(df_matched)),
    "",
    "Cohort counts:",
    capture.output(print(counts_cohort)),
    "",
    "T-cell frequency summary by cohort:",
    capture.output(print(cohort_T_cell_summary)),
    "",
    "Output files:",
    out_plot_png,
    out_plot_pdf,
    out_excel,
    out_matched_csv,
    out_association_csv,
    out_sensitivity_csv,
    out_sensitivity_counts_csv,
    out_plot_source_csv,
    "",
    "Session information:",
    capture.output(sessionInfo())
  ),
  out_session_info
)


# ----------------------------- #
# 10. Completion messages
# ----------------------------- #

message("Saved matrix PNG to: ", out_plot_png)
message("Saved matrix PDF to: ", out_plot_pdf)
message("Saved Excel results to: ", out_excel)
message("Saved session info to: ", out_session_info)
message("Mixed association matrix and cohort sensitivity workflow complete.")
