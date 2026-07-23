#!/usr/bin/env Rscript

# ==============================================================================
# Figure S1A: Spearman correlation matrix with EBV load
# ==============================================================================
#
# Run from the repository root:
#
#   Rscript scripts/figure_s1a_spearman_matrix.R
#
# Or provide custom paths:
#
#   Rscript scripts/figure_s1a_spearman_matrix.R \
#     data/correlation_variables.xlsx \
#     results/figure_s1a \
#     Full
#
# Positional arguments:
#   1. Input Excel file
#      Default: data/correlation_variables.xlsx
#   2. Output directory
#      Default: results/figure_s1a
#   3. Excel sheet name or number
#      Default: 1
#
# Required columns:
#   sample, cohort, T_cell_freq, EBNA1_IgG, VCA_IgG, EA_IgG, EBV_load
#
# Analysis:
#   - overall pairwise Spearman correlations
#   - cohort-stratified sensitivity analyses in MS and Control
#   - raw and Benjamini-Hochberg-adjusted P-values are exported
#
# Plot:
#   - lower triangle plus diagonal only
#   - tile colour represents raw P-value
#   - black text represents Spearman rho
#   - asterisks represent raw P-value significance
#   - FDR values are not used for plotting
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Dependencies
# ------------------------------------------------------------------------------

required_packages <- c(
  "readxl",
  "dplyr",
  "purrr",
  "ggplot2",
  "openxlsx",
  "scales"
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
  library(purrr)
  library(ggplot2)
  library(openxlsx)
  library(scales)
})


# ------------------------------------------------------------------------------
# 2. Command-line arguments and output paths
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

input_file <- if (length(args) >= 1) {
  args[[1]]
} else {
  file.path("data", "correlation_variables.xlsx")
}

output_dir <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path("results", "figure_s1a")
}

input_sheet <- if (length(args) >= 3) {
  if (grepl("^[0-9]+$", args[[3]])) {
    as.integer(args[[3]])
  } else {
    args[[3]]
  }
} else {
  1
}

if (!file.exists(input_file)) {
  stop(
    "Input file does not exist: ",
    input_file,
    "\nRun the script from the repository root or provide a path explicitly."
  )
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_files <- list(
  plot_png = file.path(
    output_dir,
    "overall_spearman_matrix_with_EBV_load.png"
  ),
  plot_pdf = file.path(
    output_dir,
    "overall_spearman_matrix_with_EBV_load.pdf"
  ),
  workbook = file.path(
    output_dir,
    "overall_spearman_matrix_with_EBV_load_results.xlsx"
  ),
  overall_csv = file.path(
    output_dir,
    "overall_matrix_results.csv"
  ),
  sensitivity_csv = file.path(
    output_dir,
    "cohort_sensitivity_results.csv"
  ),
  matched_csv = file.path(
    output_dir,
    "matched_samples.csv"
  ),
  session_info = file.path(
    output_dir,
    "session_info.txt"
  )
)


# ------------------------------------------------------------------------------
# 3. Read data
# ------------------------------------------------------------------------------

message("Reading: ", input_file)
message("Excel sheet: ", input_sheet)

df_raw <- readxl::read_xlsx(
  path = input_file,
  sheet = input_sheet
) %>%
  as.data.frame()

names(df_raw) <- trimws(names(df_raw))


# ------------------------------------------------------------------------------
# 4. Helper functions
# ------------------------------------------------------------------------------

find_column <- function(data, possible_names) {

  hits <- names(data)[
    tolower(names(data)) %in% tolower(possible_names)
  ]

  if (length(hits) == 0) {
    stop(
      "Missing required column. Tried: ",
      paste(possible_names, collapse = ", ")
    )
  }

  hits[[1]]
}

to_continuous <- function(x) {

  value <- trimws(as.character(x))

  value[
    value %in% c(
      "", "NA", "N/A", "na", "n/a",
      "NaN", "nan", "null"
    )
  ] <- NA_character_

  value <- gsub(",", "", value)
  value <- gsub("<", "", value, fixed = TRUE)
  value <- gsub(">", "", value, fixed = TRUE)

  suppressWarnings(as.numeric(value))
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

  value <- tolower(trimws(as.character(x)))

  value[
    value %in% c("", "na", "n/a", "nan", "null")
  ] <- NA_character_

  output <- rep(NA_real_, length(value))

  output[
    value %in% c(
      "1", "ms", "pwms", "rrms",
      "multiple sclerosis", "case", "patient"
    )
  ] <- 1

  output[
    value %in% c(
      "0", "control", "controls", "ctrl",
      "healthy control", "healthy controls",
      "hc", "nms"
    )
  ] <- 0

  output
}

format_p <- function(p) {

  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

p_to_stars <- function(p) {

  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

safe_spearman <- function(data, x, y) {

  test_data <- data %>%
    filter(
      !is.na(.data[[x]]),
      !is.na(.data[[y]])
    )

  n_complete <- nrow(test_data)

  if (
    n_complete < 3 ||
    dplyr::n_distinct(test_data[[x]]) < 2 ||
    dplyr::n_distinct(test_data[[y]]) < 2
  ) {
    return(
      tibble::tibble(
        test = "Spearman correlation",
        n = n_complete,
        rho = NA_real_,
        raw_p = NA_real_,
        note = "Not tested: too few observations or no variation"
      )
    )
  }

  test_result <- suppressWarnings(
    stats::cor.test(
      test_data[[x]],
      test_data[[y]],
      method = "spearman",
      exact = FALSE
    )
  )

  tibble::tibble(
    test = "Spearman correlation",
    n = n_complete,
    rho = unname(test_result$estimate),
    raw_p = test_result$p.value,
    note = NA_character_
  )
}


# ------------------------------------------------------------------------------
# 5. Standardise input columns
# ------------------------------------------------------------------------------

df <- df_raw %>%
  transmute(
    sample = as.character(
      .data[[
        find_column(
          df_raw,
          c("sample", "Sample", "sample_id", "Sample ID")
        )
      ]]
    ),

    cohort = .data[[
      find_column(
        df_raw,
        c("cohort", "diagnosis", "group")
      )
    ]],

    T_cell_freq = .data[[
      find_column(
        df_raw,
        c(
          "T_cell_freq",
          "T cell freq",
          "T_cell_frequency",
          "T cell frequency",
          "Latent CD8 T cell frequency",
          "Latent CD8+ T cell frequency"
        )
      )
    ]],

    EBNA1_IgG = .data[[
      find_column(
        df_raw,
        c(
          "EBNA1_IgG",
          "EBNA1 IgG",
          "EBNA1_titre",
          "EBNA1 titre"
        )
      )
    ]],

    VCA_IgG = .data[[
      find_column(
        df_raw,
        c(
          "VCA_IgG",
          "VCA IgG",
          "VCA_titre",
          "VCA titre"
        )
      )
    ]],

    EA_IgG = .data[[
      find_column(
        df_raw,
        c(
          "EA_IgG",
          "EA IgG",
          "EA_titre",
          "EA titre"
        )
      )
    ]],

    EBV_load = .data[[
      find_column(
        df_raw,
        c(
          "EBV_load",
          "EBV load",
          "EBV DNA",
          "EBV_DNA",
          "Salivary EBV load",
          "Salivary EBV DNA"
        )
      )
    ]]
  ) %>%
  mutate(
    cohort = to_binary_cohort(cohort),
    T_cell_freq = to_continuous(T_cell_freq),
    EBNA1_IgG = to_continuous(EBNA1_IgG),
    VCA_IgG = to_continuous(VCA_IgG),
    EA_IgG = to_continuous(EA_IgG),
    EBV_load = to_continuous(EBV_load),

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


# ------------------------------------------------------------------------------
# 6. Analysis variables and matched dataset
# ------------------------------------------------------------------------------

matrix_vars <- c(
  "EBNA1_IgG",
  "VCA_IgG",
  "EA_IgG",
  "EBV_load",
  "T_cell_freq"
)

required_analysis_vars <- c(
  "cohort",
  matrix_vars
)

variable_labels <- c(
  "EBNA1_IgG" = "EBNA1\nIgG",
  "VCA_IgG" = "VCA\nIgG",
  "EA_IgG" = "EA\nIgG",
  "EBV_load" = "EBV\nload",
  "T_cell_freq" = "Latent CD8+ T cell\nfrequency"
)

variable_labels <- variable_labels[matrix_vars]

df_matched <- df %>%
  filter(
    if_all(
      all_of(required_analysis_vars),
      ~ !is.na(.x)
    ),
    cohort %in% c(0, 1)
  )

message("Original samples: ", nrow(df))
message("Fully matched samples: ", nrow(df_matched))

if (nrow(df_matched) < 6) {
  stop("Too few fully matched samples remain after filtering.")
}

cohort_counts <- df_matched %>%
  count(
    cohort_label,
    name = "n"
  )


# ------------------------------------------------------------------------------
# 7. Overall Spearman correlation tests
# ------------------------------------------------------------------------------

pair_grid <- expand.grid(
  i = seq_along(matrix_vars),
  j = seq_along(matrix_vars),
  stringsAsFactors = FALSE
) %>%
  filter(i > j) %>%
  mutate(
    var1 = matrix_vars[i],
    var2 = matrix_vars[j]
  )

association_results <- purrr::map2_dfr(
  pair_grid$var1,
  pair_grid$var2,
  function(var1, var2) {

    safe_spearman(
      data = df_matched,
      x = var1,
      y = var2
    ) %>%
      mutate(
        var1 = var1,
        var2 = var2,
        .before = 1
      )
  }
) %>%
  mutate(
    FDR_BH = p.adjust(raw_p, method = "BH"),
    raw_p_label = format_p(raw_p),
    FDR_BH_label = format_p(FDR_BH),
    raw_p_stars = p_to_stars(raw_p),
    rho_label = ifelse(
      is.na(rho),
      "",
      sprintf("%.2f", rho)
    )
  ) %>%
  select(
    var1,
    var2,
    test,
    n,
    rho,
    rho_label,
    raw_p,
    raw_p_label,
    raw_p_stars,
    FDR_BH,
    FDR_BH_label,
    note
  )

message("\nOverall Spearman results:")
print(association_results)


# ------------------------------------------------------------------------------
# 8. Cohort-stratified sensitivity analyses
# ------------------------------------------------------------------------------

run_matrix_set <- function(data, stratum_name) {

  purrr::map2_dfr(
    pair_grid$var1,
    pair_grid$var2,
    function(var1, var2) {

      safe_spearman(
        data = data,
        x = var1,
        y = var2
      ) %>%
        mutate(
          stratum = stratum_name,
          var1 = var1,
          var2 = var2,
          .before = 1
        )
    }
  )
}

sensitivity_results <- bind_rows(
  run_matrix_set(
    data = df_matched %>% filter(cohort == 1),
    stratum_name = "MS only"
  ),

  run_matrix_set(
    data = df_matched %>% filter(cohort == 0),
    stratum_name = "Control only"
  )
) %>%
  group_by(stratum) %>%
  mutate(
    FDR_BH_within_stratum = p.adjust(
      raw_p,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  mutate(
    FDR_BH_all_sensitivity = p.adjust(
      raw_p,
      method = "BH"
    ),
    raw_p_label = format_p(raw_p),
    FDR_BH_within_stratum_label = format_p(
      FDR_BH_within_stratum
    ),
    FDR_BH_all_sensitivity_label = format_p(
      FDR_BH_all_sensitivity
    ),
    raw_p_stars = p_to_stars(raw_p),
    rho_label = ifelse(
      is.na(rho),
      "",
      sprintf("%.2f", rho)
    )
  ) %>%
  select(
    stratum,
    var1,
    var2,
    test,
    n,
    rho,
    rho_label,
    raw_p,
    raw_p_label,
    raw_p_stars,
    FDR_BH_within_stratum,
    FDR_BH_within_stratum_label,
    FDR_BH_all_sensitivity,
    FDR_BH_all_sensitivity_label,
    note
  )

sensitivity_counts <- df_matched %>%
  count(
    cohort_label,
    name = "n"
  ) %>%
  transmute(
    stratum = paste(cohort_label, "only"),
    n
  )


# ------------------------------------------------------------------------------
# 9. Correlation matrix plot
# ------------------------------------------------------------------------------

diagonal_results <- tibble::tibble(
  var1 = matrix_vars,
  var2 = matrix_vars,
  test = "Self-correlation",
  n = nrow(df_matched),
  rho = 1,
  rho_label = "1.00",
  raw_p = NA_real_,
  raw_p_label = "",
  raw_p_stars = "",
  FDR_BH = NA_real_,
  FDR_BH_label = "",
  note = NA_character_
)

plot_data <- bind_rows(
  association_results,
  diagonal_results
) %>%
  transmute(
    x = var2,
    y = var1,
    test,
    rho,
    rho_label,
    raw_p,
    raw_p_stars,

    # Use -log10(raw P) for the fill scale.
    # Cap values at P = 0.001 to prevent one very small P-value
    # from compressing the rest of the colour range.
    p_colour = case_when(
      test == "Self-correlation" ~ NA_real_,
      is.na(raw_p) ~ NA_real_,
      TRUE ~ pmin(-log10(raw_p), 3)
    )
  ) %>%
  mutate(
    x = factor(x, levels = matrix_vars),
    y = factor(y, levels = rev(matrix_vars))
  )

plot_diagonal <- plot_data %>%
  filter(test == "Self-correlation")

plot_off_diagonal <- plot_data %>%
  filter(test != "Self-correlation")

p <- ggplot() +

  geom_tile(
    data = plot_off_diagonal,
    aes(
      x = x,
      y = y,
      fill = p_colour
    ),
    colour = "grey45",
    linewidth = 0.55
  ) +

  geom_tile(
    data = plot_diagonal,
    aes(
      x = x,
      y = y
    ),
    fill = "grey85",
    colour = "grey45",
    linewidth = 0.55
  ) +

  geom_text(
    data = plot_data,
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
    data = plot_off_diagonal %>%
      filter(raw_p_stars != ""),
    aes(
      x = x,
      y = y,
      label = raw_p_stars
    ),
    colour = "black",
    size = 4.2,
    fontface = "bold",
    nudge_y = 0.27
  ) +

  scale_fill_gradient(
    low = "#fff5f0",
    high = "#99000d",
    limits = c(0, 3),
    oob = scales::squish,
    breaks = -log10(c(1, 0.05, 0.01, 0.001)),
    labels = c("1.00", "0.05", "0.01", "≤0.001"),
    name = "Raw P-value"
  ) +

  scale_x_discrete(
    position = "top",
    drop = FALSE,
    labels = variable_labels
  ) +

  scale_y_discrete(
    drop = FALSE,
    labels = variable_labels
  ) +

  coord_fixed() +

  labs(
    x = NULL,
    y = NULL,
    title = "Overall Spearman correlation matrix",
    subtitle = paste0(
      "Fully matched samples: n = ",
      nrow(df_matched),
      ". Tile colour represents raw P-value; ",
      "black numbers represent Spearman rho. ",
      "Asterisks indicate raw P-value significance."
    )
  ) +

  theme_minimal(base_size = 14) +

  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),

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

    legend.title = element_text(
      size = 12,
      face = "bold"
    ),

    legend.text = element_text(
      size = 11
    ),

    plot.title = element_text(
      size = 18,
      face = "bold"
    ),

    plot.subtitle = element_text(
      size = 10.5
    ),

    plot.margin = margin(
      10,
      14,
      10,
      10
    )
  )

print(p)


# ------------------------------------------------------------------------------
# 10. Save plots and tables
# ------------------------------------------------------------------------------

ggsave(
  filename = output_files$plot_png,
  plot = p,
  width = 8.4,
  height = 7.4,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = output_files$plot_pdf,
  plot = p,
  width = 8.4,
  height = 7.4,
  bg = "white"
)

write.csv(
  association_results,
  output_files$overall_csv,
  row.names = FALSE
)

write.csv(
  sensitivity_results,
  output_files$sensitivity_csv,
  row.names = FALSE
)

write.csv(
  df_matched,
  output_files$matched_csv,
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 11. Save Excel workbook
# ------------------------------------------------------------------------------

workbook <- createWorkbook()

addWorksheet(workbook, "Matched_samples")
writeData(workbook, "Matched_samples", df_matched)

addWorksheet(workbook, "Overall_matrix_results")
writeData(workbook, "Overall_matrix_results", association_results)

addWorksheet(workbook, "Cohort_sensitivity")
writeData(workbook, "Cohort_sensitivity", sensitivity_results)

addWorksheet(workbook, "Cohort_counts")
writeData(workbook, "Cohort_counts", cohort_counts)

addWorksheet(workbook, "Sensitivity_counts")
writeData(workbook, "Sensitivity_counts", sensitivity_counts)

addWorksheet(workbook, "Notes")
writeData(
  workbook,
  "Notes",
  data.frame(
    item = c(
      "Analysis",
      "Variables",
      "Matched dataset",
      "Overall FDR",
      "Sensitivity FDR",
      "Plot colour",
      "Plot text",
      "Plot asterisks"
    ),
    note = c(
      "All pairwise tests use two-sided Spearman correlation.",
      paste(matrix_vars, collapse = ", "),
      "Samples must be complete for cohort and every matrix variable.",
      "Overall raw P-values are adjusted together using Benjamini-Hochberg.",
      paste0(
        "Sensitivity P-values are adjusted within each cohort stratum ",
        "and across all sensitivity tests."
      ),
      "Tile colour represents raw P-value using -log10(P) for scaling.",
      "Black numbers inside tiles represent Spearman rho.",
      "* raw P < 0.05; ** raw P < 0.01; *** raw P < 0.001."
    )
  )
)

saveWorkbook(
  workbook,
  output_files$workbook,
  overwrite = TRUE
)


# ------------------------------------------------------------------------------
# 12. Save reproducibility information
# ------------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = output_files$session_info
)


# ------------------------------------------------------------------------------
# 13. Completion messages
# ------------------------------------------------------------------------------

message("\nAnalysis complete.")
message("Output directory: ", output_dir)
message("PNG: ", output_files$plot_png)
message("PDF: ", output_files$plot_pdf)
message("Excel workbook: ", output_files$workbook)
message("Overall results CSV: ", output_files$overall_csv)
message("Sensitivity results CSV: ", output_files$sensitivity_csv)
message("Matched samples CSV: ", output_files$matched_csv)
message("Session information: ", output_files$session_info)
