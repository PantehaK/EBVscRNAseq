#!/usr/bin/env Rscript

# ============================================================
# Age- and sex-adjusted Gamma GLMs + forest plot
# ============================================================
#
# Purpose:
#   This script tests age- and sex-adjusted associations between
#   circulating CNS injury biomarkers and EBV/MS-related predictors.
#
# Input:
#   correlation_variables.xlsx also available in source_data.xlsx sheet FigS13B-C
#
# Fully matched samples only:
#   cohort, sex, MS_signature, EBNA1_IgG, VCA_IgG, EA_IgG,
#   age, sNfL, sGFAP
#
# Models:
#   sNfL  ~ exposure_z + age_z + sex
#   sGFAP ~ exposure_z + age_z + sex
#
# Exposures:
#   MS_signature
#   EBNA1_IgG
#   VCA_IgG
#   EA_IgG
#
# Interpretation:
#   Continuous exposures are z-scored.
#   Effect ratios represent multiplicative change in sNfL or sGFAP
#   per 1 SD increase in exposure, adjusted for age and sex.
#
# Outputs:
#   - GLM_age_sex_results_table.csv
#   - GLM_with_covariates_results.xlsx
#   - GLM_results_NfL_GFAP_forest_plot_age_sex_rawP.png
#   - GLM_results_NfL_GFAP_forest_plot_age_sex_rawP.pdf
#
# ============================================================


# -----------------------------
# 0. Libraries
# -----------------------------

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(broom)
  library(ggplot2)
  library(openxlsx)
  library(tidyr)
  library(scales)
})


# -----------------------------
# 1. User settings
# -----------------------------

file_path <- "path/to/correlation_variables.xlsx"
out_dir <- dirname(file_path)

out_table_xlsx <- file.path(
  out_dir,
  "GLM_with_covariates_results.xlsx"
)

out_table_csv <- file.path(
  out_dir,
  "GLM_age_sex_results_table.csv"
)

out_forest_png <- file.path(
  out_dir,
  "GLM_results_NfL_GFAP_forest_plot_age_sex_rawP.png"
)

out_forest_pdf <- file.path(
  out_dir,
  "GLM_results_NfL_GFAP_forest_plot_age_sex_rawP.pdf"
)


# -----------------------------
# 2. Read data
# -----------------------------

df_raw <- read_xlsx(file_path) %>%
  as.data.frame()

names(df_raw) <- trimws(names(df_raw))

find_col_any <- function(possible_names) {
  
  hits <- names(df_raw)[tolower(names(df_raw)) %in% tolower(possible_names)]
  
  if (length(hits) == 0) {
    stop(
      "Missing required column. Tried: ",
      paste(possible_names, collapse = ", ")
    )
  }
  
  hits[1]
}

df_raw <- df_raw %>%
  transmute(
    sample = .data[[find_col_any(c("sample", "Sample"))]],
    
    cohort = .data[[find_col_any(c(
      "cohort",
      "diagnosis",
      "group"
    ))]],
    
    sex = .data[[find_col_any(c(
      "sex",
      "Sex",
      "gender",
      "Gender",
      "biological sex",
      "Biological sex"
    ))]],
    
    MS_signature = .data[[find_col_any(c(
      "MS_signature",
      "MS signature",
      "MS_sig",
      "MS score"
    ))]],
    
    EBNA1_IgG = .data[[find_col_any(c(
      "EBNA1_IgG",
      "EBNA1 IgG",
      "EBNA1_titre",
      "EBNA1 titre"
    ))]],
    
    VCA_IgG = .data[[find_col_any(c(
      "VCA_IgG",
      "VCA IgG",
      "VCA_titre",
      "VCA titre"
    ))]],
    
    EA_IgG = .data[[find_col_any(c(
      "EA_IgG",
      "EA IgG",
      "EA_titre",
      "EA titre"
    ))]],
    
    age = .data[[find_col_any(c(
      "age",
      "Age"
    ))]],
    
    sNfL = .data[[find_col_any(c(
      "sNfL",
      "sNFL",
      "NfL",
      "NFL"
    ))]],
    
    sGFAP = .data[[find_col_any(c(
      "sGFAP",
      "GFAP"
    ))]]
  )


# -----------------------------
# 3. Helper functions
# -----------------------------

to_continuous <- function(x) {
  
  z <- trimws(as.character(x))
  z[z %in% c("", "NA", "N/A", "na", "n/a", "NaN", "nan", "null")] <- NA
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
  z[z %in% c("", "NA", "na", "n/a", "nan", "null")] <- NA
  
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
    "ctrl",
    "healthy control",
    "hc",
    "nms"
  )] <- 0
  
  out
}

to_sex_factor <- function(x) {
  
  z <- tolower(trimws(as.character(x)))
  z[z %in% c("", "NA", "na", "n/a", "nan", "null")] <- NA
  
  out <- rep(NA_character_, length(z))
  
  out[z %in% c("f", "female", "woman")] <- "Female"
  out[z %in% c("m", "male", "man")] <- "Male"
  
  # Handles common numeric encodings.
  # Edit this if your metadata uses a different known coding scheme.
  out[z %in% c("0", "2")] <- "Female"
  out[z %in% c("1")] <- "Male"
  
  factor(out, levels = c("Female", "Male"))
}

format_p <- function(p) {
  
  case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}


# -----------------------------
# 4. Clean data
# -----------------------------

df <- df_raw %>%
  mutate(
    sample = as.character(sample),
    cohort = to_binary_cohort(cohort),
    sex_factor = to_sex_factor(sex),
    MS_signature = to_continuous(MS_signature),
    EBNA1_IgG = to_continuous(EBNA1_IgG),
    VCA_IgG = to_continuous(VCA_IgG),
    EA_IgG = to_continuous(EA_IgG),
    age = to_continuous(age),
    sNfL = to_continuous(sNfL),
    sGFAP = to_continuous(sGFAP)
  )


# -----------------------------
# 5. Fully matched samples only
# -----------------------------

matched_vars <- c(
  "cohort",
  "sex_factor",
  "MS_signature",
  "EBNA1_IgG",
  "VCA_IgG",
  "EA_IgG",
  "age",
  "sNfL",
  "sGFAP"
)

df_matched <- df %>%
  filter(if_all(all_of(matched_vars), ~ !is.na(.x))) %>%
  filter(
    cohort %in% c(0, 1),
    EBNA1_IgG >= 0,
    VCA_IgG >= 0,
    EA_IgG >= 0,
    sNfL > 0,
    sGFAP > 0
  ) %>%
  mutate(
    cohort_label = case_when(
      cohort == 1 ~ "MS",
      cohort == 0 ~ "Control"
    ),
    cohort_label = factor(cohort_label, levels = c("Control", "MS")),
    
    MS_signature_z = as.numeric(scale(MS_signature)),
    EBNA1_IgG_z = as.numeric(scale(EBNA1_IgG)),
    VCA_IgG_z = as.numeric(scale(VCA_IgG)),
    EA_IgG_z = as.numeric(scale(EA_IgG)),
    age_z = as.numeric(scale(age))
  )

message("Original samples: ", nrow(df))
message("Fully matched samples used: ", nrow(df_matched))

if (nrow(df_matched) < 6) {
  stop("Too few fully matched samples remain.")
}

if (length(unique(df_matched$sex_factor)) < 2) {
  stop("Sex has fewer than two levels in df_matched. Cannot include sex as a covariate.")
}

cohort_counts <- df_matched %>%
  count(cohort_label, name = "n")

sex_counts <- df_matched %>%
  count(sex_factor, name = "n")

cohort_sex_counts <- df_matched %>%
  count(cohort_label, sex_factor, name = "n")

print(cohort_counts)
print(sex_counts)
print(cohort_sex_counts)


# ============================================================
# PART A: Age- and sex-adjusted Gamma GLMs
# ============================================================

fit_gamma_age_sex_adjusted <- function(outcome, exposure_z, exposure_label) {
  
  needed_cols <- c(outcome, exposure_z, "age_z", "sex_factor")
  missing_cols <- setdiff(needed_cols, names(df_matched))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing column(s) in df_matched: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  model_df <- df_matched %>%
    select(all_of(needed_cols)) %>%
    filter(
      complete.cases(.),
      .data[[outcome]] > 0
    )
  
  if (nrow(model_df) < 6) {
    stop(
      "Too few complete observations for model: ",
      outcome,
      " ~ ",
      exposure_z,
      " + age_z + sex_factor"
    )
  }
  
  if (length(unique(model_df$sex_factor)) < 2) {
    stop(
      "Sex has fewer than two levels in the model data for ",
      outcome,
      " ~ ",
      exposure_z,
      ". Cannot include sex as a covariate."
    )
  }
  
  if (length(unique(model_df[[exposure_z]])) < 2) {
    stop("Exposure has no variation in the model data: ", exposure_z)
  }
  
  form <- as.formula(
    paste0(outcome, " ~ ", exposure_z, " + age_z + sex_factor")
  )
  
  mod <- glm(
    form,
    data = model_df,
    family = Gamma(link = "log")
  )
  
  broom::tidy(mod) %>%
    filter(term == exposure_z) %>%
    mutate(
      outcome = outcome,
      exposure = exposure_label,
      model = "Gamma GLM with log link",
      adjustment = "Age + sex",
      n = nrow(model_df),
      
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      
      effect_ratio = exp(estimate),
      effect_ratio_low = exp(conf.low),
      effect_ratio_high = exp(conf.high),
      
      p_label = format_p(p.value)
    )
}

glm_results <- bind_rows(
  fit_gamma_age_sex_adjusted(
    outcome = "sNfL",
    exposure_z = "MS_signature_z",
    exposure_label = "MS signature"
  ),
  fit_gamma_age_sex_adjusted(
    outcome = "sNfL",
    exposure_z = "EBNA1_IgG_z",
    exposure_label = "EBNA1 IgG"
  ),
  fit_gamma_age_sex_adjusted(
    outcome = "sNfL",
    exposure_z = "VCA_IgG_z",
    exposure_label = "VCA IgG"
  ),
  fit_gamma_age_sex_adjusted(
    outcome = "sNfL",
    exposure_z = "EA_IgG_z",
    exposure_label = "EA IgG"
  ),
  fit_gamma_age_sex_adjusted(
    outcome = "sGFAP",
    exposure_z = "MS_signature_z",
    exposure_label = "MS signature"
  ),
  fit_gamma_age_sex_adjusted(
    outcome = "sGFAP",
    exposure_z = "EBNA1_IgG_z",
    exposure_label = "EBNA1 IgG"
  ),
  fit_gamma_age_sex_adjusted(
    outcome = "sGFAP",
    exposure_z = "VCA_IgG_z",
    exposure_label = "VCA IgG"
  ),
  fit_gamma_age_sex_adjusted(
    outcome = "sGFAP",
    exposure_z = "EA_IgG_z",
    exposure_label = "EA IgG"
  )
) %>%
  group_by(outcome) %>%
  mutate(
    FDR_within_outcome = p.adjust(p.value, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    FDR_all_GLMs = p.adjust(p.value, method = "BH")
  )

publication_table <- glm_results %>%
  transmute(
    Outcome = outcome,
    Exposure = exposure,
    Model = model,
    Adjustment = adjustment,
    n = n,
    `Effect ratio` = sprintf("%.2f", effect_ratio),
    `95% CI` = paste0(
      sprintf("%.2f", effect_ratio_low),
      "",
      sprintf("%.2f", effect_ratio_high)
    ),
    `Raw P-value` = format_p(p.value),
    `FDR P-value within outcome` = format_p(FDR_within_outcome),
    `FDR P-value all GLMs` = format_p(FDR_all_GLMs)
  )

print(publication_table)


# ============================================================
# PART B: Save GLM tables
# ============================================================

write.csv(
  publication_table,
  out_table_csv,
  row.names = FALSE
)

wb <- createWorkbook()

addWorksheet(wb, "Publication_GLM_table")
writeData(wb, "Publication_GLM_table", publication_table)

addWorksheet(wb, "Raw_GLM_results")
writeData(wb, "Raw_GLM_results", glm_results)

addWorksheet(wb, "Matched_samples")
writeData(wb, "Matched_samples", df_matched)

addWorksheet(wb, "Cohort_counts")
writeData(wb, "Cohort_counts", cohort_counts)

addWorksheet(wb, "Sex_counts")
writeData(wb, "Sex_counts", sex_counts)

addWorksheet(wb, "Cohort_sex_counts")
writeData(wb, "Cohort_sex_counts", cohort_sex_counts)

addWorksheet(wb, "Notes")
writeData(
  wb,
  "Notes",
  data.frame(
    item = c(
      "GLM exposure scaling",
      "Gamma GLM adjustment",
      "Gamma GLM interpretation",
      "Sex coding"
    ),
    note = c(
      "MS_signature, EBNA1_IgG, VCA_IgG and EA_IgG are z-scored for GLM effect-ratio comparability.",
      "Gamma GLMs are adjusted for age and sex.",
      "Effect ratios represent multiplicative change in sNfL or sGFAP per 1 SD increase in exposure, adjusted for age and sex.",
      "Sex is converted to a two-level factor named sex_factor. Check the Sex_counts and Cohort_sex_counts sheets to confirm coding."
    )
  )
)

saveWorkbook(wb, out_table_xlsx, overwrite = TRUE)

message("Saved GLM table XLSX to: ", out_table_xlsx)
message("Saved GLM table CSV to: ", out_table_csv)


# ============================================================
# PART C: Forest plot
# ============================================================

plot_df <- glm_results %>%
  transmute(
    Outcome = outcome,
    Exposure = exposure,
    effect_ratio = effect_ratio,
    ci_low = effect_ratio_low,
    ci_high = effect_ratio_high,
    raw_p = p.value
  ) %>%
  mutate(
    Outcome = factor(Outcome, levels = c("sNfL", "sGFAP")),
    Exposure = factor(
      Exposure,
      levels = rev(c("MS signature", "EBNA1 IgG", "VCA IgG", "EA IgG"))
    ),
    
    p_label = paste0("P = ", format_p(raw_p)),
    
    effect_label = paste0(
      sprintf("%.2f", effect_ratio),
      " (",
      sprintf("%.2f", ci_low),
      "",
      sprintf("%.2f", ci_high),
      ")"
    )
  )

exposure_cols <- c(
  "MS signature" = "#7B3294",
  "EBNA1 IgG"    = "#0072B2",
  "VCA IgG"      = "#009E73",
  "EA IgG"       = "#fec832"
)

x_data_min <- min(plot_df$ci_low, na.rm = TRUE)
x_data_max <- max(plot_df$ci_high, na.rm = TRUE)
x_range <- x_data_max - x_data_min

if (!is.finite(x_range) || x_range <= 0) {
  x_range <- 1
}

effect_col_x <- max(1.68, x_data_max + 0.15 * x_range)
p_col_x <- effect_col_x + 0.42 * x_range

x_axis_min <- min(0.65, x_data_min - 0.08 * x_range)
x_axis_max <- p_col_x + 0.55 * x_range

p_forest <- ggplot(plot_df, aes(y = Exposure, colour = Exposure)) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.6,
    colour = "grey40"
  ) +
  
  geom_segment(
    aes(x = ci_low, xend = ci_high, yend = Exposure),
    linewidth = 1.1
  ) +
  
  geom_point(
    aes(x = effect_ratio),
    size = 3.4
  ) +
  
  geom_text(
    aes(x = effect_col_x, label = effect_label),
    colour = "black",
    hjust = 0,
    size = 3.8
  ) +
  
  geom_text(
    aes(x = p_col_x, label = p_label),
    colour = "black",
    hjust = 0,
    size = 3.8
  ) +
  
  facet_grid(
    Outcome ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  
  scale_colour_manual(values = exposure_cols) +
  
  scale_x_continuous(
    limits = c(x_axis_min, x_axis_max),
    breaks = scales::pretty_breaks(n = 5),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  labs(
    x = "Effect ratio",
    y = NULL,
    title = "Age- and sex-adjusted associations with circulating CNS injury biomarkers"
  ) +
  
  annotate(
    "text",
    x = effect_col_x,
    y = Inf,
    label = "Effect ratio (95% CI)",
    hjust = 0,
    vjust = 1.2,
    size = 3.8,
    fontface = "bold"
  ) +
  
  annotate(
    "text",
    x = p_col_x,
    y = Inf,
    label = "P-value",
    hjust = 0,
    vjust = 1.2,
    size = 3.8,
    fontface = "bold"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.background = element_blank(),
    strip.text.y = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    panel.spacing.y = unit(1.2, "lines"),
    plot.margin = margin(10, 35, 10, 10)
  )

print(p_forest)

ggsave(
  filename = out_forest_png,
  plot = p_forest,
  width = 12,
  height = 5,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = out_forest_pdf,
  plot = p_forest,
  width = 12,
  height = 5,
  bg = "white"
)

message("Saved forest plot PNG to: ", out_forest_png)
message("Saved forest plot PDF to: ", out_forest_pdf)
message("GLM and forest plot workflow complete.")
