#!/usr/bin/env Rscript

# ============================================================
# Spearman correlations + antibody IgG vs MS signature plots
# ============================================================
#
# Purpose:
#   This script tests Spearman correlations between raw EBV antibody IgG
#   values and the raw MS signature score.
#
# Input:
#   correlation_variables.xlsx Also available from source_data sheet Fig13B-C
#
# Fully matched samples only:
#   cohort, sex, MS_signature, EBNA1_IgG, VCA_IgG, EA_IgG,
#   age, sNfL, sGFAP
#
# Tests:
#   EBNA1_IgG vs MS_signature
#   VCA_IgG   vs MS_signature
#   EA_IgG    vs MS_signature
#
# Plots:
#   EBNA1_IgG vs MS_signature
#   VCA_IgG   vs MS_signature
#   EA_IgG    vs MS_signature
#
# No statistics are printed on the graphs.
#
# Outputs:
#   - Spearman_results.csv
#   - Spearman_results.xlsx
#   - EBNA1_IgG_MS_signature_scatter.png/pdf
#   - VCA_IgG_MS_signature_scatter.png/pdf
#   - EA_IgG_MS_signature_scatter.png/pdf
#
# ============================================================


# -----------------------------
# 0. Libraries
# -----------------------------

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(openxlsx)
  library(tidyr)
})


# -----------------------------
# 1. User settings
# -----------------------------

file_path <- "path/to/correlation_variables.xlsx"
out_dir <- dirname(file_path)

out_plot_ebna1_png <- file.path(
  out_dir,
  "EBNA1_IgG_MS_signature_scatter.png"
)

out_plot_ebna1_pdf <- file.path(
  out_dir,
  "EBNA1_IgG_MS_signature_scatter.pdf"
)

out_plot_vca_png <- file.path(
  out_dir,
  "VCA_IgG_MS_signature_scatter.png"
)

out_plot_vca_pdf <- file.path(
  out_dir,
  "VCA_IgG_MS_signature_scatter.pdf"
)

out_plot_ea_png <- file.path(
  out_dir,
  "EA_IgG_MS_signature_scatter.png"
)

out_plot_ea_pdf <- file.path(
  out_dir,
  "EA_IgG_MS_signature_scatter.pdf"
)

out_spearman_csv <- file.path(
  out_dir,
  "Spearman_results.csv"
)

out_spearman_xlsx <- file.path(
  out_dir,
  "Spearman_results.xlsx"
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

run_spearman <- function(x, y, comparison_label) {
  
  keep <- complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  
  n <- length(x)
  
  if (
    n < 3 ||
    length(unique(x)) < 2 ||
    length(unique(y)) < 2
  ) {
    return(
      data.frame(
        comparison = comparison_label,
        test = "Spearman correlation",
        n = n,
        rho = NA_real_,
        p.value = NA_real_,
        note = "Not tested: too few samples or no variation"
      )
    )
  }
  
  test <- suppressWarnings(
    cor.test(
      x,
      y,
      method = "spearman",
      exact = FALSE
    )
  )
  
  data.frame(
    comparison = comparison_label,
    test = "Spearman correlation",
    n = n,
    rho = unname(test$estimate),
    p.value = test$p.value,
    note = NA_character_
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
    cohort_label = factor(cohort_label, levels = c("Control", "MS"))
  )

message("Original samples: ", nrow(df))
message("Fully matched samples used: ", nrow(df_matched))

if (nrow(df_matched) < 6) {
  stop("Too few fully matched samples remain.")
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
# PART A: Spearman tests
# ============================================================

spearman_results <- bind_rows(
  run_spearman(
    df_matched$EBNA1_IgG,
    df_matched$MS_signature,
    "EBNA1_IgG vs MS_signature"
  ),
  run_spearman(
    df_matched$VCA_IgG,
    df_matched$MS_signature,
    "VCA_IgG vs MS_signature"
  ),
  run_spearman(
    df_matched$EA_IgG,
    df_matched$MS_signature,
    "EA_IgG vs MS_signature"
  )
) %>%
  mutate(
    FDR = p.adjust(p.value, method = "BH"),
    rho_label = ifelse(is.na(rho), "NA", sprintf("%.3f", rho)),
    p_label = format_p(p.value),
    FDR_label = format_p(FDR)
  )

print(spearman_results)


# ============================================================
# PART B: Publication-style antibody IgG vs MS signature plots
# ============================================================

cohort_cols <- c(
  "Control" = "#2E6F9E",
  "MS" = "#8B1E3F"
)

make_ms_signature_scatter <- function(xvar, xlab, out_png, out_pdf) {
  
  p <- ggplot(
    df_matched,
    aes(
      x = .data[[xvar]],
      y = MS_signature
    )
  ) +
    geom_point(
      aes(fill = cohort_label),
      shape = 21,
      size = 4.6,
      stroke = 0.75,
      colour = "grey15",
      alpha = 0.92
    ) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      linewidth = 0.95,
      colour = "grey25",
      linetype = "dashed"
    ) +
    scale_fill_manual(values = cohort_cols) +
    labs(
      x = xlab,
      y = "MS signature",
      fill = NULL
    ) +
    theme_classic(base_size = 16) +
    theme(
      axis.title = element_text(size = 17, face = "bold", colour = "black"),
      axis.text = element_text(size = 13, colour = "black"),
      axis.line = element_line(linewidth = 0.7, colour = "black"),
      axis.ticks = element_line(linewidth = 0.6, colour = "black"),
      
      legend.position = "top",
      legend.justification = "right",
      legend.text = element_text(size = 13),
      legend.key.size = unit(0.7, "cm"),
      
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(8, 10, 8, 8)
    )
  
  print(p)
  
  ggsave(
    filename = out_png,
    plot = p,
    width = 5.8,
    height = 5.2,
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    filename = out_pdf,
    plot = p,
    width = 5.8,
    height = 5.2,
    bg = "white"
  )
  
  return(p)
}

p_ebna1 <- make_ms_signature_scatter(
  xvar = "EBNA1_IgG",
  xlab = "EBNA1 IgG",
  out_png = out_plot_ebna1_png,
  out_pdf = out_plot_ebna1_pdf
)

p_vca <- make_ms_signature_scatter(
  xvar = "VCA_IgG",
  xlab = "VCA IgG",
  out_png = out_plot_vca_png,
  out_pdf = out_plot_vca_pdf
)

p_ea <- make_ms_signature_scatter(
  xvar = "EA_IgG",
  xlab = "EA IgG",
  out_png = out_plot_ea_png,
  out_pdf = out_plot_ea_pdf
)


# ============================================================
# PART C: Save Spearman outputs
# ============================================================

write.csv(
  spearman_results,
  out_spearman_csv,
  row.names = FALSE
)

wb <- createWorkbook()

addWorksheet(wb, "Spearman_results")
writeData(wb, "Spearman_results", spearman_results)

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
      "EBNA1 plot",
      "VCA plot",
      "EA plot",
      "Spearman tests"
    ),
    note = c(
      "The EBNA1 scatter plot uses raw EBNA1_IgG values on the x-axis and raw MS_signature values on the y-axis. No statistics are printed on the graph.",
      "The VCA scatter plot uses raw VCA_IgG values on the x-axis and raw MS_signature values on the y-axis. No statistics are printed on the graph.",
      "The EA scatter plot uses raw EA_IgG values on the x-axis and raw MS_signature values on the y-axis. No statistics are printed on the graph.",
      "Spearman correlations use raw EBNA1_IgG, raw VCA_IgG, raw EA_IgG and raw MS_signature values."
    )
  )
)

saveWorkbook(wb, out_spearman_xlsx, overwrite = TRUE)

message("Saved Spearman CSV to: ", out_spearman_csv)
message("Saved Spearman XLSX to: ", out_spearman_xlsx)

message("Saved EBNA1 scatter PNG to: ", out_plot_ebna1_png)
message("Saved EBNA1 scatter PDF to: ", out_plot_ebna1_pdf)

message("Saved VCA scatter PNG to: ", out_plot_vca_png)
message("Saved VCA scatter PDF to: ", out_plot_vca_pdf)

message("Saved EA scatter PNG to: ", out_plot_ea_png)
message("Saved EA scatter PDF to: ", out_plot_ea_pdf)

message("Spearman correlation workflow complete.")
