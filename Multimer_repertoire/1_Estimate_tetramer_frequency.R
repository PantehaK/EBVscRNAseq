#!/usr/bin/env Rscript

# ==============================================================================
# Baseline EBV TCR metadata export, CD8 proportion estimation, tetramer plot,
# and MS-vs-Control frequency statistics
# ==============================================================================
#
# Purpose:
#   This script:
#     1. Exports metadata from the final baseline EBV tetramer-enriched Seurat object.
#     2. Counts detected tetramer-specific cells per GEMEBV sample.
#     3. Estimates CD8+ T cell proportions from global PBMC CITE-seq/ADT data.
#     4. Plots corrected EBV tetramer frequencies from a supplied source-data table.
#     5. Runs Wilcoxon rank-sum tests comparing EBV multimer frequencies between
#        MS and Control groups by:
#          - latency class, and
#          - antigen.
#     6. Applies BH/FDR correction and saves all outputs.
#
# Notes:
#   - Paths are intentionally generic for GitHub.
#   - Update the config block before running.
#   - Tetramer counts are detected-only counts from the GEMEBV object.
#   - Corrected tetramer frequencies are plotted from Source_data.xlsx Fig 1B.
#   - In the manuscript/source-data workflow, corrected frequencies were calculated as:
#
#       (tetramer cells / aligned cells) *
#       (sequenced cells / (CD8 proportion * total PBMC count)) * 100
#
#   - HLA-matched but undetected tetramers were assigned 1e-5 for plotting on a
#     log10 scale.
#   - Participant HLA information can be found in Table S1.
#   - PBMC count, number of cells sequenced and aligned can be found in
#     Table S15.
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
  library(readxl)
  library(readr)
  library(ggplot2)
  library(scales)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  
  # Input RDS files
  baseline_ebv_rds = "path/to/Publication_RDS/Baseline_EBV_integrated_seurat.rds",
  global_pbmc_rds = "path/to/Global_PBMC_RDS/9_merged_all_objects_final_filtered.rds",
  
  # Output directory
  out_dir = "path/to/Publication_data/EBV_baseline/Tetramer_frequency",
  
  # Source data file containing final corrected tetramer frequencies for plotting.
  # Expected structure by default:
  #   column 2 = cohort
  #   column 3 = epitope
  #   column 6 = frequency
  source_data_file = "path/to/Publication_data/Source_data.xlsx",
  source_data_sheet = 1,
  source_data_cohort_col_index = 2,
  source_data_epitope_col_index = 3,
  source_data_frequency_col_index = 6,
  
  # Input Excel file for MS-vs-Control frequency statistics. This is available in Supplementary_inputs_for_scripts on Zenodo.org.
  # Sheet 1 should contain Latency, Cohort, Frequency.
  # Sheet 2 should contain Antigen, Cohort, Frequency.
  multimer_frequency_stats_file = "path/to/Baseline_EBV_multimer_frequency_input.xlsx",
  latency_stats_sheet = 1,
  antigen_stats_sheet = 2,
  
  # Assay and ADT feature names in the global PBMC object
  adt_assay = "ADT",
  adt_cd8_feature = "CD8",
  adt_cd3_feature = "CD3",
  adt_cd4_feature = "CD4.1",
  
  # ADT thresholds used to classify CD8+ T cells
  cd3_threshold = 3,
  cd8_threshold = 2.5,
  cd4_threshold = 2.5,
  
  # Metadata columns
  sample_col = "sample",
  cohort_col = "cohort",
  diagnosis_col = "diagnosis",
  batch_col = "batch",
  tetramer_col = "tetramer",
  
  # Statistics settings
  stats_min_n_per_cohort = 2,
  stats_control_label = "Control",
  stats_ms_label = "MS",
  p_adjust_method = "BH",
  
  # Plot settings
  frequency_threshold = 1e-5,
  frequency_plot_y_limits = c(1e-5, 10),
  frequency_plot_width = 12,
  frequency_plot_height = 4,
  frequency_plot_dpi = 600,
  
  session_info_file = "path/to/Publication_data/EBV_baseline/Tetramer_frequency/sessionInfo_baseline_EBV_tetramer_frequency_and_stats.txt"
)


# ----------------------------- #
# 2. Tetramer order and colours
# ----------------------------- #

tetramers_with_hla <- c(
  "A02_GLCT",
  "A11_ATIG",
  "A24_TYPV",
  "B08_RAKF",
  "B35_EPLP*",
  "B40_VEDL",
  "B07_RPPI",
  "B08_FLRG",
  "B08_QAKW",
  "B35_YPLH*",
  "B40_LEKA",
  "A11_AVFD*",
  "A24_TYSA",
  "A02_LLDF",
  "A02_YLQQ",
  "A02_CLGG",
  "A02_FLYA",
  "A11_SSCS*",
  "A24_PYLF",
  "A24_TYGP",
  "B35_MGSL",
  "B40_IEDP",
  "B07_RPQK*",
  "B08_YNLR*",
  "B35_HPVG",
  "B40_FENI"
)

tetramers_short <- str_remove(tetramers_with_hla, "^[^_]+_")

tetramer_lookup <- tibble(
  tetramer_hla = tetramers_with_hla,
  tetramer_clean = tetramers_short
)

hla_cols <- c(
  "A02" = "#0072B2",
  "A11" = "#D55E00",
  "A24" = "#009E73",
  "B07" = "#CC79A7",
  "B08" = "#E69F00",
  "B35" = "#56B4E9",
  "B40" = "#F0E442"
)

# Samples used for CD8 proportion estimation from the global PBMC object.
# Update this vector if your manuscript/sample set changes.
samples_keep <- c(
  "MS173", "C263", "C287", "C288", "GR051", "C307", "MS359", "MS361",
  "C246", "MS069", "C111", "C339", "MS045", "MS362", "C233", "MS352",
  "MS389", "C049", "MS394", "C343", "MS351", "MS438", "MS387", "MS476",
  "MS136", "MS138", "MS358", "C073", "C083", "C126", "C219", "MS153",
  "MS034", "MS355", "MS461", "C134"
)


# ----------------------------- #
# 3. Helper functions
# ----------------------------- #

create_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}


create_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}


load_seurat_object <- function(path, object_label = "Seurat object") {
  
  if (!file.exists(path)) {
    stop(object_label, " RDS file does not exist: ", path)
  }
  
  obj <- readRDS(path)
  
  if (!inherits(obj, "Seurat")) {
    stop(object_label, " RDS does not contain a Seurat object: ", path)
  }
  
  obj
}


get_assay_data_compat <- function(object, assay, data_layer = "data") {
  
  # Compatible with Seurat v4 and Seurat v5.
  tryCatch(
    {
      GetAssayData(object, assay = assay, layer = data_layer)
    },
    error = function(e) {
      GetAssayData(object, assay = assay, slot = data_layer)
    }
  )
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


clean_cohort <- function(x, control_label = "Control", ms_label = "MS") {
  
  x <- str_to_lower(str_trim(as.character(x)))
  
  case_when(
    x %in% c("control", "controls", "healthy control", "healthy controls", "hc", "nms", "non-ms", "non_ms") ~ control_label,
    x %in% c("ms", "pwms", "multiple sclerosis", "rrms", "ppms", "spms") ~ ms_label,
    TRUE ~ as.character(x)
  )
}


clean_tetramer_name <- function(x) {
  
  x %>%
    as.character() %>%
    str_trim() %>%
    str_remove("^[^_]+_")
}


safe_wilcox_p <- function(values, groups, control_label, ms_label) {
  
  test_df <- tibble(
    value = suppressWarnings(as.numeric(values)),
    group = as.character(groups)
  ) %>%
    filter(
      !is.na(value),
      group %in% c(control_label, ms_label)
    )
  
  if (n_distinct(test_df$group) < 2) {
    return(NA_real_)
  }
  
  if (nrow(test_df) < 2) {
    return(NA_real_)
  }
  
  if (length(unique(test_df$value)) < 2) {
    return(1)
  }
  
  suppressWarnings(
    wilcox.test(value ~ group, data = test_df, exact = FALSE)$p.value
  )
}


run_grouped_frequency_stats <- function(
    df,
    group_col,
    cohort_col = "Cohort",
    frequency_col = "Frequency",
    control_label = "Control",
    ms_label = "MS",
    min_n_per_cohort = 2,
    p_adjust_method = "BH"
) {
  
  required_cols <- c(group_col, cohort_col, frequency_col)
  check_required_columns(df, required_cols, object_name = "frequency stats input")
  
  df_clean <- df %>%
    transmute(
      group = as.character(.data[[group_col]]),
      Cohort = clean_cohort(
        .data[[cohort_col]],
        control_label = control_label,
        ms_label = ms_label
      ),
      Frequency = suppressWarnings(as.numeric(.data[[frequency_col]]))
    ) %>%
    filter(
      !is.na(group),
      group != "",
      Cohort %in% c(control_label, ms_label),
      !is.na(Frequency)
    )
  
  stats_df <- df_clean %>%
    group_by(group) %>%
    summarise(
      n_MS = sum(Cohort == ms_label),
      n_Control = sum(Cohort == control_label),
      median_MS = median(Frequency[Cohort == ms_label], na.rm = TRUE),
      median_Control = median(Frequency[Cohort == control_label], na.rm = TRUE),
      diff_median = median_MS - median_Control,
      fold_change = median_MS / (median_Control + 1e-9),
      p_value = safe_wilcox_p(
        values = Frequency,
        groups = Cohort,
        control_label = control_label,
        ms_label = ms_label
      ),
      .groups = "drop"
    ) %>%
    mutate(
      tested = n_MS >= min_n_per_cohort & n_Control >= min_n_per_cohort,
      p_value = if_else(tested, p_value, NA_real_),
      p_adj = p.adjust(p_value, method = p_adjust_method)
    ) %>%
    rename(!!group_col := group) %>%
    arrange(p_adj, p_value)
  
  stats_df
}


# ----------------------------- #
# 4. Load baseline EBV object and export metadata
# ----------------------------- #

create_dir(config$out_dir)
create_parent_dir(config$session_info_file)

message("Loading baseline EBV object...")
baseline_obj <- load_seurat_object(
  config$baseline_ebv_rds,
  object_label = "Baseline EBV"
)

message("Metadata columns in baseline EBV object:")
print(colnames(baseline_obj@meta.data))

cols_to_export <- c(
  "id",
  "tetramer_CLR",
  "cohort",
  "diagnosis",
  "age",
  "sex",
  "batch",
  "years_since_diagnosis",
  "sample",
  "TRA_v_gene",
  "TRA_d_gene",
  "TRA_j_gene",
  "TRA_c_gene",
  "TRA_cdr3",
  "TRA_cdr3_nt",
  "TRB_v_gene",
  "TRB_d_gene",
  "TRB_j_gene",
  "TRB_c_gene",
  "TRB_cdr3",
  "TRB_cdr3_nt",
  "tetramer",
  "virus",
  "lifecycle",
  "latency",
  "antigen",
  "new_cluster",
  "celltype_new"
)

required_baseline_cols <- c(
  config$sample_col,
  config$cohort_col,
  config$diagnosis_col,
  config$batch_col,
  config$tetramer_col
)

check_required_columns(
  baseline_obj@meta.data,
  required_baseline_cols,
  object_name = "baseline EBV metadata"
)

missing_optional_cols <- setdiff(cols_to_export, colnames(baseline_obj@meta.data))

if (length(missing_optional_cols) > 0) {
  message("Optional export columns not found and therefore skipped:")
  print(missing_optional_cols)
}

meta_export <- baseline_obj@meta.data %>%
  tibble::rownames_to_column("cell_barcode") %>%
  select(cell_barcode, any_of(cols_to_export))

metadata_export_file <- file.path(
  config$out_dir,
  "Baseline_EBV_TCR_metadata_export.csv"
)

write_csv(meta_export, metadata_export_file)

message("Saved metadata export: ", metadata_export_file)


# ----------------------------- #
# 5. Count detected tetramer-specific cells per GEMEBV sample
# ----------------------------- #

gemebv_meta <- meta_export %>%
  mutate(
    batch = as.character(.data[[config$batch_col]]),
    tetramer = as.character(.data[[config$tetramer_col]]),
    tetramer_clean = clean_tetramer_name(tetramer)
  ) %>%
  filter(str_starts(batch, "GEMEBV"))

message("Number of cells from batches starting with GEMEBV:")
print(nrow(gemebv_meta))

message("GEMEBV batches found:")
print(sort(unique(gemebv_meta$batch)))

message("Tetramer values found in GEMEBV metadata:")
print(sort(unique(na.omit(gemebv_meta$tetramer_clean))))

message("Expected tetramers not detected in GEMEBV metadata:")
print(setdiff(tetramers_short, unique(na.omit(gemebv_meta$tetramer_clean))))

tetramer_counts_long <- gemebv_meta %>%
  filter(tetramer_clean %in% tetramers_short) %>%
  count(
    sample,
    cohort,
    diagnosis,
    batch,
    tetramer_clean,
    name = "n_tetramer_cells"
  ) %>%
  left_join(tetramer_lookup, by = "tetramer_clean") %>%
  mutate(
    tetramer_clean = factor(tetramer_clean, levels = tetramers_short),
    tetramer_hla = factor(tetramer_hla, levels = tetramers_with_hla)
  ) %>%
  arrange(sample, tetramer_hla)

tetramer_counts_per_sample <- tetramer_counts_long %>%
  group_by(sample, cohort, diagnosis, batch) %>%
  summarise(
    total_tetramer_specific_cells = sum(n_tetramer_cells),
    n_detected_tetramers = n_distinct(tetramer_clean),
    .groups = "drop"
  ) %>%
  arrange(sample)

tetramer_counts_long_file <- file.path(
  config$out_dir,
  "GEMEBV_tetramer_cell_counts_long_detected_only.csv"
)

tetramer_counts_per_sample_file <- file.path(
  config$out_dir,
  "GEMEBV_total_tetramer_specific_cells_per_sample.csv"
)

write_csv(tetramer_counts_long, tetramer_counts_long_file)
write_csv(tetramer_counts_per_sample, tetramer_counts_per_sample_file)

message("Saved GEMEBV tetramer count outputs.")


# ----------------------------- #
# 6. Estimate CD8 proportions from global PBMC ADT data
# ----------------------------- #

message("Loading global PBMC object...")
global_obj <- load_seurat_object(
  config$global_pbmc_rds,
  object_label = "Global PBMC"
)

if (!config$adt_assay %in% Assays(global_obj)) {
  stop("ADT assay not found in global PBMC object: ", config$adt_assay)
}

check_required_columns(
  global_obj@meta.data,
  c(config$sample_col, config$batch_col),
  object_name = "global PBMC metadata"
)

DefaultAssay(global_obj) <- config$adt_assay

adt_mat <- get_assay_data_compat(
  object = global_obj,
  assay = config$adt_assay,
  data_layer = "data"
)

message("ADT features present:")
print(rownames(adt_mat))

adt_features <- c(
  config$adt_cd8_feature,
  config$adt_cd3_feature,
  config$adt_cd4_feature
)

missing_adt_features <- setdiff(adt_features, rownames(adt_mat))

if (length(missing_adt_features) > 0) {
  stop("Missing ADT features: ", paste(missing_adt_features, collapse = ", "))
}

global_obj$ADT_CD8 <- as.numeric(adt_mat[config$adt_cd8_feature, colnames(global_obj)])
global_obj$ADT_CD3 <- as.numeric(adt_mat[config$adt_cd3_feature, colnames(global_obj)])
global_obj$ADT_CD4_1 <- as.numeric(adt_mat[config$adt_cd4_feature, colnames(global_obj)])

message("ADT_CD8 summary:")
print(summary(global_obj$ADT_CD8))

message("ADT_CD3 summary:")
print(summary(global_obj$ADT_CD3))

message("ADT_CD4_1 summary:")
print(summary(global_obj$ADT_CD4_1))

cd8_scatter <- FeatureScatter(
  object = global_obj,
  feature1 = "ADT_CD8",
  feature2 = "ADT_CD4_1"
)

cd8_scatter_file <- file.path(config$out_dir, "ADT_CD8_vs_CD4_feature_scatter.png")

ggsave(
  filename = cd8_scatter_file,
  plot = cd8_scatter,
  width = 5,
  height = 4,
  dpi = 600
)

global_obj$CD8 <- ifelse(
  global_obj$ADT_CD3 > config$cd3_threshold &
    global_obj$ADT_CD8 > config$cd8_threshold &
    global_obj$ADT_CD4_1 < config$cd4_threshold,
  "YES",
  "NO"
)

message("CD8 classification table:")
print(table(global_obj$CD8, useNA = "ifany"))

cd8_props <- global_obj@meta.data %>%
  mutate(
    sample = as.character(.data[[config$sample_col]]),
    batch = as.character(.data[[config$batch_col]])
  ) %>%
  filter(
    sample %in% samples_keep,
    !str_detect(batch, "EBV")
  ) %>%
  group_by(sample) %>%
  summarise(
    total_cells = n(),
    cd8_yes = sum(CD8 == "YES", na.rm = TRUE),
    prop_cd8 = cd8_yes / total_cells,
    .groups = "drop"
  ) %>%
  arrange(desc(prop_cd8))

message("CD8 proportions:")
print(cd8_props)

cd8_props_file <- file.path(config$out_dir, "CD8_proportions_global_PBMC.csv")
write_csv(cd8_props, cd8_props_file)

message("Saved CD8 proportion output: ", cd8_props_file)


# ----------------------------- #
# 7. Plot corrected EBV tetramer frequencies from source data
# ----------------------------- #

if (!file.exists(config$source_data_file)) {
  stop(
    "Could not find source data file: ",
    config$source_data_file,
    "\nPlease update config$source_data_file."
  )
}

message("Loading source data file for plotting:")
message(config$source_data_file)

df_raw <- read_excel(
  config$source_data_file,
  sheet = config$source_data_sheet,
  .name_repair = "unique"
)

source_col_indices <- c(
  config$source_data_cohort_col_index,
  config$source_data_epitope_col_index,
  config$source_data_frequency_col_index
)

if (ncol(df_raw) < max(source_col_indices)) {
  stop(
    "Source data file has fewer columns than expected. Required column indices: ",
    paste(source_col_indices, collapse = ", ")
  )
}

df <- df_raw[, source_col_indices]
colnames(df) <- c("cohort", "epitope", "frequency")

df <- df %>%
  mutate(
    cohort = clean_cohort(
      cohort,
      control_label = config$stats_control_label,
      ms_label = config$stats_ms_label
    ),
    epitope = str_trim(as.character(epitope)),
    frequency = suppressWarnings(
      as.numeric(str_replace_all(as.character(frequency), "[,<]", ""))
    ),
    frequency_plot = case_when(
      is.na(frequency) ~ NA_real_,
      frequency <= 0 ~ config$frequency_threshold,
      TRUE ~ frequency
    ),
    detection_status = case_when(
      frequency_plot <= config$frequency_threshold ~ "Below threshold/non-detected",
      TRUE ~ "Detected"
    ),
    hla = sub("_.*$", "", epitope),
    hla = case_when(
      hla == "B7" ~ "B07",
      TRUE ~ hla
    ),
    epitope = factor(epitope, levels = tetramers_with_hla),
    cohort = factor(
      cohort,
      levels = c(config$stats_control_label, config$stats_ms_label)
    )
  ) %>%
  filter(
    !is.na(epitope),
    !is.na(cohort),
    !is.na(frequency_plot),
    frequency_plot > 0
  )

missing_epitopes <- setdiff(tetramers_with_hla, unique(as.character(df$epitope)))

message("Epitopes in desired order that are missing from plotting data:")
print(missing_epitopes)

plot_input_file <- file.path(
  config$out_dir,
  "EBV_tetramer_frequencies_plot_input_cleaned.csv"
)

write_csv(df, plot_input_file)

p_freq <- ggplot(df, aes(x = epitope, y = frequency_plot)) +
  
  geom_hline(
    yintercept = config$frequency_threshold,
    linetype = "dashed",
    linewidth = 0.3,
    colour = "grey40"
  ) +
  
  geom_point(
    aes(shape = cohort, fill = hla),
    position = position_jitter(width = 0.15, height = 0),
    size = 3,
    alpha = 0.95,
    colour = "black",
    stroke = 0.7
  ) +
  
  stat_summary(
    fun = median,
    geom = "crossbar",
    aes(group = epitope),
    width = 0.6,
    linewidth = 0.4,
    colour = "black",
    na.rm = TRUE
  ) +
  
  scale_y_log10(
    limits = config$frequency_plot_y_limits,
    breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10),
    labels = function(x) {
      sub("\\.?0+$", "", format(x, scientific = FALSE))
    }
  ) +
  
  scale_x_discrete(drop = FALSE) +
  
  scale_shape_manual(
    values = c(
      "Control" = 21,
      "MS" = 24
    ),
    drop = FALSE
  ) +
  
  scale_fill_manual(
    values = hla_cols,
    breaks = c("A02", "A11", "A24", "B07", "B08", "B35", "B40"),
    labels = c(
      "A02" = "HLA-A*02:01",
      "A11" = "HLA-A*11:01",
      "A24" = "HLA-A*24:02",
      "B07" = "HLA-B*07:02",
      "B08" = "HLA-B*08:01",
      "B35" = "HLA-B*35:01",
      "B40" = "HLA-B*40:01"
    ),
    drop = FALSE
  ) +
  
  guides(
    fill = guide_legend(
      title = "HLA",
      override.aes = list(
        shape = 21,
        colour = "black",
        size = 4,
        stroke = 0.7
      )
    ),
    shape = guide_legend(
      title = "Cohort",
      override.aes = list(
        fill = "white",
        colour = "black",
        size = 4,
        stroke = 0.7
      )
    )
  ) +
  
  labs(
    x = NULL,
    y = "Estimated tetramer-specific CD8+ T cells (%)",
    shape = "Cohort",
    fill = "HLA"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      colour = "black"
    ),
    axis.text.y = element_text(colour = "black"),
    axis.title.y = element_text(colour = "black"),
    axis.line = element_line(colour = "black"),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_freq)

frequency_plot_png <- file.path(config$out_dir, "EBV_tetramer_frequencies_colourblind.png")
frequency_plot_pdf <- file.path(config$out_dir, "EBV_tetramer_frequencies_colourblind.pdf")

ggsave(
  filename = frequency_plot_png,
  plot = p_freq,
  width = config$frequency_plot_width,
  height = config$frequency_plot_height,
  dpi = config$frequency_plot_dpi
)

ggsave(
  filename = frequency_plot_pdf,
  plot = p_freq,
  width = config$frequency_plot_width,
  height = config$frequency_plot_height
)

message("Saved EBV tetramer frequency plots.")


# ----------------------------- #
# 8. Run MS-vs-Control frequency statistics
# ----------------------------- #

if (!file.exists(config$multimer_frequency_stats_file)) {
  stop(
    "Could not find multimer frequency statistics input file: ",
    config$multimer_frequency_stats_file,
    "\nPlease update config$multimer_frequency_stats_file."
  )
}

message("Loading multimer frequency statistics input:")
message(config$multimer_frequency_stats_file)

latency_input <- read_excel(
  config$multimer_frequency_stats_file,
  sheet = config$latency_stats_sheet
)

antigen_input <- read_excel(
  config$multimer_frequency_stats_file,
  sheet = config$antigen_stats_sheet
)

latency_stats <- run_grouped_frequency_stats(
  df = latency_input,
  group_col = "Latency",
  cohort_col = "Cohort",
  frequency_col = "Frequency",
  control_label = config$stats_control_label,
  ms_label = config$stats_ms_label,
  min_n_per_cohort = config$stats_min_n_per_cohort,
  p_adjust_method = config$p_adjust_method
)

antigen_stats <- run_grouped_frequency_stats(
  df = antigen_input,
  group_col = "Antigen",
  cohort_col = "Cohort",
  frequency_col = "Frequency",
  control_label = config$stats_control_label,
  ms_label = config$stats_ms_label,
  min_n_per_cohort = config$stats_min_n_per_cohort,
  p_adjust_method = config$p_adjust_method
)

latency_stats_file <- file.path(
  config$out_dir,
  "Baseline_EBV_latency_frequencies_stats.csv"
)

antigen_stats_file <- file.path(
  config$out_dir,
  "Baseline_EBV_antigen_frequencies_stats.csv"
)

write_csv(latency_stats, latency_stats_file)
write_csv(antigen_stats, antigen_stats_file)

message("Saved latency frequency statistics: ", latency_stats_file)
message("Saved antigen frequency statistics: ", antigen_stats_file)


# ----------------------------- #
# 9. Save session information
# ----------------------------- #

writeLines(
  c(
    "Configuration:",
    capture.output(str(config)),
    "",
    "Metadata export:",
    metadata_export_file,
    "",
    "Tetramer count outputs:",
    tetramer_counts_long_file,
    tetramer_counts_per_sample_file,
    "",
    "CD8 proportion output:",
    cd8_props_file,
    "",
    "Frequency plot outputs:",
    frequency_plot_png,
    frequency_plot_pdf,
    "",
    "Frequency statistics outputs:",
    latency_stats_file,
    antigen_stats_file,
    "",
    "Latency statistics:",
    capture.output(print(latency_stats)),
    "",
    "Antigen statistics:",
    capture.output(print(antigen_stats)),
    "",
    "Session information:",
    capture.output(sessionInfo())
  ),
  config$session_info_file
)

message("\nBaseline EBV tetramer frequency and statistics workflow complete.")
message("Session info: ", config$session_info_file)
