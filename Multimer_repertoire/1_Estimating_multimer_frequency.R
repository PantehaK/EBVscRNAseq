#!/usr/bin/env Rscript

# ==============================================================================
# Baseline EBV TCR export, CD8 proportion estimation, and tetramer frequencies
# ==============================================================================
#
# Purpose:
#   This script exports metadata from the module-scored baseline EBV T cell object,
#   estimates CD8+ T cell proportions from global PBMC CITE-seq/ADT data, and
#   calculates corrected baseline EBV tetramer frequencies.
#
#   It saves:
#     1. baseline EBV TCR/module-score metadata,
#     2. global PBMC CD8 proportion estimates,
#     3. per-sample tetramer frequency estimates,
#     4. HLA-expanded tetramer frequency table with undetected values filled,
#     5. detected versus non-detected responder summaries,
#     6. stacked bar plots of detected/non-detected proportions,
#     7. session information.
#
# Expected inputs:
#   - Module-scored baseline EBV object:
#       15_baseline_EBV_module_scored.rds
#
#   - Final merged all-object Seurat object with global PBMCs and ADT assay:
#       8_merged_all_objects_final.rds
#
#   - External sample-level frequency input table containing:
#       id, sample, cohort, Total PBMC count, Cells sequenced, Cells aligned
#     and optionally:
#       proportion of CD8
#
#   - HLA table for expected tetramer/sample combinations.
#
# Notes:
#   - The estimated tetramer percentage is calculated as:
#       (tetramer cells / aligned cells) *
#       (sequenced cells / (CD8 proportion * total PBMC count)) * 100
#
#   - Undetected HLA-matched tetramers are assigned a small value so they can be
#     retained in downstream responder/non-responder plots.
#
# ==============================================================================


# ----------------------------- #
# 0. Load packages
# ----------------------------- #

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(readxl)
  library(ggplot2)
  library(scales)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  dataset_id = "baseline_ebv",
  
  baseline_module_scored_rds = "path/to/input/15_baseline_EBV_module_scored.rds",
  merged_all_objects_rds = "path/to/input/8_merged_all_objects_final.rds",
  
  output_dir = "path/to/output/Publication_data/EBV_baseline/Tetramer_frequency",
  
  baseline_tcr_export_csv = "path/to/output/Publication_data/EBV_baseline/Tetramer_frequency/Baseline_EBV_TCR_module_score.csv",
  cd8_proportions_csv = "path/to/output/Publication_data/EBV_baseline/Tetramer_frequency/CD8_proportions.csv",
  responders_tetramer_freq_csv = "path/to/output/Publication_data/EBV_baseline/Tetramer_frequency/responders_tetramer_freq.csv",
  baseline_tetramer_frequencies_csv = "path/to/output/Publication_data/EBV_baseline/Tetramer_frequency/Baseline_EBV_tetramer_frequencies.csv",
  detected_summary_csv = "path/to/output/Publication_data/EBV_baseline/Tetramer_frequency/detected_vs_non_detected.csv",
  session_info_file = "path/to/output/Publication_data/EBV_baseline/Tetramer_frequency/sessionInfo_tetramer_frequency.txt",
  
  # External sample-level input file.
  # Required columns:
  #   id, sample, Total PBMC count, Cells sequenced, Cells aligned
  # Recommended:
  #   cohort, proportion of CD8
  frequency_input_xlsx = "path/to/input/single_cell_output.xlsx",
  
  # HLA reference input.
  # This can either be:
  #   1. long/wide sample matrix with sample column and allele columns, or
  #   2. original format where first column is allele and remaining cells contain samples.
  sample_hla_xlsx = "path/to/input/sample_hla.xlsx",
  
  # ADT CD8 proportion settings.
  adt_assay = "ADT",
  adt_layer = "data",
  adt_cd8_feature = "CD8",
  adt_cd3_feature = "CD3",
  adt_cd4_feature = "CD4.1",
  adt_cd3_min = 3,
  adt_cd8_min = 2.5,
  adt_cd4_max = 2.5,
  
  global_pbmc_exclude_batch_pattern = "EBV",
  
  samples_keep = c(
    "MS173", "C263", "C287", "C288", "GR051", "C307", "MS359", "MS361",
    "C246", "MS069", "C111", "C339", "MS045", "MS362", "C233", "MS352",
    "MS389", "C049", "MS394", "C343", "MS351", "MS438", "MS387", "MS476",
    "MS136", "MS138", "MS358", "C073", "C083", "C126", "C219", "MS153",
    "MS034", "MS355", "MS461", "C134"
  ),
  
  # Optional manual CD8 proportion overrides for samples not present in global PBMC.
  manual_cd8_overrides = tibble::tribble(
    ~sample, ~prop_cd8, ~note,
    "GR001", 0.30, "Manual in-house flow cytometry estimate"
  ),
  
  # Metadata export columns from baseline module-scored object.
  cols_to_export = c(
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
    "matched_cluster",
    "celltype",
    "Cytotoxicity_score"
  ),
  
  # Baseline EBV batch pattern used for tetramer counts.
  baseline_ebv_batch_pattern = "^GEMEBV",
  
  tetramers_to_exclude_from_frequency = c(
    NA,
    "",
    "negative",
    "enriched",
    "LCL-responsive"
  ),
  
  # If tetramers in the Seurat metadata are not HLA-prefixed, this lookup adds
  # the allele prefix needed for HLA-matched frequency completion.
  tetramer_prefix_lookup = tibble::tribble(
    ~tetramer, ~tetramer_prefixed,
    "GLCT",  "A02_GLCT",
    "CLGG",  "A02_CLGG",
    "FLYA",  "A02_FLYA",
    "YLQQ",  "A02_YLQQ",
    "LLDF",  "A02_LLDF",
    "ATIG",  "A11_ATIG",
    "AVFD*", "A11_AVFD*",
    "SSCS*", "A11_SSCS*",
    "RPPI",  "B07_RPPI",
    "RPQG",  "B07_RPQG",
    "RPQK*", "B07_RPQK*",
    "LPRR",  "B07_LPRR",
    "RAKF",  "B08_RAKF",
    "FLRG",  "B08_FLRG",
    "QAKW",  "B08_QAKW",
    "ELRS",  "B08_ELRS",
    "YNLR*", "B08_YNLR*",
    "RYGF",  "A24_RYGF",
    "LPFE",  "B35_LPFE",
    "QEIR",  "B40_QEIR",
    "LEKA",  "B40_LEKA"
  ),
  
  undetected_value = 0.00001,
  
  # Tetramers where only specific samples were stained; other sample/tetramer
  # combinations are removed as illegitimate sticky calls.
  allowed_tetramer_sample_map = tibble::tribble(
    ~tetramer,  ~sample,
    "B07_LPRR", "C246",
    "B07_LPRR", "C287",
    "B07_RPQG", "C246",
    "B07_RPQG", "C287",
    "B08_ELRS", "MS359",
    "B08_ELRS", "MS394",
    "B40_QEIR", "MS173",
    "B35_LPFE", "C339",
    "B35_LPFE", "MS069",
    "A24_RYGF", "MS045",
    "A24_RYGF", "MS136",
    "B40_LEKA", "MS358",
    "B40_LEKA", "MS173"
  ),
  
  tetramer_plot_order = c(
    "A02_GLCT", "A11_ATIG", "B08_RAKF",
    "B07_RPPI", "B08_FLRG", "B08_QAKW",
    "A11_AVFD*", "A02_LLDF", "A02_YLQQ",
    "A02_CLGG", "A02_FLYA", "A11_SSCS*",
    "B07_RPQK*", "B08_YNLR*"
  ),
  
  min_samples_per_cohort_for_plot = 2,
  tetramers_per_plot = 3,
  
  detected_colour = "#f07979",
  not_detected_colour = "#e5d9cd",
  
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
  create_parent_dir(cfg$baseline_tcr_export_csv)
  create_parent_dir(cfg$cd8_proportions_csv)
  create_parent_dir(cfg$responders_tetramer_freq_csv)
  create_parent_dir(cfg$baseline_tetramer_frequencies_csv)
  create_parent_dir(cfg$detected_summary_csv)
  create_parent_dir(cfg$session_info_file)
}


load_seurat_object <- function(path) {
  if (!file.exists(path)) {
    stop("Input RDS file does not exist: ", path)
  }
  
  obj <- readRDS(path)
  
  if (!inherits(obj, "Seurat")) {
    stop("Input RDS must contain a Seurat object.")
  }
  
  obj
}


get_assay_data_compat <- function(seurat_obj, assay, layer_or_slot = "data") {
  tryCatch(
    GetAssayData(seurat_obj, assay = assay, layer = layer_or_slot),
    error = function(e) {
      GetAssayData(seurat_obj, assay = assay, slot = layer_or_slot)
    }
  )
}


standardise_tetramer_name <- function(x, cfg) {
  x <- as.character(x)
  
  # If already HLA-prefixed, keep as-is.
  already_prefixed <- grepl("^[A-Za-z0-9]+_", x)
  
  lookup <- cfg$tetramer_prefix_lookup
  
  mapped <- lookup$tetramer_prefixed[
    match(x, lookup$tetramer)
  ]
  
  ifelse(
    already_prefixed,
    x,
    ifelse(!is.na(mapped), mapped, x)
  )
}


export_baseline_metadata <- function(obj, cfg) {
  export_df <- obj@meta.data |>
    select(any_of(cfg$cols_to_export)) |>
    rownames_to_column("cell_barcode")
  
  write.csv(
    export_df,
    cfg$baseline_tcr_export_csv,
    row.names = FALSE
  )
  
  export_df
}


add_adt_cd8_metadata <- function(obj, cfg) {
  if (!cfg$adt_assay %in% Assays(obj)) {
    stop("ADT assay not found in object: ", cfg$adt_assay)
  }
  
  DefaultAssay(obj) <- cfg$adt_assay
  
  adt_mat <- get_assay_data_compat(
    seurat_obj = obj,
    assay = cfg$adt_assay,
    layer_or_slot = cfg$adt_layer
  )
  
  required_features <- c(
    cfg$adt_cd8_feature,
    cfg$adt_cd3_feature,
    cfg$adt_cd4_feature
  )
  
  missing_features <- setdiff(required_features, rownames(adt_mat))
  
  if (length(missing_features) > 0) {
    stop(
      "Missing ADT feature(s): ",
      paste(missing_features, collapse = ", ")
    )
  }
  
  obj@meta.data <- obj@meta.data |>
    mutate(
      ADT_CD8 = as.numeric(adt_mat[cfg$adt_cd8_feature, colnames(obj)]),
      ADT_CD3 = as.numeric(adt_mat[cfg$adt_cd3_feature, colnames(obj)]),
      ADT_CD4_1 = as.numeric(adt_mat[cfg$adt_cd4_feature, colnames(obj)])
    )
  
  obj@meta.data <- obj@meta.data |>
    mutate(
      CD8 = ifelse(
        ADT_CD3 > cfg$adt_cd3_min &
          ADT_CD8 > cfg$adt_cd8_min &
          ADT_CD4_1 < cfg$adt_cd4_max,
        "YES",
        "NO"
      )
    )
  
  obj
}


estimate_cd8_proportions <- function(obj, cfg) {
  required_cols <- c("sample", "batch", "CD8")
  missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing metadata column(s) for CD8 proportion estimation: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  cd8_props <- obj@meta.data |>
    filter(
      sample %in% cfg$samples_keep,
      !grepl(cfg$global_pbmc_exclude_batch_pattern, batch)
    ) |>
    group_by(sample) |>
    summarise(
      total_cells = n(),
      cd8_yes = sum(CD8 == "YES", na.rm = TRUE),
      prop_cd8 = cd8_yes / total_cells,
      source = "global_PBMC_ADT",
      note = NA_character_,
      .groups = "drop"
    )
  
  if (nrow(cfg$manual_cd8_overrides) > 0) {
    manual_rows <- cfg$manual_cd8_overrides |>
      transmute(
        sample = sample,
        total_cells = NA_integer_,
        cd8_yes = NA_integer_,
        prop_cd8 = prop_cd8,
        source = "manual_override",
        note = note
      )
    
    cd8_props <- cd8_props |>
      bind_rows(manual_rows) |>
      distinct(sample, .keep_all = TRUE)
  }
  
  write.csv(
    cd8_props,
    cfg$cd8_proportions_csv,
    row.names = FALSE
  )
  
  cd8_props
}


make_tetramer_counts <- function(baseline_export_df, cfg) {
  required_cols <- c("id", "cohort", "sample", "batch", "tetramer", "antigen")
  missing_cols <- setdiff(required_cols, colnames(baseline_export_df))
  
  if (length(missing_cols) > 0) {
    stop(
      "Baseline metadata export is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  baseline_export_df |>
    filter(
      grepl(cfg$baseline_ebv_batch_pattern, batch),
      !is.na(id),
      !is.na(tetramer),
      !tetramer %in% cfg$tetramers_to_exclude_from_frequency
    ) |>
    mutate(
      tetramer = standardise_tetramer_name(tetramer, cfg)
    ) |>
    count(id, cohort, sample, tetramer, antigen, name = "n_cells") |>
    group_by(id) |>
    mutate(
      total_cells_id = sum(n_cells),
      within_enriched_freq = n_cells / total_cells_id
    ) |>
    ungroup()
}


load_frequency_input <- function(cfg, cd8_props) {
  if (!file.exists(cfg$frequency_input_xlsx)) {
    stop("Frequency input Excel file does not exist: ", cfg$frequency_input_xlsx)
  }
  
  input_df <- read_xlsx(cfg$frequency_input_xlsx) |>
    as.data.frame()
  
  required_cols <- c(
    "id",
    "sample",
    "Total PBMC count",
    "Cells sequenced",
    "Cells aligned"
  )
  
  missing_cols <- setdiff(required_cols, colnames(input_df))
  
  if (length(missing_cols) > 0) {
    stop(
      "Frequency input is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (!"proportion of CD8" %in% colnames(input_df)) {
    input_df <- input_df |>
      left_join(
        cd8_props |> select(sample, prop_cd8),
        by = "sample"
      ) |>
      mutate(`proportion of CD8` = prop_cd8) |>
      select(-prop_cd8)
  }
  
  input_df
}


calculate_corrected_tetramer_frequencies <- function(tetramer_counts, frequency_input, cfg) {
  join_cols <- intersect(c("id", "sample"), intersect(colnames(frequency_input), colnames(tetramer_counts)))
  
  if (length(join_cols) == 0) {
    stop("No shared join columns found between frequency input and tetramer counts.")
  }
  
  df_joined <- frequency_input |>
    left_join(
      tetramer_counts |>
        select(id, cohort, sample, tetramer, antigen, n_cells),
      by = join_cols,
      suffix = c("", ".tetramer")
    )
  
  # If cohort was not present in input, use tetramer-count cohort.
  if (!"cohort" %in% colnames(frequency_input) && "cohort" %in% colnames(df_joined)) {
    df_joined$cohort <- df_joined$cohort
  }
  
  df_joined <- df_joined |>
    mutate(
      estimated_tetramer_percentage = if_else(
        !is.na(n_cells) &
          `Cells aligned` > 0 &
          `Cells sequenced` > 0 &
          `proportion of CD8` > 0 &
          `Total PBMC count` > 0,
        (n_cells / `Cells aligned`) *
          (`Cells sequenced` / (`proportion of CD8` * `Total PBMC count`)) *
          100,
        NA_real_
      )
    )
  
  write.csv(
    df_joined,
    cfg$responders_tetramer_freq_csv,
    row.names = FALSE
  )
  
  df_joined
}


format_hla_long <- function(hla_raw) {
  hla_raw <- as.data.frame(hla_raw)
  
  # Case 1: already sample-by-allele matrix with a sample column.
  if ("sample" %in% colnames(hla_raw)) {
    return(
      hla_raw |>
        pivot_longer(
          cols = -sample,
          names_to = "allele",
          values_to = "sample_hla_status"
        ) |>
        filter(sample_hla_status == "positive") |>
        transmute(
          allele = as.character(allele),
          sample = as.character(sample)
        ) |>
        distinct()
    )
  }
  
  # Case 2: original no-header style, first column = allele, remaining values = samples.
  first_col <- colnames(hla_raw)[1]
  
  hla_raw |>
    rename(allele = all_of(first_col)) |>
    pivot_longer(
      cols = -allele,
      values_to = "sample"
    ) |>
    filter(!is.na(sample), sample != "") |>
    transmute(
      allele = as.character(allele),
      sample = as.character(sample)
    ) |>
    distinct()
}


complete_hla_expected_frequencies <- function(freq_df, cfg) {
  if (!file.exists(cfg$sample_hla_xlsx)) {
    stop("HLA input Excel file does not exist: ", cfg$sample_hla_xlsx)
  }
  
  hla_raw <- read_xlsx(
    cfg$sample_hla_xlsx,
    col_names = FALSE
  )
  
  hla_long <- format_hla_long(hla_raw)
  
  tet_ref <- freq_df |>
    filter(!is.na(tetramer), tetramer != "") |>
    select(tetramer, antigen) |>
    distinct() |>
    mutate(
      allele = sub("_.*$", "", tetramer)
    )
  
  expected <- hla_long |>
    inner_join(tet_ref, by = "allele") |>
    select(sample, tetramer, antigen) |>
    distinct()
  
  sample_info <- freq_df |>
    select(any_of(c("sample", "cohort", "id"))) |>
    distinct()
  
  out <- expected |>
    left_join(sample_info, by = "sample") |>
    left_join(
      freq_df |>
        select(
          any_of(c(
            "sample",
            "id",
            "cohort",
            "tetramer",
            "antigen",
            "estimated_tetramer_percentage"
          ))
        ),
      by = c("sample", "tetramer", "antigen"),
      suffix = c("", ".freq")
    ) |>
    mutate(
      cohort = coalesce(cohort, cohort.freq),
      estimated_tetramer_percentage = ifelse(
        is.na(estimated_tetramer_percentage),
        cfg$undetected_value,
        estimated_tetramer_percentage
      )
    ) |>
    select(
      estimated_tetramer_percentage,
      tetramer,
      antigen,
      sample,
      cohort,
      any_of("id")
    ) |>
    arrange(sample, tetramer)
  
  out
}


remove_illegitimate_tetramer_rows <- function(freq_complete, cfg) {
  allowed_map <- cfg$allowed_tetramer_sample_map
  tetramers_to_check <- unique(allowed_map$tetramer)
  
  illegit_rows <- freq_complete |>
    filter(tetramer %in% tetramers_to_check) |>
    anti_join(
      allowed_map,
      by = c("tetramer", "sample")
    )
  
  out_clean <- freq_complete |>
    anti_join(
      illegit_rows |>
        select(sample, tetramer) |>
        distinct(),
      by = c("sample", "tetramer")
    )
  
  out_clean
}


make_detection_summary <- function(freq_df, cfg) {
  df2 <- freq_df |>
    mutate(
      detected = if_else(
        estimated_tetramer_percentage > cfg$undetected_value,
        "Detected",
        "Not detected"
      ),
      cohort = as.character(cohort),
      tetramer = as.character(tetramer),
      sample = as.character(sample)
    ) |>
    filter(cohort %in% c("MS", "Control"))
  
  tetramer_keep <- df2 |>
    distinct(tetramer, cohort, sample) |>
    count(tetramer, cohort, name = "n_samples") |>
    filter(n_samples >= cfg$min_samples_per_cohort_for_plot) |>
    count(tetramer, name = "n_cohorts") |>
    filter(n_cohorts == 2) |>
    pull(tetramer)
  
  prop_df <- df2 |>
    filter(tetramer %in% tetramer_keep) |>
    group_by(tetramer, cohort) |>
    summarise(
      n_total = n(),
      n_detected = sum(detected == "Detected"),
      n_not_detected = sum(detected == "Not detected"),
      prop_detected = n_detected / n_total,
      prop_not_detected = n_not_detected / n_total,
      .groups = "drop"
    )
  
  write.csv(
    prop_df,
    cfg$detected_summary_csv,
    row.names = FALSE
  )
  
  list(
    df_detected = df2,
    prop_df = prop_df,
    tetramer_keep = tetramer_keep
  )
}


plot_detection_proportions <- function(detection_result, cfg) {
  df2 <- detection_result$df_detected
  
  tet_order <- cfg$tetramer_plot_order
  
  df_sub <- df2 |>
    filter(tetramer %in% tet_order)
  
  tet_keep <- df_sub |>
    distinct(tetramer, cohort, sample) |>
    count(tetramer, cohort, name = "n_samples") |>
    pivot_wider(
      names_from = cohort,
      values_from = n_samples,
      values_fill = 0
    ) |>
    filter(
      MS >= cfg$min_samples_per_cohort_for_plot,
      Control >= cfg$min_samples_per_cohort_for_plot
    ) |>
    pull(tetramer)
  
  tet_keep_ordered <- tet_order[tet_order %in% tet_keep]
  
  if (length(tet_keep_ordered) == 0) {
    warning("No tetramers passed plotting sample-count filter.")
    return(invisible(NULL))
  }
  
  plot_df <- df_sub |>
    filter(tetramer %in% tet_keep_ordered) |>
    mutate(
      cohort = factor(cohort, levels = c("Control", "MS")),
      detected = factor(detected, levels = c("Detected", "Not detected")),
      tetramer = factor(tetramer, levels = tet_keep_ordered)
    ) |>
    count(tetramer, cohort, detected, name = "n") |>
    group_by(tetramer, cohort) |>
    mutate(prop = n / sum(n)) |>
    ungroup()
  
  tetramer_groups <- split(
    tet_keep_ordered,
    ceiling(seq_along(tet_keep_ordered) / cfg$tetramers_per_plot)
  )
  
  make_plot <- function(tets, data) {
    d <- data |>
      filter(as.character(tetramer) %in% tets) |>
      mutate(tetramer = factor(as.character(tetramer), levels = tets)) |>
      droplevels()
    
    ggplot(d, aes(x = cohort, y = prop, fill = detected)) +
      geom_col(width = 0.75, color = "white", linewidth = 0.3) +
      facet_wrap(~ tetramer, ncol = cfg$tetramers_per_plot, drop = TRUE) +
      scale_y_continuous(
        labels = percent_format(accuracy = 1),
        expand = c(0, 0)
      ) +
      scale_fill_manual(
        values = c(
          "Detected" = cfg$detected_colour,
          "Not detected" = cfg$not_detected_colour
        )
      ) +
      labs(
        x = NULL,
        y = "Proportion of samples",
        fill = NULL
      ) +
      theme_classic(base_size = 16) +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 17, face = "bold"),
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        legend.position = "top",
        panel.spacing = unit(1.3, "lines")
      )
  }
  
  plots <- lapply(tetramer_groups, make_plot, data = plot_df)
  
  for (i in seq_along(plots)) {
    ggsave(
      filename = file.path(
        cfg$output_dir,
        paste0("Tetramer_proportions_set_", i, ".png")
      ),
      plot = plots[[i]],
      width = 9,
      height = 4,
      units = "in",
      dpi = 300,
      bg = "white"
    )
  }
  
  invisible(plots)
}


# ----------------------------- #
# 3. Run pipeline
# ----------------------------- #

setup_dirs(config)

message("Loading baseline module-scored object...")
baseline_obj <- load_seurat_object(config$baseline_module_scored_rds)

message("Exporting baseline TCR/module metadata...")
baseline_export <- export_baseline_metadata(
  obj = baseline_obj,
  cfg = config
)

message("Loading final merged all-object for CD8 proportion estimation...")
merged_obj <- load_seurat_object(config$merged_all_objects_rds)

message("Adding ADT-based CD8 calls...")
merged_obj <- add_adt_cd8_metadata(
  obj = merged_obj,
  cfg = config
)

message("Estimating CD8 proportions from global PBMC batches...")
cd8_props <- estimate_cd8_proportions(
  obj = merged_obj,
  cfg = config
)

message("Counting baseline EBV tetramer-positive cells...")
tetramer_counts <- make_tetramer_counts(
  baseline_export_df = baseline_export,
  cfg = config
)

message("Loading sample-level frequency input...")
frequency_input <- load_frequency_input(
  cfg = config,
  cd8_props = cd8_props
)

message("Calculating corrected tetramer frequencies...")
responders_freq <- calculate_corrected_tetramer_frequencies(
  tetramer_counts = tetramer_counts,
  frequency_input = frequency_input,
  cfg = config
)

message("Completing expected HLA-matched tetramer/sample combinations...")
freq_complete <- complete_hla_expected_frequencies(
  freq_df = responders_freq,
  cfg = config
)

message("Removing illegitimate sticky tetramer/sample combinations...")
freq_clean <- remove_illegitimate_tetramer_rows(
  freq_complete = freq_complete,
  cfg = config
)

write.csv(
  freq_clean,
  config$baseline_tetramer_frequencies_csv,
  row.names = FALSE
)

message("Calculating detected versus non-detected summaries...")
detection_result <- make_detection_summary(
  freq_df = freq_clean,
  cfg = config
)

message("Saving detected/non-detected stacked bar plots...")
plot_detection_proportions(
  detection_result = detection_result,
  cfg = config
)

writeLines(
  capture.output(sessionInfo()),
  config$session_info_file
)

message("\nBaseline EBV tetramer-frequency workflow complete.")
message("Saved baseline TCR export to: ", config$baseline_tcr_export_csv)
message("Saved CD8 proportions to: ", config$cd8_proportions_csv)
message("Saved responder frequency table to: ", config$responders_tetramer_freq_csv)
message("Saved HLA-completed frequency table to: ", config$baseline_tetramer_frequencies_csv)
message("Saved detected/non-detected summary to: ", config$detected_summary_csv)
message("Saved session info to: ", config$session_info_file)