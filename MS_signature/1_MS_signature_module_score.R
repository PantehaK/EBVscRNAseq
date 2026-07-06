#!/usr/bin/env Rscript

# ==============================================================================
# Donor-level MS transcriptional signature score in stimulated latent EBV T cells
# ==============================================================================
#
# Purpose:
#   This script scores an MS-associated gene signature in an activated/baseline
#   integrated Seurat object, then compares donor-level mean signature scores
#   between MS and Control participants in stimulated latent EBV-specific T cells.
#
#   Workflow:
#     1. Load activated/baseline integrated Seurat RDS.
#     2. Define latent_group from batch and lifecycle:
#          GEMEBV + Latent  = Baseline Latent
#          EBVLCL + Latent  = Stimulated Latent
#          all other cells  = Non-Latent
#     3. Add MS signature module score using Seurat::AddModuleScore().
#     4. Keep stimulated latent cells.
#     5. Calculate donor-level mean MS signature score.
#     6. Compare MS vs Control using a two-sided Wilcoxon rank-sum test.
#     7. Calculate rank-biserial effect size.
#     8. Export donor-level source data, statistics table, gene list, and session info.
#
# Notes:
#   - Paths are intentionally generic for GitHub.
#   - Update the config block before running.
#   - The donor/sample column defaults to `sample`; change to `id` if needed.
#   - The module score column created by AddModuleScore is usually
#     paste0(config$score_name, "1"), e.g. MS_signature1.
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readr)
})

config <- list(
  input_rds = "path/to/Publication_RDS/Activated_EBV_integrated_seurat.rds",
  output_dir = "path/to/Publication_data/EBV_activated/MS_signature_score",
  donor_col = "sample",
  cohort_col = "cohort",
  batch_col = "batch",
  lifecycle_col = "lifecycle",
  baseline_batch_pattern = "^GEMEBV",
  stimulated_batch_pattern = "^EBVLCL",
  latent_lifecycle_label = "Latent",
  baseline_label = "Baseline Latent",
  stimulated_label = "Stimulated Latent",
  non_latent_label = "Non-Latent",
  control_label = "Control",
  ms_label = "MS",
  assay = "RNA",
  score_name = "MS_signature",
  gene_set = c(
    "FTH1", "JAK1", "ITM2B", "RPL13", "SMCHD1", "IL7R", "HLA-B", "PTMA",
    "HERC5", "TOMM7", "MT-CYB", "H3-3A", "HLA-F", "RPS10", "MT-ND4L",
    "HNRNPA1", "FYN", "PFN1", "RPS29", "MT-CO3", "RPS15", "RPS18",
    "SLFN5", "ZBTB20", "RPLP2", "MT-ND3", "ARPC2", "IFI6", "CD3D",
    "RPS27", "HLA-C", "CXCR4", "RPS27A", "CALM1", "CCSER2", "ANKRD12",
    "VIM", "RPS9", "CDC42SE2", "EEF1G", "RPS6", "PDE3B", "EEF1A1",
    "RPS16", "YWHAB", "RPL30", "SARAF", "SRSF7"
  ),
  signature_gene_list_csv = "MS_signature_gene_list.csv",
  signature_gene_missing_csv = "MS_signature_genes_missing_from_object.csv",
  cell_level_metadata_csv = "MS_signature_cell_level_metadata.csv",
  donor_mean_signature_csv = "Mean_MS_signature.csv",
  stats_table_csv = "Mean_MS_signature_stats.csv",
  session_info_file = "sessionInfo_MS_signature_stimulated_latent.txt"
)

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
  x <- str_trim(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL", "null", "None", "none")] <- NA_character_
  x
}

clean_cohort <- function(x, control_label = "Control", ms_label = "MS") {
  x_raw <- as.character(x)
  x_clean <- str_to_lower(str_trim(x_raw))
  case_when(
    x_clean %in% c(
      "control", "controls", "healthy control", "healthy controls",
      "ctrl", "hc", "nms", "non-ms", "nonms", "non_ms"
    ) ~ control_label,
    x_clean %in% c("ms", "rrms", "pwms", "multiple sclerosis") ~ ms_label,
    TRUE ~ x_raw
  )
}

safe_quantile <- function(x, prob) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  unname(quantile(x, prob, na.rm = TRUE))
}

safe_min <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

calculate_rank_biserial <- function(ms_vals, control_vals) {
  ms_vals <- ms_vals[!is.na(ms_vals)]
  control_vals <- control_vals[!is.na(control_vals)]
  n_ms <- length(ms_vals)
  n_control <- length(control_vals)
  if (n_ms == 0 || n_control == 0) {
    return(tibble(
      n_MS = n_ms,
      n_Control = n_control,
      U_MS_vs_Control = NA_real_,
      rank_biserial = NA_real_
    ))
  }
  U_ms <- sum(outer(ms_vals, control_vals, ">")) +
    0.5 * sum(outer(ms_vals, control_vals, "=="))
  rank_biserial <- (2 * U_ms / (n_ms * n_control)) - 1
  tibble(
    n_MS = n_ms,
    n_Control = n_control,
    U_MS_vs_Control = U_ms,
    rank_biserial = rank_biserial
  )
}

format_range <- function(x1, x2, digits = 4) {
  if (is.na(x1) || is.na(x2)) return(NA_character_)
  paste0(round(x1, digits), " - ", round(x2, digits))
}

create_dir(config$output_dir)

message("Loading activated EBV object:")
message(config$input_rds)

merged <- load_seurat_object(config$input_rds)

required_cols <- c(
  config$donor_col,
  config$cohort_col,
  config$batch_col,
  config$lifecycle_col
)

check_required_columns(
  df = merged@meta.data,
  required_cols = required_cols,
  object_name = "activated EBV metadata"
)

if (!config$assay %in% Assays(merged)) {
  stop(
    "Assay not found: ",
    config$assay,
    ". Available assays: ",
    paste(Assays(merged), collapse = ", ")
  )
}

DefaultAssay(merged) <- config$assay

meta_clean <- merged@meta.data %>%
  as.data.frame() %>%
  mutate(
    donor_id = clean_chr(.data[[config$donor_col]]),
    cohort_clean = clean_cohort(
      .data[[config$cohort_col]],
      control_label = config$control_label,
      ms_label = config$ms_label
    ),
    batch_clean = clean_chr(.data[[config$batch_col]]),
    lifecycle_clean = clean_chr(.data[[config$lifecycle_col]]),
    latent_group = case_when(
      !is.na(batch_clean) &
        grepl(config$baseline_batch_pattern, batch_clean) &
        lifecycle_clean == config$latent_lifecycle_label ~ config$baseline_label,
      !is.na(batch_clean) &
        grepl(config$stimulated_batch_pattern, batch_clean) &
        lifecycle_clean == config$latent_lifecycle_label ~ config$stimulated_label,
      TRUE ~ config$non_latent_label
    ),
    latent_group = factor(
      latent_group,
      levels = c(
        config$non_latent_label,
        config$baseline_label,
        config$stimulated_label
      )
    )
  )

rownames(meta_clean) <- colnames(merged)
merged@meta.data <- meta_clean

message("Latent group counts:")
print(table(merged$latent_group, useNA = "ifany"))

message("Cohort counts:")
print(table(merged$cohort_clean, useNA = "ifany"))

available_genes <- rownames(merged[[config$assay]])
genes_present <- intersect(config$gene_set, available_genes)
genes_missing <- setdiff(config$gene_set, available_genes)

if (length(genes_present) == 0) {
  stop("None of the signature genes were found in assay: ", config$assay)
}

if (length(genes_missing) > 0) {
  warning(
    "Some signature genes were not found in the object and will be omitted: ",
    paste(genes_missing, collapse = ", ")
  )
}

signature_gene_list_path <- file.path(
  config$output_dir,
  config$signature_gene_list_csv
)

signature_gene_missing_path <- file.path(
  config$output_dir,
  config$signature_gene_missing_csv
)

write_csv(
  tibble(
    gene = config$gene_set,
    present_in_object = config$gene_set %in% available_genes
  ),
  signature_gene_list_path
)

write_csv(
  tibble(gene = genes_missing),
  signature_gene_missing_path
)

merged <- AddModuleScore(
  object = merged,
  features = list(genes_present),
  assay = config$assay,
  name = config$score_name
)

score_col <- paste0(config$score_name, "1")

if (!score_col %in% colnames(merged@meta.data)) {
  stop("Expected module score column was not created: ", score_col)
}

cell_level_metadata <- merged@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell_barcode") %>%
  transmute(
    cell_barcode,
    donor_id,
    cohort = cohort_clean,
    latent_group,
    batch = batch_clean,
    lifecycle = lifecycle_clean,
    MS_signature = .data[[score_col]]
  )

cell_level_metadata_path <- file.path(
  config$output_dir,
  config$cell_level_metadata_csv
)

write_csv(cell_level_metadata, cell_level_metadata_path)

stim_latent_cells <- cell_level_metadata %>%
  filter(
    latent_group == config$stimulated_label,
    cohort %in% c(config$control_label, config$ms_label),
    !is.na(donor_id),
    donor_id != "",
    !is.na(MS_signature)
  )

if (nrow(stim_latent_cells) == 0) {
  stop("No stimulated latent cells remained for donor-level MS signature analysis.")
}

mean_signature <- stim_latent_cells %>%
  group_by(donor_id, cohort, latent_group) %>%
  summarise(
    MS_signature = mean(MS_signature, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(cohort, donor_id)

message("Donor counts in stimulated latent analysis:")
print(mean_signature %>% count(cohort))

message("Checking duplicate donor rows:")
print(
  mean_signature %>%
    count(donor_id) %>%
    filter(n > 1)
)

donor_mean_signature_path <- file.path(
  config$output_dir,
  config$donor_mean_signature_csv
)

write_csv(mean_signature, donor_mean_signature_path)

cohort_counts <- mean_signature %>%
  count(cohort, name = "n_donors")

has_both_groups <- all(
  c(config$control_label, config$ms_label) %in% cohort_counts$cohort
)

if (!has_both_groups) {
  warning("Both Control and MS cohorts are required for Wilcoxon testing.")
  wilcox_p <- NA_real_
  wilcox_W <- NA_real_
} else {
  wilcox_res <- wilcox.test(
    MS_signature ~ cohort,
    data = mean_signature,
    exact = FALSE
  )
  wilcox_p <- wilcox_res$p.value
  wilcox_W <- unname(wilcox_res$statistic)
}

control_vals <- mean_signature %>%
  filter(cohort == config$control_label) %>%
  pull(MS_signature)

ms_vals <- mean_signature %>%
  filter(cohort == config$ms_label) %>%
  pull(MS_signature)

rb <- calculate_rank_biserial(
  ms_vals = ms_vals,
  control_vals = control_vals
)

stats_table <- mean_signature %>%
  group_by(cohort) %>%
  summarise(
    group = config$stimulated_label,
    n = sum(!is.na(MS_signature)),
    mean = mean(MS_signature, na.rm = TRUE),
    sd = sd(MS_signature, na.rm = TRUE),
    sem = sd / sqrt(n),
    median = median(MS_signature, na.rm = TRUE),
    q1 = safe_quantile(MS_signature, 0.25),
    q3 = safe_quantile(MS_signature, 0.75),
    min_value = safe_min(MS_signature),
    max_value = safe_max(MS_signature),
    IQR = format_range(q1, q3),
    min_max = format_range(min_value, max_value),
    .groups = "drop"
  ) %>%
  mutate(
    comparison = paste(config$ms_label, "vs", config$control_label),
    test = "Two-sided Mann-Whitney U test / Wilcoxon rank-sum test",
    wilcoxon_W_reported_by_R = wilcox_W,
    U_MS_vs_Control = rb$U_MS_vs_Control,
    p_value = wilcox_p,
    effect_size = "Rank-biserial correlation",
    rank_biserial = rb$rank_biserial,
    interpretation = case_when(
      is.na(rank_biserial) ~ NA_character_,
      rank_biserial > 0 ~ paste("Higher values in", config$ms_label),
      rank_biserial < 0 ~ paste("Higher values in", config$control_label),
      TRUE ~ "No directional difference"
    )
  ) %>%
  select(
    group,
    cohort,
    n,
    mean,
    sd,
    sem,
    median,
    IQR,
    min_max,
    comparison,
    test,
    wilcoxon_W_reported_by_R,
    U_MS_vs_Control,
    p_value,
    effect_size,
    rank_biserial,
    interpretation
  )

stats_table_path <- file.path(
  config$output_dir,
  config$stats_table_csv
)

write_csv(stats_table, stats_table_path)

message("Statistics table:")
print(stats_table)

session_info_path <- file.path(
  config$output_dir,
  config$session_info_file
)

writeLines(
  c(
    "Configuration:",
    capture.output(str(config)),
    "",
    "Score column:",
    score_col,
    "",
    "Number of genes in signature:",
    as.character(length(config$gene_set)),
    "",
    "Number of genes present in object:",
    as.character(length(genes_present)),
    "",
    "Number of genes missing from object:",
    as.character(length(genes_missing)),
    "",
    "Latent group counts:",
    capture.output(print(table(merged$latent_group, useNA = "ifany"))),
    "",
    "Donor counts for stimulated latent analysis:",
    capture.output(print(mean_signature %>% count(cohort))),
    "",
    "Statistics table:",
    capture.output(print(stats_table)),
    "",
    "Output files:",
    signature_gene_list_path,
    signature_gene_missing_path,
    cell_level_metadata_path,
    donor_mean_signature_path,
    stats_table_path,
    "",
    "Session information:",
    capture.output(sessionInfo())
  ),
  session_info_path
)

message("\nStimulated latent MS signature scoring workflow complete.")
message("Signature gene list: ", signature_gene_list_path)
message("Missing signature genes: ", signature_gene_missing_path)
message("Cell-level metadata: ", cell_level_metadata_path)
message("Donor mean signature: ", donor_mean_signature_path)
message("Statistics table: ", stats_table_path)
message("Session info: ", session_info_path)
