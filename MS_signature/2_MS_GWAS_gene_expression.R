#!/usr/bin/env Rscript

# ==============================================================================
# Sample-level logCPM export for selected MS-associated genes
# ==============================================================================
#
# Purpose:
#   This script extracts pseudobulk sample-level expression for selected genes
#   in stimulated latent EBV-specific T cells from an activated/baseline integrated
#   Seurat object.
#
#   It:
#     1. Loads the activated/baseline integrated Seurat RDS.
#     2. Defines latent_group from batch and lifecycle.
#     3. Subsets to stimulated latent cells from MS and Control participants.
#     4. Further subsets to a target cluster, default new_cluster == "2".
#     5. Aggregates RNA counts per sample.
#     6. Calculates logCPM using edgeR.
#     7. Exports sample-level logCPM values for selected genes.
#     8. Runs per-gene Wilcoxon rank-sum tests comparing MS vs Control.
#     9. Applies BH correction across the selected genes.
#    10. Saves outputs and session information.
#
# Notes:
#   - Paths are intentionally generic for GitHub.
#   - Update the config block before running.
#   - The donor/sample column defaults to `sample`; change to `id` if needed.
#   - Control is used as the reference ordering for display/export.
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
  library(edgeR)
  library(Matrix)
  library(readr)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  
  # Input RDS
  input_rds = "path/to/Publication_RDS/Activated_EBV_integrated_seurat.rds",
  
  # Output directory
  output_dir = "path/to/Publication_data/EBV_activated/DEG/gene_panel_logCPM",
  
  # Assay and counts layer/slot
  assay = "RNA",
  counts_layer = "counts",
  
  # Metadata columns
  donor_col = "sample",
  cohort_col = "cohort",
  cluster_col = "new_cluster",
  batch_col = "batch",
  lifecycle_col = "lifecycle",
  
  # Group definitions
  baseline_batch_pattern = "^GEMEBV",
  stimulated_batch_pattern = "^EBVLCL",
  latent_lifecycle_label = "Latent",
  baseline_label = "Baseline Latent",
  stimulated_label = "Stimulated Latent",
  non_latent_label = "Non-Latent",
  
  # Analysis subset
  target_latent_group = "Stimulated Latent",
  target_cluster = "2",
  cohort_order = c("Control", "MS"),
  
  # Genes to export
  genes_to_plot = c(
    "IL7R",
    "JAK1",
    "CXCR4"
  ),
  
  # Output files
  filtered_metadata_csv = "StimulatedLatent_cluster2_gene_panel_filtered_metadata.csv",
  genes_requested_csv = "StimulatedLatent_cluster2_gene_panel_requested_genes.csv",
  sample_logcpm_csv = "StimulatedLatent_cluster2_gene_panel_logCPM.csv",
  sample_logcpm_with_pvalues_csv = "StimulatedLatent_cluster2_gene_panel_logCPM_with_pvalues.csv",
  stats_csv = "StimulatedLatent_cluster2_gene_panel_Wilcoxon_BH_stats.csv",
  session_info_file = "sessionInfo_StimulatedLatent_cluster2_gene_panel_logCPM.txt"
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


get_assay_counts_compat <- function(object, assay = "RNA", counts_layer = "counts") {
  
  # Compatible with Seurat v5 and Seurat v4.
  tryCatch(
    {
      GetAssayData(object, assay = assay, layer = counts_layer)
    },
    error = function(e) {
      GetAssayData(object, assay = assay, slot = counts_layer)
    }
  )
}


clean_chr <- function(x) {
  
  x <- str_trim(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL", "null", "None", "none")] <- NA_character_
  x
}


clean_cohort <- function(x, cohort_order = c("Control", "MS")) {
  
  x_raw <- as.character(x)
  x_clean <- str_to_lower(str_trim(x_raw))
  
  out <- case_when(
    x_clean %in% c(
      "control", "controls", "healthy control", "healthy controls",
      "ctrl", "hc", "nms", "non-ms", "nonms", "non_ms"
    ) ~ "Control",
    
    x_clean %in% c(
      "ms", "rrms", "pwms", "multiple sclerosis"
    ) ~ "MS",
    
    TRUE ~ x_raw
  )
  
  factor(out, levels = cohort_order)
}


safe_wilcox <- function(df, value_col = "logCPM", group_col = "cohort") {
  
  df_test <- df %>%
    filter(
      !is.na(.data[[value_col]]),
      !is.na(.data[[group_col]])
    )
  
  if (n_distinct(df_test[[group_col]]) < 2) {
    return(NA_real_)
  }
  
  if (length(unique(df_test[[value_col]])) < 2) {
    return(1)
  }
  
  tryCatch(
    {
      wilcox.test(
        df_test[[value_col]] ~ df_test[[group_col]],
        exact = FALSE
      )$p.value
    },
    error = function(e) NA_real_
  )
}


format_p <- function(p) {
  
  case_when(
    is.na(p) ~ NA_character_,
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
    TRUE ~ "ns"
  )
}


# ----------------------------- #
# 3. Load object and prepare metadata
# ----------------------------- #

create_dir(config$output_dir)

message("Loading activated EBV object:")
message(config$input_rds)

merged <- load_seurat_object(config$input_rds)

required_cols <- c(
  config$donor_col,
  config$cohort_col,
  config$cluster_col,
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

counts_all <- get_assay_counts_compat(
  object = merged,
  assay = config$assay,
  counts_layer = config$counts_layer
)

meta_all <- merged@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell_barcode") %>%
  mutate(
    donor_id = clean_chr(.data[[config$donor_col]]),
    cohort = clean_cohort(
      .data[[config$cohort_col]],
      cohort_order = config$cohort_order
    ),
    cluster = clean_chr(.data[[config$cluster_col]]),
    batch = clean_chr(.data[[config$batch_col]]),
    lifecycle = clean_chr(.data[[config$lifecycle_col]]),
    
    latent_group = case_when(
      !is.na(batch) &
        grepl(config$baseline_batch_pattern, batch) &
        lifecycle == config$latent_lifecycle_label ~ config$baseline_label,
      
      !is.na(batch) &
        grepl(config$stimulated_batch_pattern, batch) &
        lifecycle == config$latent_lifecycle_label ~ config$stimulated_label,
      
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

if (!identical(colnames(counts_all), meta_all$cell_barcode)) {
  stop("Counts columns and metadata cell barcodes are not aligned.")
}

message("Latent group counts:")
print(table(meta_all$latent_group, useNA = "ifany"))

message("Cluster counts:")
print(table(meta_all$cluster, useNA = "ifany"))


# ----------------------------- #
# 4. Filter to target group and cluster
# ----------------------------- #

meta_use <- meta_all %>%
  filter(
    latent_group == config$target_latent_group,
    cohort %in% config$cohort_order,
    cluster == config$target_cluster,
    !is.na(donor_id), donor_id != ""
  )

if (nrow(meta_use) == 0) {
  stop("No cells remained after stimulated latent / cohort / cluster filtering.")
}

counts_use <- counts_all[, meta_use$cell_barcode, drop = FALSE]

filtered_metadata_path <- file.path(
  config$output_dir,
  config$filtered_metadata_csv
)

write_csv(meta_use, filtered_metadata_path)

message("Cells retained for gene-panel export:")
print(table(meta_use$cohort, useNA = "ifany"))

message("Donors retained:")
print(
  meta_use %>%
    distinct(donor_id, cohort) %>%
    count(cohort, name = "n_donors")
)


# ----------------------------- #
# 5. Pseudobulk counts per sample
# ----------------------------- #

donor_factor <- factor(meta_use$donor_id)

aggregation_matrix <- Matrix::sparseMatrix(
  i = seq_along(donor_factor),
  j = as.integer(donor_factor),
  x = 1,
  dims = c(length(donor_factor), nlevels(donor_factor)),
  dimnames = list(meta_use$cell_barcode, levels(donor_factor))
)

pb_counts <- counts_use[, meta_use$cell_barcode, drop = FALSE] %*% aggregation_matrix
pb_counts <- as.matrix(pb_counts)

sample_meta <- meta_use %>%
  distinct(donor_id, cohort) %>%
  filter(donor_id %in% colnames(pb_counts)) %>%
  arrange(match(donor_id, colnames(pb_counts))) %>%
  mutate(
    cohort = factor(cohort, levels = config$cohort_order)
  )

if (!identical(sample_meta$donor_id, colnames(pb_counts))) {
  stop("Pseudobulk count columns and sample metadata are not aligned.")
}


# ----------------------------- #
# 6. CPM + log transform
# ----------------------------- #

dge <- edgeR::DGEList(counts = pb_counts)
dge <- edgeR::calcNormFactors(dge)

logCPM <- edgeR::cpm(
  dge,
  log = TRUE,
  prior.count = 1
)

genes_present <- intersect(config$genes_to_plot, rownames(logCPM))
genes_missing <- setdiff(config$genes_to_plot, rownames(logCPM))

genes_requested_path <- file.path(
  config$output_dir,
  config$genes_requested_csv
)

write_csv(
  tibble(
    gene = config$genes_to_plot,
    present_in_logCPM = config$genes_to_plot %in% rownames(logCPM)
  ),
  genes_requested_path
)

if (length(genes_present) == 0) {
  stop("None of config$genes_to_plot were found in the logCPM matrix.")
}

if (length(genes_missing) > 0) {
  warning(
    "These requested genes were not found and will be omitted: ",
    paste(genes_missing, collapse = ", ")
  )
}


# ----------------------------- #
# 7. Build sample-level long dataframe
# ----------------------------- #

df_logcpm <- as.data.frame(t(logCPM[genes_present, , drop = FALSE])) %>%
  rownames_to_column("donor_id") %>%
  pivot_longer(
    cols = -donor_id,
    names_to = "gene",
    values_to = "logCPM"
  ) %>%
  left_join(sample_meta, by = "donor_id") %>%
  filter(!is.na(cohort)) %>%
  mutate(
    cohort = factor(cohort, levels = config$cohort_order),
    latent_group = config$target_latent_group,
    cluster = config$target_cluster
  ) %>%
  select(
    latent_group,
    cluster,
    donor_id,
    cohort,
    gene,
    logCPM
  ) %>%
  arrange(gene, cohort, donor_id)

sample_logcpm_path <- file.path(
  config$output_dir,
  config$sample_logcpm_csv
)

write_csv(df_logcpm, sample_logcpm_path)


# ----------------------------- #
# 8. Per-gene Wilcoxon tests
# ----------------------------- #

stats_df <- df_logcpm %>%
  group_by(gene) %>%
  summarise(
    n_Control = sum(cohort == "Control" & !is.na(logCPM)),
    n_MS = sum(cohort == "MS" & !is.na(logCPM)),
    median_Control = median(logCPM[cohort == "Control"], na.rm = TRUE),
    median_MS = median(logCPM[cohort == "MS"], na.rm = TRUE),
    diff_median_MS_minus_Control = median_MS - median_Control,
    p_value = safe_wilcox(cur_data(), value_col = "logCPM", group_col = "cohort"),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    
    # Raw P-value stars, matching the original script.
    signif_raw_p = p_to_stars(p_value),
    
    # FDR stars are included separately for clarity.
    signif_FDR = p_to_stars(p_adj),
    
    p_value_label = format_p(p_value),
    p_adj_label = format_p(p_adj)
  ) %>%
  arrange(p_adj, p_value)

stats_path <- file.path(
  config$output_dir,
  config$stats_csv
)

write_csv(stats_df, stats_path)


# ----------------------------- #
# 9. Join p-values back onto sample-level data
# ----------------------------- #

df_export <- df_logcpm %>%
  left_join(stats_df, by = "gene") %>%
  arrange(gene, cohort, donor_id)

sample_logcpm_with_pvalues_path <- file.path(
  config$output_dir,
  config$sample_logcpm_with_pvalues_csv
)

write_csv(df_export, sample_logcpm_with_pvalues_path)


# ----------------------------- #
# 10. Save session information
# ----------------------------- #

session_info_path <- file.path(
  config$output_dir,
  config$session_info_file
)

writeLines(
  c(
    "Configuration:",
    capture.output(str(config)),
    "",
    "Latent group counts:",
    capture.output(print(table(meta_all$latent_group, useNA = "ifany"))),
    "",
    "Cells retained for analysis:",
    as.character(nrow(meta_use)),
    "",
    "Cells by cohort:",
    capture.output(print(table(meta_use$cohort, useNA = "ifany"))),
    "",
    "Donors retained:",
    capture.output(print(meta_use %>% distinct(donor_id, cohort) %>% count(cohort, name = "n_donors"))),
    "",
    "Genes present:",
    paste(genes_present, collapse = ", "),
    "",
    "Genes missing:",
    paste(genes_missing, collapse = ", "),
    "",
    "Statistics:",
    capture.output(print(stats_df)),
    "",
    "Output files:",
    filtered_metadata_path,
    genes_requested_path,
    sample_logcpm_path,
    sample_logcpm_with_pvalues_path,
    stats_path,
    "",
    "Session information:",
    capture.output(sessionInfo())
  ),
  session_info_path
)


# ----------------------------- #
# 11. Completion messages
# ----------------------------- #

message("\nStimulated latent cluster gene-panel logCPM workflow complete.")
message("Filtered metadata: ", filtered_metadata_path)
message("Requested genes: ", genes_requested_path)
message("Sample logCPM: ", sample_logcpm_path)
message("Sample logCPM with p-values: ", sample_logcpm_with_pvalues_path)
message("Stats: ", stats_path)
message("Session info: ", session_info_path)

message("\nStats:")
print(stats_df)
