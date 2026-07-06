#!/usr/bin/env Rscript

# ==============================================================================
# edgeR pseudobulk: Latent vs Lytic EBV-specific T cells by cell type
# ==============================================================================
#
# Purpose:
#   This script performs pseudobulk differential expression comparing latent-
#   directed versus lytic-directed EBV-specific T cells within each T cell
#   cluster/cell type.
#
#   Starting object:
#     5_baseline_EBV_cluster_reannotated.rds
#
#   Comparison:
#     Latent vs Lytic, within each celltype_new group
#
#   Pseudobulk unit:
#     donor/sample x celltype_new x latency_group
#
#   Model:
#     ~ donor/sample + latency_group
#
#   Because Lytic is the reference:
#     logFC > 0 = higher in Latent-directed cells
#     logFC < 0 = higher in Lytic-directed cells
#
# Notes:
#   - Paths are intentionally generic for GitHub.
#   - Update the config block before running.
#   - The donor/sample column defaults to `sample`, matching the original script.
#     Change config$donor_col to "id" if your analysis should use the `id`
#     column instead.
#   - This script uses raw RNA counts for edgeR pseudobulk.
#
# ==============================================================================


# ----------------------------- #
# 0. Load packages
# ----------------------------- #

suppressPackageStartupMessages({
  library(Seurat)
  library(edgeR)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(stringr)
  library(readr)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  
  # Input RDS
  input_rds = "path/to/Publication_RDS/Baseline_EBV_integrated_seurat.rds",
  
  # Output directory
  output_dir = "path/to/Publication_data/EBV_baseline/DEG/edgeR_pseudobulk_Latent_vs_Lytic_by_celltype_new",
  
  # Assay and count layer/slot
  assay = "RNA",
  counts_layer = "counts",
  
  # Metadata columns
  latency_col = "latency",
  celltype_col = "celltype_new",
  donor_col = "sample",
  
  # Latency grouping
  latent_levels = c("Latency I/0", "Latency II", "Latency III"),
  lytic_level = "Lytic",
  latency_reference = "Lytic",
  latency_test_level = "Latent",
  
  # Pseudobulk filters
  min_cells_per_pseudobulk = 10,
  min_paired_donors_per_celltype = 3,
  
  # edgeR significance threshold
  fdr_threshold = 0.05,
  abs_logfc_threshold = 2,
  
  # Optional fixed cell type order.
  # Leave as NULL to use alphabetical order from the data.
  celltype_order = c(
    "Naive/early TCM",
    "TEM",
    "Late TEM",
    "CTL",
    "CD69+ early activated T",
    "innate-like T",
    "CD69+ TEM"
  ),
  
  # Output files
  pseudobulk_metadata_csv = "pseudobulk_metadata_sample_celltype_latency.csv",
  pseudobulk_metadata_filtered_csv = "pseudobulk_metadata_after_min_cell_filter.csv",
  pseudobulk_summary_csv = "pseudobulk_counts_by_celltype_new.csv",
  combined_results_csv = "edgeR_pseudobulk_Latent_vs_Lytic_by_celltype_new_ALL.csv",
  significant_results_csv = "edgeR_pseudobulk_SIGNIFICANT_Latent_vs_Lytic_by_celltype_new.csv",
  significant_summary_csv = "edgeR_pseudobulk_significant_gene_summary_by_celltype_new.csv",
  skipped_celltypes_csv = "edgeR_pseudobulk_skipped_celltypes.csv",
  session_info_file = "sessionInfo_edgeR_pseudobulk_Latent_vs_Lytic_by_celltype_new.txt"
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


safe_name <- function(x) {
  
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+$", "") %>%
    str_replace_all("^_+", "")
}


assign_latency_group <- function(x, latent_levels, lytic_level) {
  
  case_when(
    x %in% latent_levels ~ "Latent",
    x == lytic_level ~ "Lytic",
    TRUE ~ NA_character_
  )
}


format_direction <- function(logfc, fdr, fdr_threshold, abs_logfc_threshold) {
  
  case_when(
    !is.na(fdr) & fdr < fdr_threshold & !is.na(logfc) & logfc > abs_logfc_threshold ~ "Latent_up",
    !is.na(fdr) & fdr < fdr_threshold & !is.na(logfc) & logfc < -abs_logfc_threshold ~ "Lytic_up",
    TRUE ~ "NS"
  )
}


# ----------------------------- #
# 3. Load object and check metadata
# ----------------------------- #

create_dir(config$output_dir)

message("Loading baseline EBV object:")
message(config$input_rds)

merged <- load_seurat_object(config$input_rds)

required_cols <- c(
  config$latency_col,
  config$celltype_col,
  config$donor_col
)

check_required_columns(
  df = merged@meta.data,
  required_cols = required_cols,
  object_name = "baseline EBV metadata"
)

message("Latency values before grouping:")
print(table(merged@meta.data[[config$latency_col]], useNA = "ifany"))


# ----------------------------- #
# 4. Collapse latency into Latent vs Lytic
# ----------------------------- #

merged$latency_group <- assign_latency_group(
  x = merged@meta.data[[config$latency_col]],
  latent_levels = config$latent_levels,
  lytic_level = config$lytic_level
)

merged$latency_group <- factor(
  merged$latency_group,
  levels = c(config$latency_reference, config$latency_test_level)
)

message("Latency group counts:")
print(table(merged$latency_group, useNA = "ifany"))


# ----------------------------- #
# 5. Subset to usable cells
# ----------------------------- #

obj <- subset(
  merged,
  subset = !is.na(latency_group) &
    !is.na(.data[[config$celltype_col]]) &
    !is.na(.data[[config$donor_col]])
)

if (ncol(obj) == 0) {
  stop("No cells remained after filtering for latency_group, celltype and donor/sample.")
}

if (!config$assay %in% Assays(obj)) {
  stop(
    "Assay not found: ",
    config$assay,
    ". Available assays: ",
    paste(Assays(obj), collapse = ", ")
  )
}

DefaultAssay(obj) <- config$assay

counts <- get_assay_counts_compat(
  object = obj,
  assay = config$assay,
  counts_layer = config$counts_layer
)

meta <- obj@meta.data

if (!identical(colnames(counts), rownames(meta))) {
  stop("Counts columns and metadata rownames are not aligned.")
}

meta$donor_id <- factor(as.character(meta[[config$donor_col]]))
meta$celltype_new_for_pb <- as.character(meta[[config$celltype_col]])
meta$latency_group <- factor(
  as.character(meta$latency_group),
  levels = c(config$latency_reference, config$latency_test_level)
)

message("Cells retained for pseudobulk:")
print(table(meta$latency_group, useNA = "ifany"))

message("Cells per celltype and latency group:")
print(table(meta$celltype_new_for_pb, meta$latency_group, useNA = "ifany"))


# ----------------------------- #
# 6. Create pseudobulk IDs
# ----------------------------- #
#
# One pseudobulk = donor/sample x celltype_new x latency_group
#
# ----------------------------- #

meta$pb_id <- paste(
  meta$donor_id,
  meta$celltype_new_for_pb,
  meta$latency_group,
  sep = "___"
)

pb_factor <- factor(meta$pb_id)

aggregation_matrix <- Matrix::sparseMatrix(
  i = seq_along(pb_factor),
  j = as.integer(pb_factor),
  x = 1,
  dims = c(length(pb_factor), nlevels(pb_factor)),
  dimnames = list(rownames(meta), levels(pb_factor))
)

pb_counts <- counts[, rownames(meta), drop = FALSE] %*% aggregation_matrix

# edgeR is happiest with a regular matrix after aggregation.
pb_counts <- as.matrix(pb_counts)

message("Pseudobulk count matrix dimensions:")
print(dim(pb_counts))


# ----------------------------- #
# 7. Pseudobulk metadata
# ----------------------------- #

pb_meta <- meta %>%
  as.data.frame() %>%
  group_by(
    pb_id,
    donor_id,
    celltype_new = celltype_new_for_pb,
    latency_group
  ) %>%
  summarise(
    n_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(match(pb_id, colnames(pb_counts)))

if (!identical(pb_meta$pb_id, colnames(pb_counts))) {
  stop("Pseudobulk metadata and count matrix columns are not aligned.")
}

pseudobulk_metadata_path <- file.path(
  config$output_dir,
  config$pseudobulk_metadata_csv
)

write_csv(pb_meta, pseudobulk_metadata_path)


# ----------------------------- #
# 8. Filter tiny pseudobulks
# ----------------------------- #

keep_pb <- pb_meta$n_cells >= config$min_cells_per_pseudobulk

pb_counts <- pb_counts[, keep_pb, drop = FALSE]
pb_meta <- pb_meta[keep_pb, , drop = FALSE]

if (!identical(pb_meta$pb_id, colnames(pb_counts))) {
  stop("Filtered pseudobulk metadata and count matrix columns are not aligned.")
}

pseudobulk_metadata_filtered_path <- file.path(
  config$output_dir,
  config$pseudobulk_metadata_filtered_csv
)

write_csv(pb_meta, pseudobulk_metadata_filtered_path)


# ----------------------------- #
# 9. Summarise pseudobulk numbers
# ----------------------------- #

pb_summary <- pb_meta %>%
  count(celltype_new, latency_group, name = "n_pseudobulks") %>%
  pivot_wider(
    names_from = latency_group,
    values_from = n_pseudobulks,
    values_fill = 0
  ) %>%
  arrange(celltype_new)

pseudobulk_summary_path <- file.path(
  config$output_dir,
  config$pseudobulk_summary_csv
)

write_csv(pb_summary, pseudobulk_summary_path)

message("Pseudobulk counts by celltype:")
print(pb_summary)


# ----------------------------- #
# 10. edgeR pseudobulk function per celltype_new
# ----------------------------- #

run_edger_one_celltype <- function(celltype_label) {
  
  message("Running edgeR pseudobulk for: ", celltype_label)
  
  meta_ct <- pb_meta %>%
    filter(celltype_new == celltype_label) %>%
    filter(
      latency_group %in% c(
        config$latency_reference,
        config$latency_test_level
      )
    )
  
  paired_donors <- meta_ct %>%
    count(donor_id, latency_group) %>%
    group_by(donor_id) %>%
    summarise(
      n_groups = n_distinct(latency_group),
      .groups = "drop"
    ) %>%
    filter(n_groups == 2) %>%
    pull(donor_id)
  
  meta_ct <- meta_ct %>%
    filter(donor_id %in% paired_donors)
  
  if (length(unique(meta_ct$donor_id)) < config$min_paired_donors_per_celltype) {
    message(
      "Skipping ",
      celltype_label,
      ": fewer than ",
      config$min_paired_donors_per_celltype,
      " paired donors."
    )
    
    return(list(
      result = NULL,
      skipped = tibble(
        celltype_new = celltype_label,
        reason = paste0(
          "fewer_than_",
          config$min_paired_donors_per_celltype,
          "_paired_donors"
        ),
        n_paired_donors = length(unique(meta_ct$donor_id)),
        n_pseudobulks = nrow(meta_ct)
      )
    ))
  }
  
  meta_ct$latency_group <- factor(
    as.character(meta_ct$latency_group),
    levels = c(config$latency_reference, config$latency_test_level)
  )
  
  meta_ct$donor_id <- factor(meta_ct$donor_id)
  
  counts_ct <- pb_counts[, meta_ct$pb_id, drop = FALSE]
  
  if (!identical(colnames(counts_ct), meta_ct$pb_id)) {
    stop("Celltype-specific counts and metadata are not aligned for: ", celltype_label)
  }
  
  design_ct <- model.matrix(
    ~ donor_id + latency_group,
    data = meta_ct
  )
  
  if (qr(design_ct)$rank < ncol(design_ct)) {
    message("Skipping ", celltype_label, ": design matrix is not full rank.")
    
    return(list(
      result = NULL,
      skipped = tibble(
        celltype_new = celltype_label,
        reason = "design_matrix_not_full_rank",
        n_paired_donors = length(unique(meta_ct$donor_id)),
        n_pseudobulks = nrow(meta_ct)
      )
    ))
  }
  
  dge_ct <- DGEList(
    counts = counts_ct,
    group = meta_ct$latency_group
  )
  
  keep_ct <- filterByExpr(
    dge_ct,
    design = design_ct
  )
  
  dge_ct <- dge_ct[keep_ct, , keep.lib.sizes = FALSE]
  
  if (nrow(dge_ct) == 0) {
    message("Skipping ", celltype_label, ": no genes retained after filterByExpr.")
    
    return(list(
      result = NULL,
      skipped = tibble(
        celltype_new = celltype_label,
        reason = "no_genes_retained_after_filterByExpr",
        n_paired_donors = length(unique(meta_ct$donor_id)),
        n_pseudobulks = nrow(meta_ct)
      )
    ))
  }
  
  dge_ct <- calcNormFactors(dge_ct)
  dge_ct <- estimateDisp(dge_ct, design_ct)
  
  fit_ct <- glmQLFit(dge_ct, design_ct)
  
  coef_name <- paste0("latency_group", config$latency_test_level)
  
  if (!(coef_name %in% colnames(design_ct))) {
    message("Skipping ", celltype_label, ": coefficient ", coef_name, " not found.")
    
    return(list(
      result = NULL,
      skipped = tibble(
        celltype_new = celltype_label,
        reason = paste0("coefficient_not_found_", coef_name),
        n_paired_donors = length(unique(meta_ct$donor_id)),
        n_pseudobulks = nrow(meta_ct)
      )
    ))
  }
  
  res_ct <- glmQLFTest(
    fit_ct,
    coef = which(colnames(design_ct) == coef_name)
  )
  
  results_ct <- topTags(res_ct, n = Inf)$table %>%
    rownames_to_column("gene") %>%
    mutate(
      celltype_new = celltype_label,
      comparison = paste0(config$latency_test_level, "_vs_", config$latency_reference),
      n_pseudobulks = nrow(meta_ct),
      n_samples = length(unique(meta_ct$donor_id)),
      n_latent_pseudobulks = sum(meta_ct$latency_group == config$latency_test_level),
      n_lytic_pseudobulks = sum(meta_ct$latency_group == config$latency_reference),
      direction = format_direction(
        logfc = logFC,
        fdr = FDR,
        fdr_threshold = config$fdr_threshold,
        abs_logfc_threshold = config$abs_logfc_threshold
      )
    ) %>%
    arrange(FDR)
  
  output_file <- file.path(
    config$output_dir,
    paste0(
      "edgeR_pseudobulk_",
      config$latency_test_level,
      "_vs_",
      config$latency_reference,
      "_",
      safe_name(celltype_label),
      ".csv"
    )
  )
  
  write_csv(results_ct, output_file)
  
  list(
    result = results_ct,
    skipped = tibble()
  )
}


# ----------------------------- #
# 11. Run edgeR across all celltype_new groups
# ----------------------------- #

if (is.null(config$celltype_order)) {
  celltypes_to_test <- pb_meta %>%
    distinct(celltype_new) %>%
    arrange(celltype_new) %>%
    pull(celltype_new)
} else {
  celltypes_to_test <- config$celltype_order[
    config$celltype_order %in% unique(pb_meta$celltype_new)
  ]
  
  # Include any unexpected cell types at the end.
  unexpected_celltypes <- setdiff(
    sort(unique(pb_meta$celltype_new)),
    config$celltype_order
  )
  
  celltypes_to_test <- c(celltypes_to_test, unexpected_celltypes)
}

edger_run_list <- lapply(
  celltypes_to_test,
  run_edger_one_celltype
)

edger_by_celltype <- bind_rows(
  lapply(edger_run_list, `[[`, "result")
)

skipped_celltypes <- bind_rows(
  lapply(edger_run_list, `[[`, "skipped")
)

if (nrow(edger_by_celltype) == 0) {
  warning("No edgeR results were generated for any cell type.")
  
  edger_by_celltype <- tibble(
    gene = character(),
    logFC = numeric(),
    logCPM = numeric(),
    F = numeric(),
    PValue = numeric(),
    FDR = numeric(),
    celltype_new = character(),
    comparison = character(),
    n_pseudobulks = integer(),
    n_samples = integer(),
    n_latent_pseudobulks = integer(),
    n_lytic_pseudobulks = integer(),
    direction = character()
  )
}

if (nrow(skipped_celltypes) == 0) {
  skipped_celltypes <- tibble(
    celltype_new = character(),
    reason = character(),
    n_paired_donors = integer(),
    n_pseudobulks = integer()
  )
}

combined_results_path <- file.path(
  config$output_dir,
  config$combined_results_csv
)

write_csv(edger_by_celltype, combined_results_path)

skipped_celltypes_path <- file.path(
  config$output_dir,
  config$skipped_celltypes_csv
)

write_csv(skipped_celltypes, skipped_celltypes_path)

message("Top edgeR results:")
print(head(edger_by_celltype, 20))


# ----------------------------- #
# 12. Significant edgeR genes
# ----------------------------- #

sig_edger_by_celltype <- edger_by_celltype %>%
  filter(
    !is.na(FDR),
    FDR < config$fdr_threshold,
    !is.na(logFC),
    abs(logFC) > config$abs_logfc_threshold
  )

significant_results_path <- file.path(
  config$output_dir,
  config$significant_results_csv
)

write_csv(sig_edger_by_celltype, significant_results_path)

sig_edger_summary <- sig_edger_by_celltype %>%
  count(
    celltype_new,
    direction,
    name = "n_genes"
  ) %>%
  arrange(celltype_new, direction)

significant_summary_path <- file.path(
  config$output_dir,
  config$significant_summary_csv
)

write_csv(sig_edger_summary, significant_summary_path)

message("Significant gene summary:")
print(sig_edger_summary)


# ----------------------------- #
# 13. Save session information
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
    "Latency values before grouping:",
    capture.output(print(table(merged@meta.data[[config$latency_col]], useNA = "ifany"))),
    "",
    "Latency group counts:",
    capture.output(print(table(merged$latency_group, useNA = "ifany"))),
    "",
    "Pseudobulk summary:",
    capture.output(print(pb_summary)),
    "",
    "Skipped cell types:",
    capture.output(print(skipped_celltypes)),
    "",
    "Significant gene summary:",
    capture.output(print(sig_edger_summary)),
    "",
    "Output files:",
    pseudobulk_metadata_path,
    pseudobulk_metadata_filtered_path,
    pseudobulk_summary_path,
    combined_results_path,
    significant_results_path,
    significant_summary_path,
    skipped_celltypes_path,
    "",
    "Session information:",
    capture.output(sessionInfo())
  ),
  session_info_path
)


# ----------------------------- #
# 14. Completion messages
# ----------------------------- #

message("\nedgeR pseudobulk Latent vs Lytic workflow complete.")
message("Pseudobulk metadata: ", pseudobulk_metadata_path)
message("Filtered pseudobulk metadata: ", pseudobulk_metadata_filtered_path)
message("Pseudobulk summary: ", pseudobulk_summary_path)
message("Combined edgeR results: ", combined_results_path)
message("Significant edgeR results: ", significant_results_path)
message("Significant gene summary: ", significant_summary_path)
message("Skipped cell types: ", skipped_celltypes_path)
message("Session info: ", session_info_path)
