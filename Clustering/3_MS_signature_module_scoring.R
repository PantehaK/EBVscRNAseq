#!/usr/bin/env Rscript

# ==============================================================================
# MS stimulated-latent signature module scoring in activated EBV CD8+ T cells
# ==============================================================================
#
# Purpose:
#   This script scores cells using a gene signature derived from pseudobulk
#   differential expression comparing MS versus Control stimulated latent
#   EBV-specific T cells.
#
#   It saves:
#     1. a Seurat object with the signature score added to metadata,
#     2. cell-level signature metadata,
#     3. per-sample mean signature scores by latent group,
#     4. group-level summary statistics,
#     5. missing/present gene list,
#     6. session information.
#
# Expected input:
#   A clustered activated EBV CD8+ T cell Seurat object, e.g.
#     14_activated_EBV_clustered_module_scored.rds
#
# Notes:
#   - AddModuleScore() is run on the RNA assay by default.
#   - The output score column is renamed to MS_stimLat_signature_score.
#   - If your object uses "non-EBV" instead of "Non-Latent", edit
#     config$latent_groups_to_keep.
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
  library(readr)
  library(ggplot2)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  dataset_id = "activated_ebv",
  
  input_rds = "path/to/input/14_activated_EBV_clustered_module_scored.rds",
  output_rds = "path/to/output/15_activated_EBV_MS_stimLat_signature_scored.rds",
  
  output_dir = "path/to/output/Publication_data/EBV_activated/MS_stimLat_signature",
  
  per_sample_mean_csv = "path/to/output/Publication_data/EBV_activated/MS_stimLat_signature/mean_MS_stimLat_signature_by_sample.csv",
  cell_level_csv = "path/to/output/Publication_data/EBV_activated/MS_stimLat_signature/cell_level_MS_stimLat_signature.csv",
  group_summary_csv = "path/to/output/Publication_data/EBV_activated/MS_stimLat_signature/group_summary_MS_stimLat_signature.csv",
  gene_presence_csv = "path/to/output/Publication_data/EBV_activated/MS_stimLat_signature/MS_stimLat_signature_gene_presence.csv",
  boxplot_png = "path/to/output/Publication_data/EBV_activated/MS_stimLat_signature/MS_stimLat_signature_by_group.png",
  session_info_file = "path/to/output/Publication_data/EBV_activated/MS_stimLat_signature/sessionInfo_MS_stimLat_signature.txt",
  
  assay = "RNA",
  
  module_name = "MS_stimLat_signature",
  score_column = "MS_stimLat_signature_score",
  
  group_col = "latent_group",
  sample_col = "sample",
  cohort_col = "cohort",
  
  latent_groups_to_keep = c(
    "Stimulated Latent",
    "Baseline Latent",
    "Non-Latent"
  ),
  
  signature_genes = c(
    "FTH1", "JAK1", "ITM2B", "RPL13", "SMCHD1", "IL7R", "HLA-B", "PTMA",
    "HERC5", "TOMM7", "MT-CYB", "H3-3A", "HLA-F", "RPS10", "MT-ND4L",
    "HNRNPA1", "FYN", "PFN1", "RPS29", "MT-CO3", "RPS15", "RPS18",
    "SLFN5", "ZBTB20", "RPLP2", "MT-ND3", "ARPC2", "IFI6", "CD3D",
    "RPS27", "HLA-C", "CXCR4", "RPS27A", "CALM1", "CCSER2", "ANKRD12",
    "VIM", "RPS9", "CDC42SE2", "EEF1G", "RPS6", "PDE3B", "EEF1A1",
    "RPS16", "YWHAB", "RPL30", "SARAF", "SRSF7"
  ),
  
  make_boxplot = TRUE,
  plot_width = 9,
  plot_height = 5,
  plot_dpi = 300,
  
  verbose = TRUE
)


# ----------------------------- #
# 2. Helper functions
# ----------------------------- #

create_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}


create_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}


setup_dirs <- function(cfg) {
  create_dir(cfg$output_dir)
  create_parent_dir(cfg$output_rds)
  create_parent_dir(cfg$per_sample_mean_csv)
  create_parent_dir(cfg$cell_level_csv)
  create_parent_dir(cfg$group_summary_csv)
  create_parent_dir(cfg$gene_presence_csv)
  create_parent_dir(cfg$boxplot_png)
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


remove_existing_signature_scores <- function(obj, cfg) {
  cols_to_remove <- grep(
    paste0("^", cfg$module_name),
    colnames(obj@meta.data),
    value = TRUE
  )
  
  cols_to_remove <- union(cols_to_remove, cfg$score_column)
  cols_to_remove <- intersect(cols_to_remove, colnames(obj@meta.data))
  
  if (length(cols_to_remove) > 0) {
    obj@meta.data <- obj@meta.data |>
      select(-all_of(cols_to_remove))
  }
  
  obj
}


add_signature_score <- function(obj, cfg) {
  if (!cfg$assay %in% Assays(obj)) {
    stop("Assay not found in object: ", cfg$assay)
  }
  
  DefaultAssay(obj) <- cfg$assay
  
  genes_present <- intersect(cfg$signature_genes, rownames(obj))
  genes_missing <- setdiff(cfg$signature_genes, genes_present)
  
  gene_presence <- tibble(
    gene = cfg$signature_genes,
    present_in_object = gene %in% genes_present
  )
  
  write.csv(
    gene_presence,
    cfg$gene_presence_csv,
    row.names = FALSE
  )
  
  if (length(genes_present) == 0) {
    stop("None of the signature genes were found in the object.")
  }
  
  if (length(genes_missing) > 0) {
    message(
      "Missing signature genes: ",
      paste(genes_missing, collapse = ", ")
    )
  }
  
  obj <- remove_existing_signature_scores(obj, cfg)
  
  obj <- AddModuleScore(
    object = obj,
    features = list(genes_present),
    assay = cfg$assay,
    name = cfg$module_name
  )
  
  score_cols <- grep(
    paste0("^", cfg$module_name),
    colnames(obj@meta.data),
    value = TRUE
  )
  
  newest_score_col <- tail(score_cols, 1)
  
  obj[[cfg$score_column]] <- obj[[newest_score_col, drop = TRUE]]
  
  if (newest_score_col != cfg$score_column) {
    obj@meta.data <- obj@meta.data |>
      select(-all_of(newest_score_col))
  }
  
  obj
}


make_signature_dataframes <- function(obj, cfg) {
  required_cols <- c(
    cfg$sample_col,
    cfg$cohort_col,
    cfg$group_col,
    cfg$score_column
  )
  
  missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required metadata column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  cell_df <- obj@meta.data |>
    rownames_to_column("cell_barcode") |>
    filter(.data[[cfg$group_col]] %in% cfg$latent_groups_to_keep) |>
    transmute(
      cell_barcode = cell_barcode,
      sample = .data[[cfg$sample_col]],
      cohort = .data[[cfg$cohort_col]],
      latent_group = .data[[cfg$group_col]],
      MS_stimLat_signature_score = .data[[cfg$score_column]]
    )
  
  per_sample_mean <- cell_df |>
    group_by(sample, cohort, latent_group) |>
    summarise(
      mean_MS_stimLat_signature_score = mean(MS_stimLat_signature_score, na.rm = TRUE),
      median_MS_stimLat_signature_score = median(MS_stimLat_signature_score, na.rm = TRUE),
      n_cells = n(),
      .groups = "drop"
    )
  
  group_summary <- per_sample_mean |>
    group_by(cohort, latent_group) |>
    summarise(
      n_samples = n_distinct(sample),
      mean_of_sample_means = mean(mean_MS_stimLat_signature_score, na.rm = TRUE),
      median_of_sample_means = median(mean_MS_stimLat_signature_score, na.rm = TRUE),
      sd_of_sample_means = sd(mean_MS_stimLat_signature_score, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(
    cell_df = cell_df,
    per_sample_mean = per_sample_mean,
    group_summary = group_summary
  )
}


save_signature_boxplot <- function(per_sample_mean, cfg) {
  if (!isTRUE(cfg$make_boxplot)) {
    return(invisible(NULL))
  }
  
  p <- ggplot(
    per_sample_mean,
    aes(
      x = latent_group,
      y = mean_MS_stimLat_signature_score,
      fill = cohort
    )
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(
      aes(group = cohort),
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
      size = 2,
      alpha = 0.8
    ) +
    theme_classic(base_size = 13) +
    labs(
      x = NULL,
      y = "Mean MS stimulated-latent signature score",
      fill = "Cohort"
    ) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      axis.title.y = element_text(face = "bold")
    )
  
  ggsave(
    plot = p,
    filename = cfg$boxplot_png,
    width = cfg$plot_width,
    height = cfg$plot_height,
    dpi = cfg$plot_dpi,
    bg = "transparent"
  )
  
  invisible(p)
}


# ----------------------------- #
# 3. Run pipeline
# ----------------------------- #

setup_dirs(config)

message("Loading activated EBV object...")
merged <- load_seurat_object(config$input_rds)

message("Adding MS stimulated-latent signature score...")
merged <- add_signature_score(merged, config)

message("Creating signature summary tables...")
signature_tables <- make_signature_dataframes(merged, config)

write.csv(
  signature_tables$cell_df,
  config$cell_level_csv,
  row.names = FALSE
)

write.csv(
  signature_tables$per_sample_mean,
  config$per_sample_mean_csv,
  row.names = FALSE
)

write.csv(
  signature_tables$group_summary,
  config$group_summary_csv,
  row.names = FALSE
)

message("Saving signature plot...")
save_signature_boxplot(
  per_sample_mean = signature_tables$per_sample_mean,
  cfg = config
)

message("Saving scored Seurat object...")
saveRDS(
  merged,
  config$output_rds
)

writeLines(
  capture.output(sessionInfo()),
  config$session_info_file
)

message("\nMS stimulated-latent signature scoring complete.")
message("Saved scored object to: ", config$output_rds)
message("Saved per-sample means to: ", config$per_sample_mean_csv)
message("Saved group summary to: ", config$group_summary_csv)
message("Saved session info to: ", config$session_info_file)