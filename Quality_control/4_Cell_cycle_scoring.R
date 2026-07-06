#!/usr/bin/env Rscript

# ==============================================================================
# Cell cycle scoring for Seurat objects
# ==============================================================================
#
# Purpose:
#   This script loads one or more RDS files containing Seurat objects, performs
#   cell cycle scoring using Seurat's built-in cell cycle gene lists, and saves:
#     1. updated Seurat object list with S.Score, G2M.Score and Phase metadata,
#     2. per-sample cell cycle summary,
#     3. failed sample log, and
#     4. session information.
#
# Expected input:
#   An .rds file containing either:
#     - a named list of Seurat objects, or
#     - a single Seurat object.
#
# Notes:
#   - This script scores cell cycle phase only.
#   - It does not remove cells.
#   - It does not regress out cell cycle effects.
#   - By default, scoring is performed on the RNA assay after log-normalisation,
#     following the Seurat cell cycle scoring vignette.
#
# ==============================================================================


# ----------------------------- #
# 0. Load packages
# ----------------------------- #

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

# Add one list entry per dataset/batch you want to process.
# Keep paths generic in the repository; users should update these for their own
# project structure.

datasets <- list(
  list(
    dataset_id = "example_dataset",
    
    input_rds = "path/to/input/3_SCT_doublet_removed_seurat_objects.rds",
    output_rds = "path/to/output/4_SCT_cellcycle_scored_seurat_objects.rds",
    
    cellcycle_summary_csv = "path/to/output/cellcycle_summary.csv",
    failed_samples_csv = "path/to/output/failed_samples.csv",
    session_info_file = "path/to/output/sessionInfo_cellcycle.txt",
    
    # Assay used for cell cycle scoring.
    assay = "RNA",
    
    # CellCycleScoring expects normalised data.
    normalize_before_scoring = TRUE,
    normalization_method = "LogNormalize",
    scale_factor = 10000,
    
    # Whether to set cell identities to cell cycle phase.
    set_ident = FALSE,
    
    # Restore the original default assay after scoring.
    restore_default_assay = TRUE,
    
    # Minimum number of matched cell cycle genes required.
    min_s_genes = 10,
    min_g2m_genes = 10,
    
    verbose = TRUE
  )
)


# ----------------------------- #
# 2. Helper functions
# ----------------------------- #

create_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}


load_seurat_object_list <- function(input_rds) {
  if (!file.exists(input_rds)) {
    stop("Input RDS file does not exist: ", input_rds)
  }
  
  x <- readRDS(input_rds)
  
  if (inherits(x, "Seurat")) {
    x <- list(sample_1 = x)
  }
  
  if (!is.list(x) || !all(vapply(x, inherits, logical(1), what = "Seurat"))) {
    stop("Input RDS must contain either a Seurat object or a list of Seurat objects.")
  }
  
  if (is.null(names(x)) || any(names(x) == "")) {
    names(x) <- paste0("sample_", seq_along(x))
  }
  
  x
}


match_features_to_object <- function(features, object_features) {
  # Case-insensitive feature matching while preserving the feature names
  # actually present in the Seurat object.
  
  feature_upper <- toupper(features)
  object_upper <- toupper(object_features)
  
  matched_indices <- match(feature_upper, object_upper)
  matched_features <- object_features[matched_indices[!is.na(matched_indices)]]
  
  unique(matched_features)
}


score_cellcycle_one_sample <- function(seurat_obj, sample_name, cfg, s_genes, g2m_genes) {
  message("Processing cell cycle scoring for: ", sample_name)
  
  if (!cfg$assay %in% Assays(seurat_obj)) {
    stop("Assay '", cfg$assay, "' was not found in sample '", sample_name, "'.")
  }
  
  original_default_assay <- DefaultAssay(seurat_obj)
  DefaultAssay(seurat_obj) <- cfg$assay
  
  object_features <- rownames(seurat_obj[[cfg$assay]])
  
  matched_s_genes <- match_features_to_object(
    features = s_genes,
    object_features = object_features
  )
  
  matched_g2m_genes <- match_features_to_object(
    features = g2m_genes,
    object_features = object_features
  )
  
  if (length(matched_s_genes) < cfg$min_s_genes) {
    stop(
      "Too few S-phase genes found in sample '", sample_name,
      "'. Found ", length(matched_s_genes), "."
    )
  }
  
  if (length(matched_g2m_genes) < cfg$min_g2m_genes) {
    stop(
      "Too few G2M-phase genes found in sample '", sample_name,
      "'. Found ", length(matched_g2m_genes), "."
    )
  }
  
  if (isTRUE(cfg$normalize_before_scoring)) {
    seurat_obj <- NormalizeData(
      object = seurat_obj,
      assay = cfg$assay,
      normalization.method = cfg$normalization_method,
      scale.factor = cfg$scale_factor,
      verbose = cfg$verbose
    )
  }
  
  seurat_obj <- CellCycleScoring(
    object = seurat_obj,
    s.features = matched_s_genes,
    g2m.features = matched_g2m_genes,
    set.ident = cfg$set_ident
  )
  
  if (isTRUE(cfg$restore_default_assay) && original_default_assay %in% Assays(seurat_obj)) {
    DefaultAssay(seurat_obj) <- original_default_assay
  }
  
  list(
    seurat_obj = seurat_obj,
    n_s_genes = length(matched_s_genes),
    n_g2m_genes = length(matched_g2m_genes)
  )
}


make_cellcycle_summary <- function(seurat_obj, sample_name, dataset_id, n_s_genes, n_g2m_genes) {
  phase_counts <- table(
    factor(seurat_obj$Phase, levels = c("G1", "S", "G2M"))
  )
  
  n_cells <- ncol(seurat_obj)
  
  tibble(
    dataset_id = dataset_id,
    sample = sample_name,
    status = "processed",
    n_cells = n_cells,
    n_s_genes_used = n_s_genes,
    n_g2m_genes_used = n_g2m_genes,
    median_S_score = median(seurat_obj$S.Score, na.rm = TRUE),
    median_G2M_score = median(seurat_obj$G2M.Score, na.rm = TRUE),
    n_G1 = as.integer(phase_counts[["G1"]]),
    n_S = as.integer(phase_counts[["S"]]),
    n_G2M = as.integer(phase_counts[["G2M"]]),
    pct_G1 = as.numeric(phase_counts[["G1"]]) / n_cells * 100,
    pct_S = as.numeric(phase_counts[["S"]]) / n_cells * 100,
    pct_G2M = as.numeric(phase_counts[["G2M"]]) / n_cells * 100
  )
}


process_dataset <- function(cfg) {
  message("\n==============================")
  message("Processing dataset: ", cfg$dataset_id)
  message("==============================")
  
  create_parent_dir(cfg$output_rds)
  create_parent_dir(cfg$cellcycle_summary_csv)
  create_parent_dir(cfg$failed_samples_csv)
  create_parent_dir(cfg$session_info_file)
  
  seurat_objects <- load_seurat_object_list(cfg$input_rds)
  
  message("Loaded ", length(seurat_objects), " Seurat object(s):")
  message(paste(names(seurat_objects), collapse = ", "))
  
  # Load Seurat's built-in cell cycle gene lists.
  cc_genes <- Seurat::cc.genes
  s_genes <- cc_genes$s.genes
  g2m_genes <- cc_genes$g2m.genes
  
  cellcycle_summaries <- list()
  failed_samples <- list()
  
  for (sample_name in names(seurat_objects)) {
    result <- tryCatch(
      {
        scored <- score_cellcycle_one_sample(
          seurat_obj = seurat_objects[[sample_name]],
          sample_name = sample_name,
          cfg = cfg,
          s_genes = s_genes,
          g2m_genes = g2m_genes
        )
        
        seurat_objects[[sample_name]] <- scored$seurat_obj
        
        cellcycle_summaries[[sample_name]] <- make_cellcycle_summary(
          seurat_obj = scored$seurat_obj,
          sample_name = sample_name,
          dataset_id = cfg$dataset_id,
          n_s_genes = scored$n_s_genes,
          n_g2m_genes = scored$n_g2m_genes
        )
        
        TRUE
      },
      error = function(e) {
        message("Failed sample: ", sample_name, " | ", e$message)
        
        failed_samples[[sample_name]] <<- tibble(
          dataset_id = cfg$dataset_id,
          sample = sample_name,
          error_message = e$message
        )
        
        cellcycle_summaries[[sample_name]] <<- tibble(
          dataset_id = cfg$dataset_id,
          sample = sample_name,
          status = "failed",
          n_cells = ncol(seurat_objects[[sample_name]]),
          n_s_genes_used = NA_integer_,
          n_g2m_genes_used = NA_integer_,
          median_S_score = NA_real_,
          median_G2M_score = NA_real_,
          n_G1 = NA_integer_,
          n_S = NA_integer_,
          n_G2M = NA_integer_,
          pct_G1 = NA_real_,
          pct_S = NA_real_,
          pct_G2M = NA_real_
        )
        
        FALSE
      }
    )
  }
  
  cellcycle_summary_df <- bind_rows(cellcycle_summaries)
  failed_samples_df <- bind_rows(failed_samples)
  
  if (nrow(failed_samples_df) == 0) {
    failed_samples_df <- tibble(
      dataset_id = character(),
      sample = character(),
      error_message = character()
    )
  }
  
  saveRDS(seurat_objects, cfg$output_rds)
  write.csv(cellcycle_summary_df, cfg$cellcycle_summary_csv, row.names = FALSE)
  write.csv(failed_samples_df, cfg$failed_samples_csv, row.names = FALSE)
  
  writeLines(
    capture.output(sessionInfo()),
    cfg$session_info_file
  )
  
  message("\nSaved cell cycle scored objects to: ", cfg$output_rds)
  message("Saved cell cycle summary to: ", cfg$cellcycle_summary_csv)
  message("Saved failed sample log to: ", cfg$failed_samples_csv)
  message("Saved session info to: ", cfg$session_info_file)
  
  invisible(list(
    seurat_objects = seurat_objects,
    cellcycle_summary = cellcycle_summary_df,
    failed_samples = failed_samples_df
  ))
}


# ----------------------------- #
# 3. Run pipeline
# ----------------------------- #

results <- lapply(datasets, process_dataset)

message("\nAll cell cycle scoring complete.")