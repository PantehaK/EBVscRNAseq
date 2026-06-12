#!/usr/bin/env Rscript

# ==============================================================================
# Merge VDJ metadata into Seurat objects and merge samples
# ==============================================================================
#
# Purpose:
#   This script loads one or more RDS files containing Seurat objects, searches
#   Cell Ranger VDJ output directories for each sample, merges TCR clonotype
#   metadata into each Seurat object, and saves:
#     1. a Seurat object list with VDJ metadata,
#     2. a merged Seurat object with unique sample-prefixed cell barcodes,
#     3. a per-sample VDJ merge summary,
#     4. a failed/skipped sample log, and
#     5. session information.
#
# Expected input:
#   An .rds file containing either:
#     - a named list of Seurat objects, or
#     - a single Seurat object.
#
# Expected VDJ files:
#   For each sample, this script searches for:
#     - filtered_contig_annotations.csv
#     - clonotypes.csv
#
# Notes:
#   - This script assumes Cell Ranger VDJ-T output.
#   - Cell barcodes in VDJ files must match the Seurat cell names before sample
#     prefixing.
#   - Samples are renamed only immediately before merging.
#   - Keep project-specific paths in the CONFIG section only.
#
# ==============================================================================


# ----------------------------- #
# 0. Load packages
# ----------------------------- #

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
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
    
    # If this dataset has ADT multimer calling, use:
    #   input_rds = "path/to/input/5_ADT_multimer_called_seurat_objects.rds"
    #
    # If this dataset does not have ADT multimer calling, use:
    #   input_rds = "path/to/input/4_SCT_cellcycle_scored_seurat_objects.rds"
    
    input_rds = "path/to/input/5_ADT_multimer_called_seurat_objects.rds",
    
    output_vdj_rds = "path/to/output/6_SCT_VDJ_seurat_objects.rds",
    output_merged_rds = "path/to/output/7_SCT_merged_seurat_object.rds",
    
    vdj_summary_csv = "path/to/output/VDJ_merge_summary.csv",
    failed_samples_csv = "path/to/output/failed_samples.csv",
    session_info_file = "path/to/output/sessionInfo_VDJ_merge.txt",
    
    # Base directory containing Cell Ranger VDJ outputs.
    vdj_base = "path/to/cellranger_vdj_outputs",
    
    # Glob pattern used to find each sample's VDJ folder.
    # Examples:
    #   "GEM*/GEM*/outs/per_sample_outs/{sample}/vdj_t"
    #   "GEMEBV*/GEMEBV*/outs/per_sample_outs/{sample}/vdj_t"
    #   "EBVLCL*/EBVLCL*/outs/per_sample_outs/{sample}/vdj_t"
    vdj_folder_pattern = "GEM*/GEM*/outs/per_sample_outs/{sample}/vdj_t",
    
    vdj_contig_file = "filtered_contig_annotations.csv",
    clonotype_file = "clonotypes.csv",
    
    # Optional samples to skip.
    excluded_samples = character(),
    
    # Metadata columns.
    sample_id_column = "id",
    
    # Optional correction for old/mislabelled tetramer names.
    # Set apply_tetramer_name_fix = TRUE only for datasets that need it.
    apply_tetramer_name_fix = FALSE,
    tetramer_column = "tetramer",
    tetramer_fix_pattern = "^RPPK",
    tetramer_fix_replacement = "RPQK*",
    
    # If multiple VDJ folders match a sample, should the sample fail?
    fail_if_multiple_vdj_folders = TRUE
  )
)


# ----------------------------- #
# 2. Helper functions
# ----------------------------- #

create_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}


safe_sample_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
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


build_vdj_search_pattern <- function(vdj_base, vdj_folder_pattern, sample_name) {
  pattern <- str_replace_all(
    vdj_folder_pattern,
    fixed("{sample}"),
    sample_name
  )
  
  file.path(vdj_base, pattern)
}


find_vdj_folder <- function(sample_name, cfg) {
  search_pattern <- build_vdj_search_pattern(
    vdj_base = cfg$vdj_base,
    vdj_folder_pattern = cfg$vdj_folder_pattern,
    sample_name = sample_name
  )
  
  vdj_folders <- Sys.glob(search_pattern)
  vdj_folders <- vdj_folders[dir.exists(vdj_folders)]
  
  if (length(vdj_folders) == 0) {
    stop("Could not find VDJ folder using pattern: ", search_pattern)
  }
  
  if (length(vdj_folders) > 1 && isTRUE(cfg$fail_if_multiple_vdj_folders)) {
    stop(
      "Multiple VDJ folders found for sample '", sample_name, "': ",
      paste(vdj_folders, collapse = "; ")
    )
  }
  
  vdj_folders[[1]]
}


read_and_format_vdj_metadata <- function(vdj_folder, cfg) {
  vdj_file <- file.path(vdj_folder, cfg$vdj_contig_file)
  clonotype_file <- file.path(vdj_folder, cfg$clonotype_file)
  
  if (!file.exists(vdj_file)) {
    stop("Missing VDJ contig file: ", vdj_file)
  }
  
  if (!file.exists(clonotype_file)) {
    stop("Missing clonotype file: ", clonotype_file)
  }
  
  vdj_data <- read_csv(vdj_file, show_col_types = FALSE)
  clonotype_data <- read_csv(clonotype_file, show_col_types = FALSE)
  
  required_vdj_cols <- c("barcode", "raw_clonotype_id")
  missing_vdj_cols <- setdiff(required_vdj_cols, colnames(vdj_data))
  
  if (length(missing_vdj_cols) > 0) {
    stop(
      "VDJ file is missing required column(s): ",
      paste(missing_vdj_cols, collapse = ", ")
    )
  }
  
  required_clonotype_cols <- c("clonotype_id", "cdr3s_aa")
  missing_clonotype_cols <- setdiff(required_clonotype_cols, colnames(clonotype_data))
  
  if (length(missing_clonotype_cols) > 0) {
    stop(
      "Clonotype file is missing required column(s): ",
      paste(missing_clonotype_cols, collapse = ", ")
    )
  }
  
  vdj_data <- vdj_data |>
    distinct(barcode, .keep_all = TRUE) |>
    rename(clonotype_id = raw_clonotype_id) |>
    left_join(
      clonotype_data |>
        select(clonotype_id, cdr3s_aa),
      by = "clonotype_id"
    )
  
  vdj_data <- as.data.frame(vdj_data)
  rownames(vdj_data) <- vdj_data$barcode
  
  vdj_data
}


align_vdj_metadata_to_seurat <- function(vdj_data, seurat_obj) {
  seurat_cells <- Cells(seurat_obj)
  
  missing_barcodes <- setdiff(seurat_cells, rownames(vdj_data))
  
  if (length(missing_barcodes) > 0) {
    missing_data <- as.data.frame(
      matrix(
        NA,
        nrow = length(missing_barcodes),
        ncol = ncol(vdj_data)
      )
    )
    
    colnames(missing_data) <- colnames(vdj_data)
    rownames(missing_data) <- missing_barcodes
    missing_data$barcode <- missing_barcodes
    
    vdj_data <- bind_rows(
      vdj_data |> rownames_to_column("rownames_tmp"),
      missing_data |> rownames_to_column("rownames_tmp")
    )
    
    rownames(vdj_data) <- vdj_data$rownames_tmp
    vdj_data$rownames_tmp <- NULL
  }
  
  # Keep only Seurat cells and force matching order.
  vdj_data <- vdj_data[seurat_cells, , drop = FALSE]
  
  vdj_data
}


merge_vdj_into_sample <- function(seurat_obj, sample_name, cfg) {
  message("Searching VDJ data for: ", sample_name)
  
  vdj_folder <- find_vdj_folder(sample_name, cfg)
  
  message("Reading VDJ files for: ", sample_name)
  
  vdj_data <- read_and_format_vdj_metadata(vdj_folder, cfg)
  
  n_vdj_barcodes <- nrow(vdj_data)
  n_cells <- ncol(seurat_obj)
  n_matched_cells <- sum(Cells(seurat_obj) %in% rownames(vdj_data))
  
  vdj_data <- align_vdj_metadata_to_seurat(vdj_data, seurat_obj)
  
  seurat_obj <- AddMetaData(seurat_obj, metadata = vdj_data)
  
  message("Successfully merged VDJ data for: ", sample_name)
  
  list(
    seurat_obj = seurat_obj,
    summary = tibble(
      sample = sample_name,
      status = "processed",
      vdj_folder = vdj_folder,
      n_cells = n_cells,
      n_vdj_barcodes = n_vdj_barcodes,
      n_cells_with_vdj_match = n_matched_cells,
      pct_cells_with_vdj_match = n_matched_cells / n_cells * 100
    )
  )
}


apply_optional_tetramer_fix <- function(seurat_obj, cfg) {
  if (!isTRUE(cfg$apply_tetramer_name_fix)) {
    return(seurat_obj)
  }
  
  if (!cfg$tetramer_column %in% colnames(seurat_obj@meta.data)) {
    return(seurat_obj)
  }
  
  values <- seurat_obj[[cfg$tetramer_column, drop = TRUE]]
  
  values[str_detect(values, cfg$tetramer_fix_pattern)] <- cfg$tetramer_fix_replacement
  
  seurat_obj[[cfg$tetramer_column]] <- values
  
  seurat_obj
}


merge_seurat_object_list <- function(seurat_objects, cfg) {
  if (length(seurat_objects) == 0) {
    stop("No Seurat objects available for merging.")
  }
  
  for (sample_name in names(seurat_objects)) {
    seurat_obj <- seurat_objects[[sample_name]]
    
    seurat_obj[[cfg$sample_id_column]] <- sample_name
    
    seurat_obj <- RenameCells(
      object = seurat_obj,
      add.cell.id = safe_sample_name(sample_name)
    )
    
    seurat_obj[[cfg$sample_id_column]] <- sample_name
    
    seurat_objects[[sample_name]] <- seurat_obj
  }
  
  if (length(seurat_objects) == 1) {
    return(seurat_objects[[1]])
  }
  
  merge(
    x = seurat_objects[[1]],
    y = seurat_objects[-1],
    add.cell.ids = NULL
  )
}


process_dataset <- function(cfg) {
  message("\n==============================")
  message("Processing dataset: ", cfg$dataset_id)
  message("==============================")
  
  create_parent_dir(cfg$output_vdj_rds)
  create_parent_dir(cfg$output_merged_rds)
  create_parent_dir(cfg$vdj_summary_csv)
  create_parent_dir(cfg$failed_samples_csv)
  create_parent_dir(cfg$session_info_file)
  
  seurat_objects <- load_seurat_object_list(cfg$input_rds)
  
  message("Loaded ", length(seurat_objects), " Seurat object(s):")
  message(paste(names(seurat_objects), collapse = ", "))
  
  vdj_summaries <- list()
  failed_samples <- list()
  
  for (sample_name in names(seurat_objects)) {
    if (sample_name %in% cfg$excluded_samples) {
      message("Skipping excluded sample: ", sample_name)
      
      failed_samples[[sample_name]] <- tibble(
        dataset_id = cfg$dataset_id,
        sample = sample_name,
        status = "skipped_excluded_sample",
        error_message = NA_character_
      )
      
      next
    }
    
    result <- tryCatch(
      {
        merged <- merge_vdj_into_sample(
          seurat_obj = seurat_objects[[sample_name]],
          sample_name = sample_name,
          cfg = cfg
        )
        
        seurat_objects[[sample_name]] <- merged$seurat_obj
        seurat_objects[[sample_name]] <- apply_optional_tetramer_fix(
          seurat_obj = seurat_objects[[sample_name]],
          cfg = cfg
        )
        
        vdj_summaries[[sample_name]] <- merged$summary |>
          mutate(dataset_id = cfg$dataset_id, .before = sample)
        
        TRUE
      },
      error = function(e) {
        message("Failed sample: ", sample_name, " | ", e$message)
        
        failed_samples[[sample_name]] <<- tibble(
          dataset_id = cfg$dataset_id,
          sample = sample_name,
          status = "failed",
          error_message = e$message
        )
        
        vdj_summaries[[sample_name]] <<- tibble(
          dataset_id = cfg$dataset_id,
          sample = sample_name,
          status = "failed",
          vdj_folder = NA_character_,
          n_cells = ncol(seurat_objects[[sample_name]]),
          n_vdj_barcodes = NA_integer_,
          n_cells_with_vdj_match = NA_integer_,
          pct_cells_with_vdj_match = NA_real_
        )
        
        FALSE
      }
    )
  }
  
  vdj_summary_df <- bind_rows(vdj_summaries)
  failed_samples_df <- bind_rows(failed_samples)
  
  if (nrow(failed_samples_df) == 0) {
    failed_samples_df <- tibble(
      dataset_id = character(),
      sample = character(),
      status = character(),
      error_message = character()
    )
  }
  
  saveRDS(seurat_objects, cfg$output_vdj_rds)
  
  message("Merging all processed Seurat objects into one dataset...")
  
  merged_seurat <- merge_seurat_object_list(seurat_objects, cfg)
  
  saveRDS(merged_seurat, cfg$output_merged_rds)
  
  write.csv(vdj_summary_df, cfg$vdj_summary_csv, row.names = FALSE)
  write.csv(failed_samples_df, cfg$failed_samples_csv, row.names = FALSE)
  
  writeLines(
    capture.output(sessionInfo()),
    cfg$session_info_file
  )
  
  message("\nSaved VDJ-annotated object list to: ", cfg$output_vdj_rds)
  message("Saved merged Seurat object to: ", cfg$output_merged_rds)
  message("Saved VDJ summary to: ", cfg$vdj_summary_csv)
  message("Saved failed/skipped sample log to: ", cfg$failed_samples_csv)
  message("Saved session info to: ", cfg$session_info_file)
  
  invisible(list(
    seurat_objects = seurat_objects,
    merged_seurat = merged_seurat,
    vdj_summary = vdj_summary_df,
    failed_samples = failed_samples_df
  ))
}


# ----------------------------- #
# 3. Run pipeline
# ----------------------------- #

results <- lapply(datasets, process_dataset)

message("\nAll VDJ merging and Seurat object merging complete.")