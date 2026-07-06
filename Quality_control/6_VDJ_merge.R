#!/usr/bin/env Rscript

# ==============================================================================
# VDJ metadata merge and Seurat object merging
# ==============================================================================
#
# Purpose:
#   This script loads one or more RDS files containing Seurat objects, searches
#   for matching 10x Genomics VDJ-T output folders, merges VDJ metadata into each
#   Seurat object, and then merges all processed Seurat objects into one object.
#
#   It saves:
#     1. Seurat object list with VDJ metadata added,
#     2. merged Seurat object,
#     3. per-sample VDJ merge summary,
#     4. failed/missing VDJ sample log, and
#     5. session information.
#
# Expected input:
#   An .rds file containing either:
#     - a named list of Seurat objects, or
#     - a single Seurat object.
#
# Expected VDJ input:
#   For each sample, this script expects a 10x VDJ-T folder containing:
#     - filtered_contig_annotations.csv
#     - clonotypes.csv
#
# Notes:
#   - The script keeps one VDJ row per cell barcode, matching the original
#     workflow where duplicated contig rows are collapsed by barcode.
#   - Samples listed in exclude_samples_regex are retained in the Seurat list and
#     final merged object, but VDJ metadata is not added for those samples.
#   - Cell names are renamed with sample IDs before final merging so that barcodes
#     are unique across samples.
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
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

# Add one list entry per dataset/batch you want to process.
# Keep paths generic in the repository. Users should update these paths for their
# own project structure.

datasets <- list(
  
  # --------------------------------------------------------------------------
  # Global PBMC batch 1
  # --------------------------------------------------------------------------
  list(
    dataset_id = "global_pbmc_batch1",
    
    input_rds = "path/to/input/4_SCT_cellcycle_filtered_seurat_objects_all_samples.rds",
    
    output_vdj_rds = "path/to/output/5_SCT_VDJ_all_samples.rds",
    output_merged_rds = "path/to/output/6_SCT_merged_all_samples.rds",
    output_dir = "path/to/output/VDJ_merge_global_batch1",
    
    vdj_base = "path/to/cellranger/output/global_pbmc_batch1",
    vdj_search_patterns = c(
      "GEM*/GEM*/outs/per_sample_outs/{sample}/vdj_t"
    ),
    
    exclude_samples_regex = NULL,
    tetramer_replacements = NULL
  ),
  
  # --------------------------------------------------------------------------
  # Global PBMC batch 2
  # --------------------------------------------------------------------------
  list(
    dataset_id = "global_pbmc_batch2",
    
    input_rds = "path/to/input/4_SCT_cellcycle_filtered_seurat_objects_all_samples_GEM11.rds",
    
    output_vdj_rds = "path/to/output/5_SCT_VDJ_all_samples_GEM11.rds",
    output_merged_rds = "path/to/output/6_SCT_merged_all_samples_GEM11.rds",
    output_dir = "path/to/output/VDJ_merge_global_batch2",
    
    vdj_base = "path/to/cellranger/output/global_pbmc_batch2",
    vdj_search_patterns = c(
      "GEM*/GEM*/outs/per_sample_outs/{sample}/vdj_t"
    ),
    
    exclude_samples_regex = NULL,
    tetramer_replacements = NULL
  ),
  
  # --------------------------------------------------------------------------
  # EBV baseline batch 1
  # Uses the output from the ADT/tetramer calling script.
  # --------------------------------------------------------------------------
  list(
    dataset_id = "ebv_baseline_batch1",
    
    input_rds = "path/to/input/5_ADT_filtered_all_samples.rds",
    
    output_vdj_rds = "path/to/output/5_SCT_VDJ_all_samples.rds",
    output_merged_rds = "path/to/output/6_SCT_merged_all_samples.rds",
    output_dir = "path/to/output/VDJ_merge_EBV_baseline_batch1",
    
    vdj_base = "path/to/cellranger/output/ebv_baseline",
    vdj_search_patterns = c(
      "GEMEBV*/GEMEBV*/outs/per_sample_outs/{sample}/vdj_t"
    ),
    
    exclude_samples_regex = "GEMEBV3NMSBE006|GEMEBV5NMSBE245|GEMEBV5MSBE443",
    
    # Correct known tetramer naming issue.
    # Named vector: names are regex patterns, values are replacements.
    tetramer_replacements = c("^RPPK" = "RPQK*")
  ),
  
  # --------------------------------------------------------------------------
  # EBV baseline batch 2
  # Uses the output from the ADT/tetramer calling script.
  # --------------------------------------------------------------------------
  list(
    dataset_id = "ebv_baseline_batch2",
    
    input_rds = "path/to/input/5_ADT_filtered_all_samples_EBV2.rds",
    
    output_vdj_rds = "path/to/output/5_SCT_VDJ_all_samples_EBV2.rds",
    output_merged_rds = "path/to/output/6_SCT_merged_all_samples_EBV2.rds",
    output_dir = "path/to/output/VDJ_merge_EBV_baseline_batch2",
    
    vdj_base = "path/to/cellranger/output/ebv_baseline",
    vdj_search_patterns = c(
      "GEMEBV*/GEMEBV*/outs/per_sample_outs/{sample}/vdj_t"
    ),
    
    exclude_samples_regex = "GEMEBV3NMSBE006|GEMEBV5NMSBE245|GEMEBV5MSBE443",
    tetramer_replacements = c("^RPPK" = "RPQK*")
  ),
  
  # --------------------------------------------------------------------------
  # EBV activated/LCL-stimulated samples
  # Uses the output from the ADT/tetramer calling script.
  # --------------------------------------------------------------------------
  list(
    dataset_id = "ebv_activated",
    
    input_rds = "path/to/input/5_ADT_filtered_all_samples.rds",
    
    output_vdj_rds = "path/to/output/5_SCT_VDJ_all_samples.rds",
    output_merged_rds = "path/to/output/6_SCT_merged_all_samples.rds",
    output_dir = "path/to/output/VDJ_merge_EBV_activated",
    
    vdj_base = "path/to/cellranger/output/ebv_activated",
    vdj_search_patterns = c(
      "EBVLCL*/EBVLCL*/outs/per_sample_outs/{sample}/vdj_t"
    ),
    
    exclude_samples_regex = NULL,
    tetramer_replacements = NULL
  )
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


normalise_dataset_config <- function(cfg) {
  
  defaults <- list(
    vdj_contig_file = "filtered_contig_annotations.csv",
    vdj_clonotype_file = "clonotypes.csv",
    
    barcode_col = "barcode",
    raw_clonotype_col = "raw_clonotype_id",
    clonotype_col = "clonotype_id",
    clonotype_cdr3_col = "cdr3s_aa",
    
    id_col = "id",
    sample_col = "sample",
    
    # If Seurat cells have already been renamed as sample_barcode, set this to
    # TRUE so the sample prefix is removed before matching to VDJ barcodes.
    strip_sample_prefix_for_vdj_match = FALSE,
    sample_prefix_separator = "_",
    
    # If multiple matching VDJ folders are found, the default is to fail the
    # sample because this is ambiguous.
    allow_multiple_vdj_folders = FALSE,
    
    # Add a prefix to all VDJ metadata columns if desired.
    # Leave as NULL to preserve the original column names used downstream.
    vdj_column_prefix = NULL,
    
    # Optional tetramer correction.
    tetramer_col = "tetramer",
    tetramer_replacements = NULL,
    
    verbose = TRUE
  )
  
  cfg <- modifyList(defaults, cfg)
  
  required_fields <- c(
    "dataset_id",
    "input_rds",
    "output_vdj_rds",
    "output_merged_rds",
    "output_dir",
    "vdj_base",
    "vdj_search_patterns"
  )
  
  for (field in required_fields) {
    if (is.null(cfg[[field]]) || length(cfg[[field]]) == 0) {
      stop("Missing required dataset configuration field: ", field)
    }
  }
  
  cfg$vdj_merge_summary_csv <- file.path(
    cfg$output_dir,
    paste0(cfg$dataset_id, "_vdj_merge_summary.csv")
  )
  
  cfg$failed_samples_csv <- file.path(
    cfg$output_dir,
    paste0(cfg$dataset_id, "_failed_vdj_samples.csv")
  )
  
  cfg$session_info_file <- file.path(
    cfg$output_dir,
    paste0(cfg$dataset_id, "_sessionInfo_VDJ_merge.txt")
  )
  
  cfg
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


is_excluded_sample <- function(sample_name, cfg) {
  
  if (is.null(cfg$exclude_samples_regex) || !nzchar(cfg$exclude_samples_regex)) {
    return(FALSE)
  }
  
  grepl(cfg$exclude_samples_regex, sample_name)
}


make_vdj_search_paths <- function(sample_name, cfg) {
  
  resolved_patterns <- gsub(
    pattern = "\\{sample\\}",
    replacement = sample_name,
    x = cfg$vdj_search_patterns
  )
  
  file.path(cfg$vdj_base, resolved_patterns)
}


find_vdj_folder <- function(sample_name, cfg) {
  
  search_paths <- make_vdj_search_paths(
    sample_name = sample_name,
    cfg = cfg
  )
  
  vdj_folders <- unlist(lapply(search_paths, Sys.glob), use.names = FALSE)
  vdj_folders <- unique(vdj_folders[dir.exists(vdj_folders)])
  
  if (length(vdj_folders) == 0) {
    stop("Could not find VDJ folder for sample: ", sample_name)
  }
  
  if (length(vdj_folders) > 1 && !isTRUE(cfg$allow_multiple_vdj_folders)) {
    stop(
      "Multiple VDJ folders found for sample '", sample_name,
      "'. Matching folders: ",
      paste(vdj_folders, collapse = "; ")
    )
  }
  
  vdj_folders[[1]]
}


read_vdj_tables <- function(vdj_folder, cfg) {
  
  vdj_file <- file.path(vdj_folder, cfg$vdj_contig_file)
  clono_file <- file.path(vdj_folder, cfg$vdj_clonotype_file)
  
  if (!file.exists(vdj_file)) {
    stop("Missing VDJ contig file: ", vdj_file)
  }
  
  if (!file.exists(clono_file)) {
    stop("Missing VDJ clonotype file: ", clono_file)
  }
  
  vdj_data <- readr::read_csv(vdj_file, show_col_types = FALSE)
  clono_data <- readr::read_csv(clono_file, show_col_types = FALSE)
  
  if (!cfg$barcode_col %in% colnames(vdj_data)) {
    stop("VDJ contig file does not contain barcode column: ", cfg$barcode_col)
  }
  
  if (
    cfg$raw_clonotype_col %in% colnames(vdj_data) &&
    !cfg$clonotype_col %in% colnames(vdj_data)
  ) {
    colnames(vdj_data)[colnames(vdj_data) == cfg$raw_clonotype_col] <- cfg$clonotype_col
  }
  
  if (!cfg$clonotype_col %in% colnames(vdj_data)) {
    stop("VDJ contig file does not contain clonotype column: ", cfg$clonotype_col)
  }
  
  if (!cfg$clonotype_col %in% colnames(clono_data)) {
    stop("Clonotype file does not contain clonotype column: ", cfg$clonotype_col)
  }
  
  if (!cfg$clonotype_cdr3_col %in% colnames(clono_data)) {
    warning("Clonotype file does not contain column: ", cfg$clonotype_cdr3_col)
    clono_data[[cfg$clonotype_cdr3_col]] <- NA_character_
  }
  
  clono_keep <- clono_data %>%
    select(
      all_of(c(cfg$clonotype_col, cfg$clonotype_cdr3_col))
    )
  
  vdj_data <- vdj_data %>%
    distinct(.data[[cfg$barcode_col]], .keep_all = TRUE) %>%
    left_join(clono_keep, by = cfg$clonotype_col)
  
  vdj_data
}


make_match_barcodes <- function(seurat_cells, sample_name, cfg) {
  
  match_barcodes <- seurat_cells
  
  if (isTRUE(cfg$strip_sample_prefix_for_vdj_match)) {
    prefix <- paste0(sample_name, cfg$sample_prefix_separator)
    match_barcodes <- sub(
      pattern = paste0("^", prefix),
      replacement = "",
      x = match_barcodes
    )
  }
  
  match_barcodes
}


prefix_vdj_columns <- function(vdj_data, cfg) {
  
  if (is.null(cfg$vdj_column_prefix) || !nzchar(cfg$vdj_column_prefix)) {
    return(vdj_data)
  }
  
  colnames(vdj_data) <- paste0(cfg$vdj_column_prefix, colnames(vdj_data))
  
  vdj_data
}


make_vdj_metadata_for_seurat <- function(seurat_obj, vdj_data, sample_name, cfg) {
  
  seurat_cells <- Cells(seurat_obj)
  
  match_barcodes <- make_match_barcodes(
    seurat_cells = seurat_cells,
    sample_name = sample_name,
    cfg = cfg
  )
  
  vdj_barcodes <- vdj_data[[cfg$barcode_col]]
  match_idx <- match(match_barcodes, vdj_barcodes)
  
  # Create an all-NA metadata table with exactly the Seurat cell names as rows.
  metadata <- as.data.frame(
    matrix(
      NA,
      nrow = length(seurat_cells),
      ncol = ncol(vdj_data)
    ),
    stringsAsFactors = FALSE
  )
  
  colnames(metadata) <- colnames(vdj_data)
  rownames(metadata) <- seurat_cells
  
  matched_cells <- which(!is.na(match_idx))
  
  if (length(matched_cells) > 0) {
    metadata[matched_cells, ] <- vdj_data[match_idx[matched_cells], , drop = FALSE]
  }
  
  metadata <- prefix_vdj_columns(metadata, cfg)
  
  list(
    metadata = metadata,
    n_vdj_barcodes = length(unique(vdj_barcodes)),
    n_cells = length(seurat_cells),
    n_cells_with_vdj = length(matched_cells),
    pct_cells_with_vdj = length(matched_cells) / length(seurat_cells) * 100
  )
}


apply_tetramer_replacements <- function(seurat_obj, cfg) {
  
  if (is.null(cfg$tetramer_replacements) || length(cfg$tetramer_replacements) == 0) {
    return(seurat_obj)
  }
  
  if (!cfg$tetramer_col %in% colnames(seurat_obj@meta.data)) {
    return(seurat_obj)
  }
  
  tetramer_values <- as.character(seurat_obj@meta.data[[cfg$tetramer_col]])
  
  for (pattern in names(cfg$tetramer_replacements)) {
    replacement <- cfg$tetramer_replacements[[pattern]]
    
    tetramer_values <- ifelse(
      grepl(pattern, tetramer_values),
      replacement,
      tetramer_values
    )
  }
  
  seurat_obj@meta.data[[cfg$tetramer_col]] <- tetramer_values
  
  seurat_obj
}


make_success_summary <- function(sample_name, cfg, vdj_folder, merge_result) {
  
  tibble(
    dataset_id = cfg$dataset_id,
    sample = sample_name,
    status = "merged",
    vdj_folder = vdj_folder,
    n_cells = merge_result$n_cells,
    n_vdj_barcodes = merge_result$n_vdj_barcodes,
    n_cells_with_vdj = merge_result$n_cells_with_vdj,
    pct_cells_with_vdj = merge_result$pct_cells_with_vdj,
    error_message = NA_character_
  )
}


make_skipped_summary <- function(sample_name, cfg) {
  
  tibble(
    dataset_id = cfg$dataset_id,
    sample = sample_name,
    status = "skipped_excluded_sample",
    vdj_folder = NA_character_,
    n_cells = NA_integer_,
    n_vdj_barcodes = NA_integer_,
    n_cells_with_vdj = NA_integer_,
    pct_cells_with_vdj = NA_real_,
    error_message = "Sample matched exclude_samples_regex; retained without VDJ metadata merge."
  )
}


make_failed_summary <- function(sample_name, cfg, error_message) {
  
  tibble(
    dataset_id = cfg$dataset_id,
    sample = sample_name,
    status = "failed",
    vdj_folder = NA_character_,
    n_cells = NA_integer_,
    n_vdj_barcodes = NA_integer_,
    n_cells_with_vdj = NA_integer_,
    pct_cells_with_vdj = NA_real_,
    error_message = error_message
  )
}


merge_vdj_one_sample <- function(seurat_obj, sample_name, cfg) {
  
  if (is_excluded_sample(sample_name, cfg)) {
    message("Skipping excluded sample: ", sample_name)
    
    seurat_obj <- apply_tetramer_replacements(
      seurat_obj = seurat_obj,
      cfg = cfg
    )
    
    return(list(
      seurat_obj = seurat_obj,
      summary = make_skipped_summary(sample_name, cfg)
    ))
  }
  
  message("Searching VDJ data for: ", sample_name)
  
  vdj_folder <- find_vdj_folder(
    sample_name = sample_name,
    cfg = cfg
  )
  
  message("Reading VDJ files from: ", vdj_folder)
  
  vdj_data <- read_vdj_tables(
    vdj_folder = vdj_folder,
    cfg = cfg
  )
  
  merge_result <- make_vdj_metadata_for_seurat(
    seurat_obj = seurat_obj,
    vdj_data = vdj_data,
    sample_name = sample_name,
    cfg = cfg
  )
  
  seurat_obj <- AddMetaData(
    object = seurat_obj,
    metadata = merge_result$metadata
  )
  
  seurat_obj <- apply_tetramer_replacements(
    seurat_obj = seurat_obj,
    cfg = cfg
  )
  
  message(
    "Successfully merged VDJ data for: ", sample_name,
    " | cells with VDJ: ", merge_result$n_cells_with_vdj,
    "/", merge_result$n_cells,
    " (", round(merge_result$pct_cells_with_vdj, 2), "%)"
  )
  
  list(
    seurat_obj = seurat_obj,
    summary = make_success_summary(
      sample_name = sample_name,
      cfg = cfg,
      vdj_folder = vdj_folder,
      merge_result = merge_result
    )
  )
}


prepare_objects_for_merge <- function(seurat_objects, cfg) {
  
  for (sample_name in names(seurat_objects)) {
    
    seurat_obj <- seurat_objects[[sample_name]]
    
    # Add sample ID columns before renaming cells.
    seurat_obj[[cfg$id_col]] <- sample_name
    
    if (!cfg$sample_col %in% colnames(seurat_obj@meta.data)) {
      seurat_obj[[cfg$sample_col]] <- sample_name
    }
    
    expected_prefix <- paste0(sample_name, cfg$sample_prefix_separator)
    
    already_prefixed <- all(startsWith(Cells(seurat_obj), expected_prefix))
    
    if (!already_prefixed) {
      seurat_obj <- RenameCells(
        object = seurat_obj,
        add.cell.id = sample_name
      )
    }
    
    # Re-add after renaming for safety.
    seurat_obj[[cfg$id_col]] <- sample_name
    
    if (!cfg$sample_col %in% colnames(seurat_obj@meta.data)) {
      seurat_obj[[cfg$sample_col]] <- sample_name
    }
    
    seurat_objects[[sample_name]] <- seurat_obj
  }
  
  seurat_objects
}


merge_seurat_object_list <- function(seurat_objects, cfg) {
  
  seurat_objects <- prepare_objects_for_merge(
    seurat_objects = seurat_objects,
    cfg = cfg
  )
  
  if (length(seurat_objects) == 1) {
    return(seurat_objects[[1]])
  }
  
  merge(
    x = seurat_objects[[1]],
    y = seurat_objects[-1],
    add.cell.ids = NULL,
    project = cfg$dataset_id
  )
}


empty_summary_table <- function() {
  
  tibble(
    dataset_id = character(),
    sample = character(),
    status = character(),
    vdj_folder = character(),
    n_cells = integer(),
    n_vdj_barcodes = integer(),
    n_cells_with_vdj = integer(),
    pct_cells_with_vdj = numeric(),
    error_message = character()
  )
}


process_dataset <- function(cfg_raw) {
  
  cfg <- normalise_dataset_config(cfg_raw)
  
  message("\n==============================")
  message("Processing dataset: ", cfg$dataset_id)
  message("==============================")
  
  create_dir(cfg$output_dir)
  create_parent_dir(cfg$output_vdj_rds)
  create_parent_dir(cfg$output_merged_rds)
  create_parent_dir(cfg$vdj_merge_summary_csv)
  create_parent_dir(cfg$failed_samples_csv)
  create_parent_dir(cfg$session_info_file)
  
  seurat_objects <- load_seurat_object_list(cfg$input_rds)
  
  message("Loaded ", length(seurat_objects), " Seurat object(s):")
  message(paste(names(seurat_objects), collapse = ", "))
  
  merge_summaries <- list()
  
  for (sample_name in names(seurat_objects)) {
    
    result <- tryCatch(
      {
        processed <- merge_vdj_one_sample(
          seurat_obj = seurat_objects[[sample_name]],
          sample_name = sample_name,
          cfg = cfg
        )
        
        seurat_objects[[sample_name]] <- processed$seurat_obj
        processed$summary
      },
      error = function(e) {
        message("Failed VDJ merge for sample: ", sample_name, " | ", e$message)
        
        make_failed_summary(
          sample_name = sample_name,
          cfg = cfg,
          error_message = e$message
        )
      }
    )
    
    merge_summaries[[sample_name]] <- result
  }
  
  merge_summary_df <- if (length(merge_summaries) > 0) {
    bind_rows(merge_summaries)
  } else {
    empty_summary_table()
  }
  
  failed_samples_df <- merge_summary_df %>%
    filter(status %in% c("failed", "skipped_excluded_sample"))
  
  saveRDS(seurat_objects, cfg$output_vdj_rds)
  
  message("Merging all processed Seurat objects into one dataset...")
  
  merged_seurat <- merge_seurat_object_list(
    seurat_objects = seurat_objects,
    cfg = cfg
  )
  
  saveRDS(merged_seurat, cfg$output_merged_rds)
  
  write.csv(merge_summary_df, cfg$vdj_merge_summary_csv, row.names = FALSE)
  write.csv(failed_samples_df, cfg$failed_samples_csv, row.names = FALSE)
  
  writeLines(
    c(
      "Dataset configuration:",
      capture.output(str(cfg)),
      "",
      "Session information:",
      capture.output(sessionInfo())
    ),
    cfg$session_info_file
  )
  
  message("\nSaved VDJ-merged Seurat object list to: ", cfg$output_vdj_rds)
  message("Saved final merged Seurat object to: ", cfg$output_merged_rds)
  message("Saved VDJ merge summary to: ", cfg$vdj_merge_summary_csv)
  message("Saved failed/skipped sample log to: ", cfg$failed_samples_csv)
  message("Saved session info to: ", cfg$session_info_file)
  
  invisible(list(
    seurat_objects = seurat_objects,
    merged_seurat = merged_seurat,
    vdj_merge_summary = merge_summary_df,
    failed_samples = failed_samples_df
  ))
}


# ----------------------------- #
# 3. Run pipeline
# ----------------------------- #

results <- lapply(datasets, process_dataset)

message("\nAll VDJ metadata merging and Seurat object merging complete.")