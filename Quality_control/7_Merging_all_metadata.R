#!/usr/bin/env Rscript

# ==============================================================================
# Final merge of curated single-cell Seurat objects
# ==============================================================================
#
# Purpose:
#   This script merges curated single-cell objects from multiple datasets/batches,
#   including global PBMC, EBV-enriched baseline T cells, and EBV-activated T cells.
#
#   It performs:
#     1. loading of merged per-dataset Seurat objects,
#     2. optional ADT layer joining,
#     3. ADT feature alignment across objects,
#     4. preservation of existing tetramer calls,
#     5. assignment of negative tetramer labels to non-multimer datasets,
#     6. removal of DoubletFinder metadata columns,
#     7. optional removal of selected assays before merging,
#     8. metadata column alignment,
#     9. final object merging, and
#    10. joining of sample-level metadata.
#
# Expected input:
#   Each input RDS should contain a merged Seurat object from an earlier pipeline
#   stage.
#
# Pipeline naming:
#   For datasets without ADT multimer calling, e.g. global PBMC:
#     6_SCT_merged_seurat_object.rds
#
#   For datasets with ADT multimer calling + VDJ merging, e.g. EBV datasets:
#     7_SCT_merged_seurat_object.rds
#
#   Final output:
#     8_merged_all_objects_final.rds
#
# Notes:
#   - This script does not re-call tetramers by default.
#   - Existing tetramer and tetramer_CLR metadata from the ADT multimer-calling
#     step are preserved.
#   - Global PBMC objects are assigned tetramer = "negative" unless those columns
#     already exist.
#   - Keep project-specific paths in the CONFIG section only.
#
# ==============================================================================


# ----------------------------- #
# 0. Load packages
# ----------------------------- #

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(readxl)
  library(readr)
  library(stringr)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  datasets = list(
    list(
      dataset_id = "global_pbmc_batch1",
      input_rds = "path/to/input/global_pbmc_batch1/6_SCT_merged_seurat_object.rds",
      has_multimer_calls = FALSE
    ),
    
    list(
      dataset_id = "global_pbmc_batch2",
      input_rds = "path/to/input/global_pbmc_batch2/6_SCT_merged_seurat_object.rds",
      has_multimer_calls = FALSE
    ),
    
    list(
      dataset_id = "ebv_enrichment_batch1",
      input_rds = "path/to/input/ebv_enrichment_batch1/7_SCT_merged_seurat_object.rds",
      has_multimer_calls = TRUE
    ),
    
    list(
      dataset_id = "ebv_enrichment_batch2",
      input_rds = "path/to/input/ebv_enrichment_batch2/7_SCT_merged_seurat_object.rds",
      has_multimer_calls = TRUE
    ),
    
    list(
      dataset_id = "ebv_activated_batch1",
      input_rds = "path/to/input/ebv_activated_batch1/7_SCT_merged_seurat_object.rds",
      has_multimer_calls = TRUE
    )
  ),
  
  output_rds = "path/to/output/8_merged_all_objects_final.rds",
  merge_summary_csv = "path/to/output/final_merge_summary.csv",
  failed_objects_csv = "path/to/output/final_merge_failed_objects.csv",
  metadata_join_summary_csv = "path/to/output/sample_metadata_join_summary.csv",
  session_info_file = "path/to/output/sessionInfo_final_merge.txt",
  
  # Optional sample-level metadata file.
  use_sample_metadata = TRUE,
  sample_metadata_file = "path/to/input/Sample_information.xlsx",
  sample_id_column = "id",
  sample_metadata_overwrite_existing = TRUE,
  
  # ADT assay settings.
  adt_assay = "ADT",
  align_adt_assays = TRUE,
  join_adt_layers = TRUE,
  normalize_adt_after_alignment = TRUE,
  adt_normalization_method = "CLR",
  adt_clr_margin = 2,
  
  # Tetramer metadata settings.
  tetramer_column = "tetramer",
  tetramer_clr_column = "tetramer_CLR",
  negative_tetramer_label = "negative",
  
  # Recalculate tetramer calls after ADT alignment?
  # Usually keep FALSE because tetramer calls should already have been assigned
  # in step 5_ADT_multimer_called_seurat_objects.rds.
  recompute_tetramer_calls = FALSE,
  tetramer_exclude_patterns = c(
    "HTO",
    "Hashtag",
    "^HTO_",
    "anti",
    "CD69",
    "CD137"
  ),
  
  # Metadata cleanup.
  doubletfinder_column_pattern = "^(DF\\.class|DF\\.classifications|pANN)",
  metadata_columns_to_drop = c("old.ident", "cohort"),
  
  # Assays to remove before final merge.
  # Removing SCT can make final merging simpler if SCT models differ across objects.
  assays_to_drop_before_merge = c("SCT"),
  
  verbose = TRUE
)


# ----------------------------- #
# 2. Helper functions
# ----------------------------- #

create_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}


load_merged_object <- function(dataset_cfg) {
  if (!file.exists(dataset_cfg$input_rds)) {
    stop("Input RDS file does not exist: ", dataset_cfg$input_rds)
  }
  
  obj <- readRDS(dataset_cfg$input_rds)
  
  if (!inherits(obj, "Seurat")) {
    stop("Input RDS must contain a merged Seurat object: ", dataset_cfg$input_rds)
  }
  
  obj$dataset_id <- dataset_cfg$dataset_id
  
  obj
}


get_assay_data_compat <- function(seurat_obj, assay, layer_or_slot = "counts") {
  tryCatch(
    GetAssayData(seurat_obj, assay = assay, layer = layer_or_slot),
    error = function(e) {
      GetAssayData(seurat_obj, assay = assay, slot = layer_or_slot)
    }
  )
}


create_assay_from_counts <- function(counts) {
  if ("CreateAssay5Object" %in% getNamespaceExports("Seurat")) {
    CreateAssay5Object(counts = counts)
  } else {
    CreateAssayObject(counts = counts)
  }
}


join_adt_layers_if_needed <- function(obj, cfg) {
  if (!isTRUE(cfg$join_adt_layers)) {
    return(obj)
  }
  
  if (!cfg$adt_assay %in% Assays(obj)) {
    return(obj)
  }
  
  obj[[cfg$adt_assay]] <- tryCatch(
    JoinLayers(obj[[cfg$adt_assay]]),
    error = function(e) {
      obj[[cfg$adt_assay]]
    }
  )
  
  obj
}


get_adt_features <- function(obj, cfg) {
  if (!cfg$adt_assay %in% Assays(obj)) {
    return(character())
  }
  
  rownames(obj[[cfg$adt_assay]])
}


fill_missing_adt_features <- function(obj, all_adt_features, cfg) {
  if (!cfg$adt_assay %in% Assays(obj)) {
    stop("ADT assay was not found in object: ", unique(obj$dataset_id))
  }
  
  adt_counts <- get_assay_data_compat(
    seurat_obj = obj,
    assay = cfg$adt_assay,
    layer_or_slot = "counts"
  )
  
  current_features <- rownames(adt_counts)
  missing_features <- setdiff(all_adt_features, current_features)
  
  if (length(missing_features) > 0) {
    zero_mat <- Matrix(
      0,
      nrow = length(missing_features),
      ncol = ncol(adt_counts),
      sparse = TRUE,
      dimnames = list(missing_features, colnames(adt_counts))
    )
    
    adt_counts <- rbind(adt_counts, zero_mat)
  }
  
  adt_counts <- adt_counts[all_adt_features, , drop = FALSE]
  
  obj[[cfg$adt_assay]] <- create_assay_from_counts(adt_counts)
  
  obj
}


normalize_adt <- function(obj, cfg) {
  if (!isTRUE(cfg$normalize_adt_after_alignment)) {
    return(obj)
  }
  
  if (!cfg$adt_assay %in% Assays(obj)) {
    return(obj)
  }
  
  obj <- NormalizeData(
    object = obj,
    assay = cfg$adt_assay,
    normalization.method = cfg$adt_normalization_method,
    margin = cfg$adt_clr_margin,
    verbose = cfg$verbose
  )
  
  obj
}


ensure_tetramer_metadata <- function(obj, has_multimer_calls, cfg) {
  if (!cfg$tetramer_column %in% colnames(obj@meta.data)) {
    obj[[cfg$tetramer_column]] <- cfg$negative_tetramer_label
  }
  
  if (!cfg$tetramer_clr_column %in% colnames(obj@meta.data)) {
    obj[[cfg$tetramer_clr_column]] <- 0
  }
  
  if (!isTRUE(has_multimer_calls)) {
    obj[[cfg$tetramer_column]] <- cfg$negative_tetramer_label
    obj[[cfg$tetramer_clr_column]] <- 0
  }
  
  obj
}


assign_tetramer_from_adt <- function(obj, cfg) {
  if (!cfg$adt_assay %in% Assays(obj)) {
    stop("Cannot recompute tetramer calls because ADT assay is missing.")
  }
  
  adt_data <- get_assay_data_compat(
    seurat_obj = obj,
    assay = cfg$adt_assay,
    layer_or_slot = "data"
  )
  
  features <- rownames(adt_data)
  
  exclude_idx <- rep(FALSE, length(features))
  
  for (pattern in cfg$tetramer_exclude_patterns) {
    exclude_idx <- exclude_idx | grepl(pattern, features, ignore.case = TRUE)
  }
  
  keep_features <- features[!exclude_idx]
  
  if (length(keep_features) == 0) {
    stop("No ADT features left after applying tetramer exclusion patterns.")
  }
  
  adt_mat <- as.matrix(adt_data[keep_features, , drop = FALSE])
  
  max_idx <- apply(adt_mat, 2, which.max)
  max_value <- apply(adt_mat, 2, max, na.rm = TRUE)
  
  tetramer_call <- keep_features[max_idx]
  
  all_zero <- colSums(adt_mat != 0, na.rm = TRUE) == 0
  tetramer_call[all_zero] <- cfg$negative_tetramer_label
  
  tetramer_clr <- max_value
  tetramer_clr[tetramer_call == cfg$negative_tetramer_label] <- NA_real_
  
  names(tetramer_call) <- colnames(adt_mat)
  names(tetramer_clr) <- colnames(adt_mat)
  
  obj[[cfg$tetramer_column]] <- tetramer_call[colnames(obj)]
  obj[[cfg$tetramer_clr_column]] <- tetramer_clr[colnames(obj)]
  
  obj
}


remove_doubletfinder_columns <- function(obj, cfg) {
  cols_to_remove <- grep(
    cfg$doubletfinder_column_pattern,
    colnames(obj@meta.data),
    value = TRUE
  )
  
  if (length(cols_to_remove) > 0) {
    obj@meta.data <- obj@meta.data[
      ,
      !colnames(obj@meta.data) %in% cols_to_remove,
      drop = FALSE
    ]
  }
  
  obj
}


drop_metadata_columns <- function(obj, cfg) {
  cols_to_drop <- intersect(cfg$metadata_columns_to_drop, colnames(obj@meta.data))
  
  if (length(cols_to_drop) > 0) {
    obj@meta.data <- obj@meta.data[
      ,
      !colnames(obj@meta.data) %in% cols_to_drop,
      drop = FALSE
    ]
  }
  
  obj
}


drop_assays_before_merge <- function(obj, cfg) {
  assays_to_drop <- intersect(cfg$assays_to_drop_before_merge, Assays(obj))
  
  for (assay_name in assays_to_drop) {
    obj[[assay_name]] <- NULL
  }
  
  obj
}


align_metadata_columns <- function(obj_list) {
  all_meta_cols <- Reduce(
    union,
    lapply(obj_list, function(obj) colnames(obj@meta.data))
  )
  
  obj_list <- lapply(obj_list, function(obj) {
    missing_cols <- setdiff(all_meta_cols, colnames(obj@meta.data))
    
    for (col in missing_cols) {
      obj@meta.data[[col]] <- NA
    }
    
    obj@meta.data <- obj@meta.data[, all_meta_cols, drop = FALSE]
    
    obj
  })
  
  obj_list
}


merge_object_list <- function(obj_list) {
  if (length(obj_list) == 0) {
    stop("No objects available for final merge.")
  }
  
  if (length(obj_list) == 1) {
    return(obj_list[[1]])
  }
  
  merge(
    x = obj_list[[1]],
    y = obj_list[-1]
  )
}


join_sample_metadata <- function(merged_obj, cfg) {
  if (!isTRUE(cfg$use_sample_metadata)) {
    return(list(
      merged_obj = merged_obj,
      summary = tibble(
        sample_metadata_used = FALSE,
        n_cells = ncol(merged_obj),
        n_cells_with_sample_metadata = NA_integer_,
        pct_cells_with_sample_metadata = NA_real_
      )
    ))
  }
  
  if (!file.exists(cfg$sample_metadata_file)) {
    stop("Sample metadata file does not exist: ", cfg$sample_metadata_file)
  }
  
  sample_info <- read_xlsx(cfg$sample_metadata_file)
  
  if (!cfg$sample_id_column %in% colnames(sample_info)) {
    stop(
      "Sample metadata file must contain sample ID column: ",
      cfg$sample_id_column
    )
  }
  
  if (!cfg$sample_id_column %in% colnames(merged_obj@meta.data)) {
    stop(
      "Merged Seurat metadata must contain sample ID column: ",
      cfg$sample_id_column
    )
  }
  
  sample_info_unique <- sample_info |>
    group_by(.data[[cfg$sample_id_column]]) |>
    slice(1) |>
    ungroup()
  
  meta <- merged_obj@meta.data |>
    rownames_to_column("cell")
  
  if (isTRUE(cfg$sample_metadata_overwrite_existing)) {
    overlapping_cols <- intersect(colnames(meta), colnames(sample_info_unique))
    overlapping_cols <- setdiff(overlapping_cols, c("cell", cfg$sample_id_column))
    
    if (length(overlapping_cols) > 0) {
      meta <- meta |>
        select(-all_of(overlapping_cols))
    }
  }
  
  meta_joined <- meta |>
    left_join(sample_info_unique, by = cfg$sample_id_column)
  
  n_cells <- nrow(meta_joined)
  
  sample_info_cols <- setdiff(
    colnames(sample_info_unique),
    cfg$sample_id_column
  )
  
  if (length(sample_info_cols) > 0) {
    n_cells_with_sample_metadata <- sum(
      complete.cases(meta_joined[, sample_info_cols, drop = FALSE])
    )
  } else {
    n_cells_with_sample_metadata <- n_cells
  }
  
  meta_joined <- meta_joined |>
    column_to_rownames("cell")
  
  merged_obj@meta.data <- meta_joined
  
  summary <- tibble(
    sample_metadata_used = TRUE,
    sample_metadata_file = cfg$sample_metadata_file,
    n_cells = n_cells,
    n_cells_with_sample_metadata = n_cells_with_sample_metadata,
    pct_cells_with_sample_metadata = n_cells_with_sample_metadata / n_cells * 100,
    n_unique_sample_ids_in_object = n_distinct(merged_obj@meta.data[[cfg$sample_id_column]]),
    n_unique_sample_ids_in_sample_metadata = n_distinct(sample_info_unique[[cfg$sample_id_column]])
  )
  
  list(
    merged_obj = merged_obj,
    summary = summary
  )
}


make_object_summary <- function(obj, dataset_cfg, cfg, status = "processed", message = NA_character_) {
  tibble(
    dataset_id = dataset_cfg$dataset_id,
    input_rds = dataset_cfg$input_rds,
    status = status,
    message = message,
    n_cells = if (inherits(obj, "Seurat")) ncol(obj) else NA_integer_,
    n_features_RNA = if (inherits(obj, "Seurat") && "RNA" %in% Assays(obj)) nrow(obj[["RNA"]]) else NA_integer_,
    has_ADT = if (inherits(obj, "Seurat")) cfg$adt_assay %in% Assays(obj) else FALSE,
    n_features_ADT = if (inherits(obj, "Seurat") && cfg$adt_assay %in% Assays(obj)) nrow(obj[[cfg$adt_assay]]) else NA_integer_,
    has_multimer_calls = dataset_cfg$has_multimer_calls
  )
}


# ----------------------------- #
# 3. Run final integration
# ----------------------------- #

create_parent_dir(config$output_rds)
create_parent_dir(config$merge_summary_csv)
create_parent_dir(config$failed_objects_csv)
create_parent_dir(config$metadata_join_summary_csv)
create_parent_dir(config$session_info_file)

object_list <- list()
object_summaries <- list()
failed_objects <- list()

for (dataset_cfg in config$datasets) {
  dataset_id <- dataset_cfg$dataset_id
  
  message("\n==============================")
  message("Loading dataset: ", dataset_id)
  message("==============================")
  
  result <- tryCatch(
    {
      obj <- load_merged_object(dataset_cfg)
      
      obj <- join_adt_layers_if_needed(obj, config)
      
      obj <- ensure_tetramer_metadata(
        obj = obj,
        has_multimer_calls = dataset_cfg$has_multimer_calls,
        cfg = config
      )
      
      obj <- remove_doubletfinder_columns(obj, config)
      obj <- drop_metadata_columns(obj, config)
      
      object_list[[dataset_id]] <- obj
      object_summaries[[dataset_id]] <- make_object_summary(
        obj = obj,
        dataset_cfg = dataset_cfg,
        cfg = config
      )
      
      TRUE
    },
    error = function(e) {
      message("Failed to load/process dataset: ", dataset_id, " | ", e$message)
      
      failed_objects[[dataset_id]] <<- tibble(
        dataset_id = dataset_id,
        input_rds = dataset_cfg$input_rds,
        status = "failed",
        error_message = e$message
      )
      
      object_summaries[[dataset_id]] <<- make_object_summary(
        obj = NULL,
        dataset_cfg = dataset_cfg,
        cfg = config,
        status = "failed",
        message = e$message
      )
      
      FALSE
    }
  )
}

if (length(object_list) == 0) {
  stop("No objects were successfully loaded. Cannot perform final merge.")
}


# ----------------------------- #
# 4. Align ADT assays
# ----------------------------- #

if (isTRUE(config$align_adt_assays)) {
  message("\nAligning ADT assays across objects...")
  
  adt_feature_list <- lapply(object_list, get_adt_features, cfg = config)
  adt_union <- Reduce(union, adt_feature_list)
  
  if (length(adt_union) == 0) {
    warning("No ADT features found across objects. Skipping ADT alignment.")
  } else {
    object_list <- lapply(
      object_list,
      fill_missing_adt_features,
      all_adt_features = adt_union,
      cfg = config
    )
    
    object_list <- lapply(
      object_list,
      normalize_adt,
      cfg = config
    )
  }
}


# ----------------------------- #
# 5. Optional tetramer re-calling
# ----------------------------- #

if (isTRUE(config$recompute_tetramer_calls)) {
  message("\nRecomputing tetramer calls from aligned ADT assay...")
  
  for (dataset_cfg in config$datasets) {
    dataset_id <- dataset_cfg$dataset_id
    
    if (!dataset_id %in% names(object_list)) {
      next
    }
    
    if (isTRUE(dataset_cfg$has_multimer_calls)) {
      object_list[[dataset_id]] <- assign_tetramer_from_adt(
        obj = object_list[[dataset_id]],
        cfg = config
      )
    } else {
      object_list[[dataset_id]] <- ensure_tetramer_metadata(
        obj = object_list[[dataset_id]],
        has_multimer_calls = FALSE,
        cfg = config
      )
    }
  }
}


# ----------------------------- #
# 6. Final cleanup before merging
# ----------------------------- #

message("\nCleaning assays and metadata before final merge...")

object_list <- lapply(object_list, drop_assays_before_merge, cfg = config)
object_list <- align_metadata_columns(object_list)


# ----------------------------- #
# 7. Merge all objects
# ----------------------------- #

message("\nMerging all objects...")

merged_obj <- merge_object_list(object_list)


# ----------------------------- #
# 8. Join sample-level metadata
# ----------------------------- #

message("\nJoining sample-level metadata...")

metadata_join <- join_sample_metadata(
  merged_obj = merged_obj,
  cfg = config
)

merged_obj <- metadata_join$merged_obj
metadata_join_summary <- metadata_join$summary


# ----------------------------- #
# 9. Save outputs
# ----------------------------- #

merge_summary_df <- bind_rows(object_summaries)
failed_objects_df <- bind_rows(failed_objects)

if (nrow(failed_objects_df) == 0) {
  failed_objects_df <- tibble(
    dataset_id = character(),
    input_rds = character(),
    status = character(),
    error_message = character()
  )
}

saveRDS(merged_obj, config$output_rds)

write.csv(
  merge_summary_df,
  config$merge_summary_csv,
  row.names = FALSE
)

write.csv(
  failed_objects_df,
  config$failed_objects_csv,
  row.names = FALSE
)

write.csv(
  metadata_join_summary,
  config$metadata_join_summary_csv,
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  config$session_info_file
)

message("\nFinal merge complete.")
message("Saved final merged object to: ", config$output_rds)
message("Saved merge summary to: ", config$merge_summary_csv)
message("Saved failed object log to: ", config$failed_objects_csv)
message("Saved metadata join summary to: ", config$metadata_join_summary_csv)
message("Saved session info to: ", config$session_info_file)