#!/usr/bin/env Rscript

# ==============================================================================
# ADT normalisation and tetramer calling for Seurat objects
# ==============================================================================
#
# Purpose:
#   This script loads one or more RDS files containing Seurat objects, performs
#   CLR normalisation and scaling of ADT features, calls the dominant tetramer
#   signal per cell, and saves:
#     1. updated Seurat object list with tetramer and tetramer_CLR metadata,
#     2. per-sample tetramer count summaries,
#     3. ADT ridge plots,
#     4. ADT violin plots,
#     5. failed sample log, and
#     6. session information.
#
# Expected input:
#   An .rds file containing either:
#     - a named list of Seurat objects, or
#     - a single Seurat object.
#
# Notes:
#   - This script is intended to be run after cell cycle scoring/filtering.
#   - ADT features matching Hashtag, HTO, CD69 or CD137 are excluded before
#     tetramer calling by default.
#   - Tetramer identity is assigned as the ADT feature with the highest
#     CLR-normalised signal per cell.
#   - Cells with all-zero tetramer ADT signal are assigned as "negative".
#   - By default, no positive-signal threshold is applied. A threshold can be
#     added in the dataset configuration if required.
#
# ==============================================================================


# ----------------------------- #
# 0. Load packages
# ----------------------------- #

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(purrr)
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
    dataset_id = "activated_EBV",
    
    input_rds = "path/to/input/activated/4_SCT_cellcycle_filtered_seurat_objects_all_samples.rds",
    output_rds = "path/to/output/activated/5_ADT_filtered_all_samples.rds",
    
    plot_dir = "path/to/output/activated/ADT_plots",
    tetramer_summary_csv = "path/to/output/activated/tetramer_summary.csv",
    failed_samples_csv = "path/to/output/activated/failed_samples.csv",
    session_info_file = "path/to/output/activated/sessionInfo_ADT_tetramer_calling.txt",
    
    assay = "ADT",
    group_by = "sample",
    
    # Features matching this pattern are excluded from tetramer calling.
    # Adjust if your experiment uses different ADT/hashtag naming.
    exclude_feature_pattern = "Hashtag|HTO|CD69|CD137",
    
    normalization_method = "CLR",
    ridge_y_max = 5,
    plot_ncol = 5,
    save_plots = TRUE,
    
    # Set to NULL to preserve the original behaviour:
    # call the highest ADT feature unless all features are zero.
    # Example threshold use: positive_call_threshold = 1
    positive_call_threshold = NULL,
    
    verbose = TRUE
  ),
  
  list(
    dataset_id = "baseline_EBV_1",
    
    input_rds = "path/to/input/baseline/4_SCT_cellcycle_filtered_seurat_objects_all_samples.rds",
    output_rds = "path/to/output/baseline/5_ADT_filtered_all_samples.rds",
    
    plot_dir = "path/to/output/baseline/ADT_plots",
    tetramer_summary_csv = "path/to/output/baseline/tetramer_summary.csv",
    failed_samples_csv = "path/to/output/baseline/failed_samples.csv",
    session_info_file = "path/to/output/baseline/sessionInfo_ADT_tetramer_calling.txt",
    
    assay = "ADT",
    group_by = "sample",
    exclude_feature_pattern = "Hashtag|HTO|CD69|CD137",
    normalization_method = "CLR",
    ridge_y_max = 5,
    plot_ncol = 5,
    save_plots = TRUE,
    positive_call_threshold = NULL,
    verbose = TRUE
  ),
  
  list(
    dataset_id = "baseline_EBV_2",
    
    input_rds = "path/to/input/baseline/4_SCT_cellcycle_filtered_seurat_objects_all_samples_EBV2.rds",
    output_rds = "path/to/output/baseline/5_ADT_filtered_all_samples_EBV2.rds",
    
    plot_dir = "path/to/output/baseline/ADT_plots_EBV2",
    tetramer_summary_csv = "path/to/output/baseline/tetramer_summary_EBV2.csv",
    failed_samples_csv = "path/to/output/baseline/failed_samples_EBV2.csv",
    session_info_file = "path/to/output/baseline/sessionInfo_ADT_tetramer_calling_EBV2.txt",
    
    assay = "ADT",
    group_by = "sample",
    exclude_feature_pattern = "Hashtag|HTO|CD69|CD137",
    normalization_method = "CLR",
    ridge_y_max = 5,
    plot_ncol = 5,
    save_plots = TRUE,
    positive_call_threshold = NULL,
    verbose = TRUE
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


sanitize_filename <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9_\\-]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "")
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


get_assay_data_layer <- function(seurat_obj, assay, layer_name) {
  # Seurat v5 uses layer; Seurat v4 used slot.
  # This wrapper keeps the script more portable.
  
  out <- tryCatch(
    {
      GetAssayData(
        object = seurat_obj,
        assay = assay,
        layer = layer_name
      )
    },
    error = function(e) {
      GetAssayData(
        object = seurat_obj,
        assay = assay,
        slot = layer_name
      )
    }
  )
  
  out
}


make_ridge_plot <- function(seurat_obj, features, cfg, sample_name) {
  args <- list(
    object = seurat_obj,
    assay = cfg$assay,
    features = features,
    y.max = cfg$ridge_y_max,
    ncol = cfg$plot_ncol,
    group.by = cfg$group_by
  )
  
  if (packageVersion("Seurat") >= "5.0.0") {
    args$layer <- "scale.data"
  } else {
    args$slot <- "scale.data"
  }
  
  do.call(RidgePlot, args) +
    ggtitle(paste("ADT staining ridge plot -", sample_name))
}


make_violin_plot <- function(seurat_obj, features, cfg, sample_name) {
  args <- list(
    object = seurat_obj,
    assay = cfg$assay,
    features = features,
    group.by = cfg$group_by,
    pt.size = 0.001,
    ncol = cfg$plot_ncol
  )
  
  if (packageVersion("Seurat") >= "5.0.0") {
    args$layer <- "data"
  } else {
    args$slot <- "data"
  }
  
  do.call(VlnPlot, args) +
    ggtitle(paste("ADT staining violin plot -", sample_name))
}


get_tetramer_features <- function(seurat_obj, cfg, sample_name) {
  if (!cfg$assay %in% Assays(seurat_obj)) {
    stop("Assay '", cfg$assay, "' was not found in sample '", sample_name, "'.")
  }
  
  adt_all <- rownames(seurat_obj[[cfg$assay]])
  
  adt_features <- adt_all[
    !grepl(
      pattern = cfg$exclude_feature_pattern,
      x = adt_all,
      ignore.case = TRUE
    )
  ]
  
  if (length(adt_features) == 0) {
    stop(
      "No ADT features remained after filtering in sample '",
      sample_name,
      "'. Check exclude_feature_pattern."
    )
  }
  
  adt_features
}


call_tetramer_from_adt <- function(adt_mat, positive_call_threshold = NULL) {
  apply(adt_mat, 2, function(x) {
    max_value <- max(x, na.rm = TRUE)
    
    if (all(x == 0, na.rm = TRUE)) {
      return("negative")
    }
    
    if (!is.null(positive_call_threshold) && max_value < positive_call_threshold) {
      return("negative")
    }
    
    names(x)[which.max(x)]
  })
}


add_tetramer_clr_metadata <- function(seurat_obj, cfg) {
  adt_data <- get_assay_data_layer(
    seurat_obj = seurat_obj,
    assay = cfg$assay,
    layer_name = "data"
  )
  
  tetramer_clr <- map2_dbl(
    .x = colnames(adt_data),
    .y = seurat_obj$tetramer,
    .f = function(cell, tetramer_feature) {
      if (is.na(tetramer_feature) || tetramer_feature == "negative") {
        return(NA_real_)
      }
      
      as.numeric(adt_data[tetramer_feature, cell])
    }
  )
  
  seurat_obj$tetramer_CLR <- tetramer_clr
  
  seurat_obj
}


make_tetramer_summary <- function(seurat_obj, sample_name, dataset_id) {
  meta_df <- seurat_obj@meta.data %>%
    rownames_to_column("barcode") %>%
    mutate(
      tetramer = if_else(
        is.na(tetramer) | tetramer == "",
        "negative",
        as.character(tetramer)
      )
    )
  
  total_cells <- nrow(meta_df)
  
  meta_df %>%
    count(tetramer, name = "n_cells") %>%
    mutate(
      dataset_id = dataset_id,
      sample = sample_name,
      total_cells = total_cells,
      pct_cells = n_cells / total_cells * 100
    ) %>%
    select(
      dataset_id,
      sample,
      tetramer,
      n_cells,
      total_cells,
      pct_cells
    ) %>%
    arrange(sample, desc(n_cells))
}


process_one_sample <- function(seurat_obj, sample_name, cfg) {
  message("Processing ADT/tetramer calling for: ", sample_name)
  
  original_default_assay <- DefaultAssay(seurat_obj)
  DefaultAssay(seurat_obj) <- cfg$assay
  
  adt_features <- get_tetramer_features(
    seurat_obj = seurat_obj,
    cfg = cfg,
    sample_name = sample_name
  )
  
  seurat_obj <- NormalizeData(
    object = seurat_obj,
    assay = cfg$assay,
    normalization.method = cfg$normalization_method,
    features = adt_features,
    verbose = cfg$verbose
  )
  
  seurat_obj <- ScaleData(
    object = seurat_obj,
    assay = cfg$assay,
    features = adt_features,
    verbose = cfg$verbose
  )
  
  if (isTRUE(cfg$save_plots)) {
    sample_file <- sanitize_filename(sample_name)
    
    ridge_plot <- make_ridge_plot(
      seurat_obj = seurat_obj,
      features = adt_features,
      cfg = cfg,
      sample_name = sample_name
    )
    
    ggsave(
      filename = file.path(cfg$plot_dir, paste0(sample_file, "_ADT_ridge_plot.jpeg")),
      plot = ridge_plot,
      width = 25,
      height = 25,
      units = "in",
      dpi = 300
    )
    
    violin_plot <- make_violin_plot(
      seurat_obj = seurat_obj,
      features = adt_features,
      cfg = cfg,
      sample_name = sample_name
    )
    
    ggsave(
      filename = file.path(cfg$plot_dir, paste0(sample_file, "_ADT_violin_plot.jpeg")),
      plot = violin_plot,
      width = 30,
      height = 30,
      units = "in",
      dpi = 300
    )
  }
  
  adt_data <- get_assay_data_layer(
    seurat_obj = seurat_obj,
    assay = cfg$assay,
    layer_name = "data"
  )
  
  adt_mat <- as.matrix(adt_data[adt_features, , drop = FALSE])
  
  tetramer_call <- call_tetramer_from_adt(
    adt_mat = adt_mat,
    positive_call_threshold = cfg$positive_call_threshold
  )
  
  seurat_obj <- AddMetaData(
    object = seurat_obj,
    metadata = tetramer_call,
    col.name = "tetramer"
  )
  
  seurat_obj <- add_tetramer_clr_metadata(
    seurat_obj = seurat_obj,
    cfg = cfg
  )
  
  if (original_default_assay %in% Assays(seurat_obj)) {
    DefaultAssay(seurat_obj) <- original_default_assay
  }
  
  list(
    seurat_obj = seurat_obj,
    adt_features = adt_features
  )
}


process_dataset <- function(cfg) {
  message("\n==============================")
  message("Processing dataset: ", cfg$dataset_id)
  message("==============================")
  
  create_parent_dir(cfg$output_rds)
  create_parent_dir(cfg$tetramer_summary_csv)
  create_parent_dir(cfg$failed_samples_csv)
  create_parent_dir(cfg$session_info_file)
  create_dir(cfg$plot_dir)
  
  seurat_objects <- load_seurat_object_list(cfg$input_rds)
  
  message("Loaded ", length(seurat_objects), " Seurat object(s):")
  message(paste(names(seurat_objects), collapse = ", "))
  
  tetramer_summaries <- list()
  failed_samples <- list()
  
  for (sample_name in names(seurat_objects)) {
    result <- tryCatch(
      {
        processed <- process_one_sample(
          seurat_obj = seurat_objects[[sample_name]],
          sample_name = sample_name,
          cfg = cfg
        )
        
        seurat_objects[[sample_name]] <- processed$seurat_obj
        
        tetramer_summaries[[sample_name]] <- make_tetramer_summary(
          seurat_obj = processed$seurat_obj,
          sample_name = sample_name,
          dataset_id = cfg$dataset_id
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
        
        tetramer_summaries[[sample_name]] <<- tibble(
          dataset_id = cfg$dataset_id,
          sample = sample_name,
          tetramer = NA_character_,
          n_cells = NA_integer_,
          total_cells = ncol(seurat_objects[[sample_name]]),
          pct_cells = NA_real_
        )
        
        FALSE
      }
    )
  }
  
  tetramer_summary_df <- bind_rows(tetramer_summaries)
  failed_samples_df <- bind_rows(failed_samples)
  
  if (nrow(failed_samples_df) == 0) {
    failed_samples_df <- tibble(
      dataset_id = character(),
      sample = character(),
      error_message = character()
    )
  }
  
  saveRDS(seurat_objects, cfg$output_rds)
  write.csv(tetramer_summary_df, cfg$tetramer_summary_csv, row.names = FALSE)
  write.csv(failed_samples_df, cfg$failed_samples_csv, row.names = FALSE)
  
  writeLines(
    capture.output(sessionInfo()),
    cfg$session_info_file
  )
  
  message("\nSaved ADT/tetramer-called objects to: ", cfg$output_rds)
  message("Saved tetramer summary to: ", cfg$tetramer_summary_csv)
  message("Saved failed sample log to: ", cfg$failed_samples_csv)
  message("Saved session info to: ", cfg$session_info_file)
  
  invisible(list(
    seurat_objects = seurat_objects,
    tetramer_summary = tetramer_summary_df,
    failed_samples = failed_samples_df
  ))
}


# ----------------------------- #
# 3. Run pipeline
# ----------------------------- #

results <- lapply(datasets, process_dataset)

message("\nAll ADT normalisation and tetramer calling complete.")