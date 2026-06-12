#!/usr/bin/env Rscript

# ==============================================================================
# ADT multimer normalisation and tetramer calling for Seurat objects
# ==============================================================================
#
# Purpose:
#   This script loads one or more RDS files containing Seurat objects, normalises
#   ADT features, excludes hashtag/non-multimer ADTs, assigns each cell to the
#   highest CLR-normalised multimer signal, and saves:
#     1. updated Seurat object list with tetramer metadata,
#     2. ADT ridge plots,
#     3. ADT violin plots,
#     4. per-sample tetramer call summaries,
#     5. failed sample log, and
#     6. session information.
#
# Expected input:
#   An .rds file containing either:
#     - a named list of Seurat objects, or
#     - a single Seurat object.
#
# Notes:
#   - This script assumes an ADT assay is present.
#   - Multimer identity is assigned as the ADT feature with the highest
#     CLR-normalised signal per cell.
#   - Hashtag/HTO features are excluded before calling multimers.
#   - For EBV activated datasets, activation markers such as CD69 and CD137 can
#     also be excluded using config$exclude_feature_pattern.
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
    dataset_id = "example_dataset",
    
    input_rds = "path/to/input/4_SCT_cellcycle_scored_seurat_objects.rds",
    output_rds = "path/to/output/5_ADT_multimer_called_seurat_objects.rds",
    
    plot_dir = "path/to/output/ADT_multimer_QC_plots",
    tetramer_summary_csv = "path/to/output/tetramer_call_summary.csv",
    failed_samples_csv = "path/to/output/failed_samples.csv",
    session_info_file = "path/to/output/sessionInfo_ADT_multimer_calling.txt",
    
    # ADT assay settings.
    assay = "ADT",
    
    # Features matching this pattern are excluded before multimer calling.
    # For EBV enrichment datasets:
    #   "Hashtag|HTO"
    #
    # For EBV activated datasets where CD69/CD137 are also in the ADT assay:
    #   "Hashtag|HTO|CD69|CD137"
    exclude_feature_pattern = "Hashtag|HTO",
    
    # ADT normalisation settings.
    normalization_method = "CLR",
    clr_margin = 1,
    scale_features = TRUE,
    
    # Multimer-calling metadata column names.
    tetramer_metadata_col = "tetramer",
    tetramer_clr_metadata_col = "tetramer_CLR",
    negative_label = "negative",
    
    # Optional threshold.
    # If NULL, cells are only labelled negative when all included ADT values are 0.
    # If numeric, cells with max CLR signal below this value are labelled negative.
    min_call_value = NULL,
    
    # Plot settings.
    make_plots = TRUE,
    group_by_column = "sample",
    plot_extension = "jpeg",
    plot_width = 30,
    plot_height = 30,
    ridge_plot_width = 25,
    ridge_plot_height = 25,
    plot_dpi = 300,
    plot_ncol = 5,
    ridge_y_max = 5,
    violin_point_size = 0.001,
    
    # Restore original default assay after processing each sample.
    restore_default_assay = TRUE,
    
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


safe_filename <- function(x) {
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


get_assay_data_compat <- function(seurat_obj, assay, layer_or_slot) {
  # Seurat v5 uses layer; older Seurat versions use slot.
  tryCatch(
    GetAssayData(seurat_obj, assay = assay, layer = layer_or_slot),
    error = function(e) {
      GetAssayData(seurat_obj, assay = assay, slot = layer_or_slot)
    }
  )
}


get_multimer_features <- function(seurat_obj, cfg) {
  if (!cfg$assay %in% Assays(seurat_obj)) {
    stop("Assay '", cfg$assay, "' was not found.")
  }
  
  all_features <- rownames(seurat_obj[[cfg$assay]])
  
  multimer_features <- all_features[
    !str_detect(
      string = all_features,
      pattern = regex(cfg$exclude_feature_pattern, ignore_case = TRUE)
    )
  ]
  
  if (length(multimer_features) == 0) {
    stop(
      "No ADT features remained after applying exclude_feature_pattern: ",
      cfg$exclude_feature_pattern
    )
  }
  
  multimer_features
}


choose_group_by <- function(seurat_obj, cfg) {
  if (
    !is.null(cfg$group_by_column) &&
    cfg$group_by_column %in% colnames(seurat_obj@meta.data)
  ) {
    cfg$group_by_column
  } else {
    NULL
  }
}


save_adt_ridge_plot <- function(seurat_obj, sample_name, adt_features, cfg) {
  group_by <- choose_group_by(seurat_obj, cfg)
  
  p <- tryCatch(
    {
      RidgePlot(
        object = seurat_obj,
        assay = cfg$assay,
        layer = "scale.data",
        features = adt_features,
        y.max = cfg$ridge_y_max,
        ncol = cfg$plot_ncol,
        group.by = group_by
      )
    },
    error = function(e) {
      RidgePlot(
        object = seurat_obj,
        assay = cfg$assay,
        slot = "scale.data",
        features = adt_features,
        y.max = cfg$ridge_y_max,
        ncol = cfg$plot_ncol,
        group.by = group_by
      )
    }
  ) +
    ggtitle(paste("ADT staining ridge plot -", sample_name))
  
  plot_file <- file.path(
    cfg$plot_dir,
    paste0(safe_filename(sample_name), "_ADT_ridge_plot.", cfg$plot_extension)
  )
  
  ggsave(
    filename = plot_file,
    plot = p,
    width = cfg$ridge_plot_width,
    height = cfg$ridge_plot_height,
    dpi = cfg$plot_dpi
  )
  
  plot_file
}


save_adt_violin_plot <- function(seurat_obj, sample_name, adt_features, cfg) {
  group_by <- choose_group_by(seurat_obj, cfg)
  
  p <- tryCatch(
    {
      VlnPlot(
        object = seurat_obj,
        features = adt_features,
        assay = cfg$assay,
        layer = "data",
        group.by = group_by,
        pt.size = cfg$violin_point_size,
        ncol = cfg$plot_ncol
      )
    },
    error = function(e) {
      VlnPlot(
        object = seurat_obj,
        features = adt_features,
        assay = cfg$assay,
        slot = "data",
        group.by = group_by,
        pt.size = cfg$violin_point_size,
        ncol = cfg$plot_ncol
      )
    }
  ) +
    ggtitle(paste("ADT staining violin plot -", sample_name))
  
  plot_file <- file.path(
    cfg$plot_dir,
    paste0(safe_filename(sample_name), "_ADT_violin_plot.", cfg$plot_extension)
  )
  
  ggsave(
    filename = plot_file,
    plot = p,
    width = cfg$plot_width,
    height = cfg$plot_height,
    dpi = cfg$plot_dpi
  )
  
  plot_file
}


call_tetramers <- function(seurat_obj, adt_features, cfg) {
  adt_data <- get_assay_data_compat(
    seurat_obj = seurat_obj,
    assay = cfg$assay,
    layer_or_slot = "data"
  )
  
  adt_mat <- as.matrix(adt_data[adt_features, , drop = FALSE])
  
  max_index <- apply(adt_mat, 2, which.max)
  max_value <- apply(adt_mat, 2, max, na.rm = TRUE)
  
  tetramer_call <- rownames(adt_mat)[max_index]
  
  all_zero <- colSums(adt_mat != 0, na.rm = TRUE) == 0
  tetramer_call[all_zero] <- cfg$negative_label
  
  if (!is.null(cfg$min_call_value)) {
    below_threshold <- max_value < cfg$min_call_value
    tetramer_call[below_threshold] <- cfg$negative_label
  }
  
  tetramer_clr <- max_value
  tetramer_clr[tetramer_call == cfg$negative_label] <- NA_real_
  
  names(tetramer_call) <- colnames(adt_mat)
  names(tetramer_clr) <- colnames(adt_mat)
  
  list(
    tetramer_call = tetramer_call,
    tetramer_clr = tetramer_clr
  )
}


add_tetramer_metadata <- function(seurat_obj, tetramer_calls, cfg) {
  metadata_df <- data.frame(
    tetramer = tetramer_calls$tetramer_call,
    tetramer_CLR = tetramer_calls$tetramer_clr,
    row.names = names(tetramer_calls$tetramer_call)
  )
  
  colnames(metadata_df) <- c(
    cfg$tetramer_metadata_col,
    cfg$tetramer_clr_metadata_col
  )
  
  AddMetaData(seurat_obj, metadata = metadata_df)
}


make_tetramer_summary <- function(seurat_obj, sample_name, dataset_id, adt_features, cfg) {
  metadata <- seurat_obj@meta.data
  
  tetramer_col <- cfg$tetramer_metadata_col
  clr_col <- cfg$tetramer_clr_metadata_col
  
  metadata |>
    rownames_to_column("cell_barcode") |>
    count(.data[[tetramer_col]], name = "n_cells") |>
    mutate(
      dataset_id = dataset_id,
      sample = sample_name,
      total_cells = ncol(seurat_obj),
      pct_cells = n_cells / total_cells * 100,
      n_adt_features_used = length(adt_features),
      excluded_feature_pattern = cfg$exclude_feature_pattern,
      .before = 1
    ) |>
    rename(tetramer = .data[[tetramer_col]])
}


make_sample_summary <- function(seurat_obj, sample_name, dataset_id, adt_features, cfg) {
  tetramer_col <- cfg$tetramer_metadata_col
  clr_col <- cfg$tetramer_clr_metadata_col
  
  tibble(
    dataset_id = dataset_id,
    sample = sample_name,
    status = "processed",
    n_cells = ncol(seurat_obj),
    n_adt_features_used = length(adt_features),
    n_negative = sum(seurat_obj[[tetramer_col, drop = TRUE]] == cfg$negative_label, na.rm = TRUE),
    n_positive = sum(seurat_obj[[tetramer_col, drop = TRUE]] != cfg$negative_label, na.rm = TRUE),
    pct_positive = n_positive / n_cells * 100,
    median_positive_CLR = median(
      seurat_obj[[clr_col, drop = TRUE]],
      na.rm = TRUE
    )
  )
}


process_one_sample <- function(seurat_obj, sample_name, cfg) {
  message("Processing sample: ", sample_name)
  
  original_default_assay <- DefaultAssay(seurat_obj)
  
  if (!cfg$assay %in% Assays(seurat_obj)) {
    stop("Assay '", cfg$assay, "' was not found in sample '", sample_name, "'.")
  }
  
  DefaultAssay(seurat_obj) <- cfg$assay
  
  adt_features <- get_multimer_features(seurat_obj, cfg)
  
  seurat_obj <- NormalizeData(
    object = seurat_obj,
    assay = cfg$assay,
    normalization.method = cfg$normalization_method,
    margin = cfg$clr_margin,
    features = adt_features,
    verbose = cfg$verbose
  )
  
  if (isTRUE(cfg$scale_features)) {
    seurat_obj <- ScaleData(
      object = seurat_obj,
      assay = cfg$assay,
      features = adt_features,
      verbose = cfg$verbose
    )
  }
  
  ridge_plot_file <- NA_character_
  violin_plot_file <- NA_character_
  
  if (isTRUE(cfg$make_plots)) {
    ridge_plot_file <- save_adt_ridge_plot(
      seurat_obj = seurat_obj,
      sample_name = sample_name,
      adt_features = adt_features,
      cfg = cfg
    )
    
    violin_plot_file <- save_adt_violin_plot(
      seurat_obj = seurat_obj,
      sample_name = sample_name,
      adt_features = adt_features,
      cfg = cfg
    )
  }
  
  tetramer_calls <- call_tetramers(
    seurat_obj = seurat_obj,
    adt_features = adt_features,
    cfg = cfg
  )
  
  seurat_obj <- add_tetramer_metadata(
    seurat_obj = seurat_obj,
    tetramer_calls = tetramer_calls,
    cfg = cfg
  )
  
  if (isTRUE(cfg$restore_default_assay) && original_default_assay %in% Assays(seurat_obj)) {
    DefaultAssay(seurat_obj) <- original_default_assay
  }
  
  sample_summary <- make_sample_summary(
    seurat_obj = seurat_obj,
    sample_name = sample_name,
    dataset_id = cfg$dataset_id,
    adt_features = adt_features,
    cfg = cfg
  ) |>
    mutate(
      ridge_plot = ridge_plot_file,
      violin_plot = violin_plot_file
    )
  
  tetramer_summary <- make_tetramer_summary(
    seurat_obj = seurat_obj,
    sample_name = sample_name,
    dataset_id = cfg$dataset_id,
    adt_features = adt_features,
    cfg = cfg
  )
  
  list(
    seurat_obj = seurat_obj,
    sample_summary = sample_summary,
    tetramer_summary = tetramer_summary
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
  
  seurat_list <- load_seurat_object_list(cfg$input_rds)
  
  message("Loaded ", length(seurat_list), " Seurat object(s):")
  message(paste(names(seurat_list), collapse = ", "))
  
  sample_summaries <- list()
  tetramer_summaries <- list()
  failed_samples <- list()
  
  for (sample_name in names(seurat_list)) {
    result <- tryCatch(
      {
        processed <- process_one_sample(
          seurat_obj = seurat_list[[sample_name]],
          sample_name = sample_name,
          cfg = cfg
        )
        
        seurat_list[[sample_name]] <- processed$seurat_obj
        sample_summaries[[sample_name]] <- processed$sample_summary
        tetramer_summaries[[sample_name]] <- processed$tetramer_summary
        
        TRUE
      },
      error = function(e) {
        message("Failed sample: ", sample_name, " | ", e$message)
        
        failed_samples[[sample_name]] <<- tibble(
          dataset_id = cfg$dataset_id,
          sample = sample_name,
          error_message = e$message
        )
        
        sample_summaries[[sample_name]] <<- tibble(
          dataset_id = cfg$dataset_id,
          sample = sample_name,
          status = "failed",
          n_cells = ncol(seurat_list[[sample_name]]),
          n_adt_features_used = NA_integer_,
          n_negative = NA_integer_,
          n_positive = NA_integer_,
          pct_positive = NA_real_,
          median_positive_CLR = NA_real_,
          ridge_plot = NA_character_,
          violin_plot = NA_character_
        )
        
        FALSE
      }
    )
  }
  
  sample_summary_df <- bind_rows(sample_summaries)
  tetramer_summary_df <- bind_rows(tetramer_summaries)
  failed_samples_df <- bind_rows(failed_samples)
  
  if (nrow(failed_samples_df) == 0) {
    failed_samples_df <- tibble(
      dataset_id = character(),
      sample = character(),
      error_message = character()
    )
  }
  
  saveRDS(seurat_list, cfg$output_rds)
  
  write.csv(
    sample_summary_df,
    cfg$tetramer_summary_csv,
    row.names = FALSE
  )
  
  write.csv(
    failed_samples_df,
    cfg$failed_samples_csv,
    row.names = FALSE
  )
  
  write.csv(
    tetramer_summary_df,
    sub("\\.csv$", "_by_tetramer.csv", cfg$tetramer_summary_csv),
    row.names = FALSE
  )
  
  writeLines(
    capture.output(sessionInfo()),
    cfg$session_info_file
  )
  
  message("\nSaved ADT multimer-called objects to: ", cfg$output_rds)
  message("Saved sample summary to: ", cfg$tetramer_summary_csv)
  message(
    "Saved tetramer-level summary to: ",
    sub("\\.csv$", "_by_tetramer.csv", cfg$tetramer_summary_csv)
  )
  message("Saved failed sample log to: ", cfg$failed_samples_csv)
  message("Saved session info to: ", cfg$session_info_file)
  
  invisible(list(
    seurat_list = seurat_list,
    sample_summary = sample_summary_df,
    tetramer_summary = tetramer_summary_df,
    failed_samples = failed_samples_df
  ))
}


# ----------------------------- #
# 3. Run pipeline
# ----------------------------- #

results <- lapply(datasets, process_dataset)

message("\nAll ADT multimer calling complete.")