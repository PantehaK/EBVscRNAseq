#!/usr/bin/env Rscript

# ==============================================================================
# DoubletFinder doublet removal for a list of Seurat objects
# ==============================================================================
#
# Purpose:
#   This script loads one or more RDS files containing Seurat objects, runs
#   DoubletFinder per sample, removes predicted doublets, and saves:
#     1. a filtered Seurat object list,
#     2. post-doublet-removal QC summaries,
#     3. per-sample DoubletFinder summaries,
#     4. pK selection plots, and
#     5. a failed/skipped sample log.
#
# Expected input:
#   An .rds file containing either:
#     - a named list of Seurat objects, or
#     - a single Seurat object.
#
# Notes:
#   - This script assumes the objects have already been normalised, usually with
#     SCTransform, and contain an SCT assay by default.
#   - Adjust the configuration section below before running.
#
# ===============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(ggplot2)
  library(dplyr)
  library(tibble)
})

# ------------------------------------------------------------------------------
# 1. User configuration
# ------------------------------------------------------------------------------

# Add one list entry per dataset/batch you want to process.
# Keep paths generic in the repository; users should update these for their own
# project structure.

datasets <- list(
  list(
    dataset_id = "example_dataset",
    
    input_rds = "path/to/input/2_SCTransform_normalised_seurat_objects.rds",
    output_rds = "path/to/output/3_SCT_doublet_removed_seurat_objects.rds",
    
    qc_summary_csv = "path/to/output/post_doublet_removal_summary.csv",
    doubletfinder_summary_csv = "path/to/output/doubletfinder_summary.csv",
    failed_samples_csv = "path/to/output/failed_samples.csv",
    pk_plot_dir = "path/to/output/pK_selection_plots",
    
    # Optional filename suffix for plots, e.g. "GEM11" or "EBV2".
    plot_suffix = "",
    
    # DoubletFinder parameters.
    # Set expected_doublet_rate according to your experiment/loading chemistry.
    expected_doublet_rate = 0.08,
    pN = 0.25,
    sct = TRUE,
    
    # Dimensionality and clustering parameters used before DoubletFinder.
    assay = "SCT",
    dims_neighbors = 1:15,
    dims_umap = 1:10,
    dims_doubletfinder = 1:10,
    clustering_resolution = 0.5,
    
    # Samples below this cell count are skipped.
    min_cells = 100,
    
    # Increase if large Seurat objects trigger future/global size errors.
    future_globals_max_size_gb = 2
  )
)

# ------------------------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------------------------

create_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

create_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

safe_dataset_id <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
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

get_latest_doubletfinder_column <- function(seurat_obj) {
  doublet_cols <- grep("^DF.classifications", colnames(seurat_obj@meta.data), value = TRUE)
  
  if (length(doublet_cols) == 0) {
    return(NA_character_)
  }
  
  # Return the latest classification column if DoubletFinder was run more than once.
  tail(doublet_cols, 1)
}

make_post_qc_summary <- function(seurat_obj, sample_name) {
  tibble(
    sample = sample_name,
    nCells = ncol(seurat_obj),
    median_nCount_RNA = median(seurat_obj$nCount_RNA, na.rm = TRUE),
    median_nFeature_RNA = median(seurat_obj$nFeature_RNA, na.rm = TRUE),
    median_log10_UMI = median(log10(seurat_obj$nCount_RNA + 1), na.rm = TRUE)
  )
}

save_pk_selection_plot <- function(bcmvn, sample_name, cfg) {
  create_dir(cfg$pk_plot_dir)
  
  suffix <- cfg$plot_suffix
  suffix <- ifelse(is.null(suffix) || suffix == "", "", paste0("_", suffix))
  
  plot_path <- file.path(
    cfg$pk_plot_dir,
    paste0(safe_sample_name(sample_name), "_pK_selection", suffix, ".png")
  )
  
  bcmvn_plot <- bcmvn %>%
    mutate(pK_numeric = as.numeric(as.character(pK)))
  
  p <- ggplot(bcmvn_plot, aes(x = pK_numeric, y = BCmetric)) +
    geom_point() +
    geom_line() +
    labs(
      title = paste("pK selection:", sample_name),
      x = "pK",
      y = "BCmvn metric"
    ) +
    theme_bw(base_size = 12)
  
  ggsave(plot_path, plot = p, height = 5, width = 5, dpi = 300)
  
  plot_path
}

run_doubletfinder_one_sample <- function(seurat_obj, sample_name, cfg) {
  start_cells <- ncol(seurat_obj)
  
  if (start_cells < cfg$min_cells) {
    return(list(
      seurat_obj = seurat_obj,
      summary = tibble(
        sample = sample_name,
        status = "skipped_too_few_cells",
        start_cells = start_cells,
        retained_cells = start_cells,
        predicted_doublets = NA_integer_,
        predicted_singlets = NA_integer_,
        expected_doublet_rate = cfg$expected_doublet_rate,
        expected_doublets = NA_integer_,
        homotypic_proportion = NA_real_,
        adjusted_expected_doublets = NA_integer_,
        optimal_pK = NA_real_,
        pK_plot = NA_character_,
        message = paste0("Sample had fewer than ", cfg$min_cells, " cells.")
      )
    ))
  }
  
  if (!cfg$assay %in% Assays(seurat_obj)) {
    stop("Assay '", cfg$assay, "' was not found in sample '", sample_name, "'.")
  }
  
  DefaultAssay(seurat_obj) <- cfg$assay
  
  options(future.globals.maxSize = cfg$future_globals_max_size_gb * 1024^3)
  
  max_dims <- max(c(cfg$dims_neighbors, cfg$dims_umap, cfg$dims_doubletfinder))
  
  if (!"pca" %in% Reductions(seurat_obj)) {
    message("Running PCA for ", sample_name, " because no PCA reduction was found.")
    seurat_obj <- RunPCA(seurat_obj, npcs = max_dims, verbose = FALSE)
  }
  
  seurat_obj <- FindNeighbors(
    seurat_obj,
    dims = cfg$dims_neighbors,
    reduction = "pca",
    verbose = FALSE
  )
  
  seurat_obj <- FindClusters(
    seurat_obj,
    resolution = cfg$clustering_resolution,
    verbose = FALSE
  )
  
  seurat_obj <- RunUMAP(
    seurat_obj,
    dims = cfg$dims_umap,
    reduction = "pca",
    verbose = FALSE
  )
  
  sweep_res <- paramSweep(
    seurat_obj,
    PCs = cfg$dims_doubletfinder,
    sct = cfg$sct
  )
  
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)
  
  if (nrow(bcmvn) == 0 || all(is.na(bcmvn$BCmetric))) {
    stop("Could not identify an optimal pK for sample '", sample_name, "'.")
  }
  
  pK_value <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  pk_plot <- save_pk_selection_plot(bcmvn, sample_name, cfg)
  
  annotations <- seurat_obj@meta.data$seurat_clusters
  homotypic_proportion <- modelHomotypic(annotations)
  
  expected_doublets <- round(cfg$expected_doublet_rate * ncol(seurat_obj))
  adjusted_expected_doublets <- round(expected_doublets * (1 - homotypic_proportion))
  adjusted_expected_doublets <- max(adjusted_expected_doublets, 1)
  
  seurat_obj <- doubletFinder(
    seurat_obj,
    PCs = cfg$dims_doubletfinder,
    pN = cfg$pN,
    pK = pK_value,
    nExp = adjusted_expected_doublets,
    reuse.pANN = NULL,
    sct = cfg$sct
  )
  
  doublet_col <- get_latest_doubletfinder_column(seurat_obj)
  
  if (is.na(doublet_col)) {
    stop("DoubletFinder classification column was not found after running DoubletFinder.")
  }
  
  classifications <- seurat_obj[[doublet_col, drop = TRUE]]
  predicted_doublets <- sum(classifications == "Doublet", na.rm = TRUE)
  predicted_singlets <- sum(classifications == "Singlet", na.rm = TRUE)
  
  singlet_cells <- colnames(seurat_obj)[classifications == "Singlet"]
  seurat_obj_filtered <- subset(seurat_obj, cells = singlet_cells)
  
  list(
    seurat_obj = seurat_obj_filtered,
    summary = tibble(
      sample = sample_name,
      status = "processed",
      start_cells = start_cells,
      retained_cells = ncol(seurat_obj_filtered),
      predicted_doublets = predicted_doublets,
      predicted_singlets = predicted_singlets,
      expected_doublet_rate = cfg$expected_doublet_rate,
      expected_doublets = expected_doublets,
      homotypic_proportion = homotypic_proportion,
      adjusted_expected_doublets = adjusted_expected_doublets,
      optimal_pK = pK_value,
      pK_plot = pk_plot,
      message = NA_character_
    )
  )
}

process_dataset <- function(cfg) {
  message("\n==============================")
  message("Processing dataset: ", cfg$dataset_id)
  message("==============================")
  
  create_parent_dir(cfg$output_rds)
  create_parent_dir(cfg$qc_summary_csv)
  create_parent_dir(cfg$doubletfinder_summary_csv)
  create_parent_dir(cfg$failed_samples_csv)
  create_dir(cfg$pk_plot_dir)
  
  seurat_objects <- load_seurat_object_list(cfg$input_rds)
  
  message("Loaded ", length(seurat_objects), " Seurat object(s):")
  message(paste(names(seurat_objects), collapse = ", "))
  
  doubletfinder_summaries <- list()
  qc_summaries <- list()
  failed_samples <- list()
  
  for (sample_name in names(seurat_objects)) {
    message("\nProcessing sample: ", sample_name)
    
    result <- tryCatch(
      {
        run_doubletfinder_one_sample(seurat_objects[[sample_name]], sample_name, cfg)
      },
      error = function(e) {
        message("Failed sample: ", sample_name, " | ", e$message)
        
        failed_samples[[sample_name]] <<- tibble(
          dataset_id = cfg$dataset_id,
          sample = sample_name,
          error_message = e$message
        )
        
        list(
          seurat_obj = seurat_objects[[sample_name]],
          summary = tibble(
            sample = sample_name,
            status = "failed",
            start_cells = ncol(seurat_objects[[sample_name]]),
            retained_cells = ncol(seurat_objects[[sample_name]]),
            predicted_doublets = NA_integer_,
            predicted_singlets = NA_integer_,
            expected_doublet_rate = cfg$expected_doublet_rate,
            expected_doublets = NA_integer_,
            homotypic_proportion = NA_real_,
            adjusted_expected_doublets = NA_integer_,
            optimal_pK = NA_real_,
            pK_plot = NA_character_,
            message = e$message
          )
        )
      }
    )
    
    seurat_objects[[sample_name]] <- result$seurat_obj
    doubletfinder_summaries[[sample_name]] <- result$summary
    qc_summaries[[sample_name]] <- make_post_qc_summary(result$seurat_obj, sample_name)
  }
  
  doubletfinder_summary_df <- bind_rows(doubletfinder_summaries) %>%
    mutate(dataset_id = cfg$dataset_id, .before = sample)
  
  post_qc_summary_df <- bind_rows(qc_summaries) %>%
    mutate(dataset_id = cfg$dataset_id, .before = sample)
  
  failed_samples_df <- bind_rows(failed_samples)
  
  write.csv(doubletfinder_summary_df, cfg$doubletfinder_summary_csv, row.names = FALSE)
  write.csv(post_qc_summary_df, cfg$qc_summary_csv, row.names = FALSE)
  write.csv(failed_samples_df, cfg$failed_samples_csv, row.names = FALSE)
  
  saveRDS(seurat_objects, cfg$output_rds)
  
  message("\nSaved filtered objects to: ", cfg$output_rds)
  message("Saved post-QC summary to: ", cfg$qc_summary_csv)
  message("Saved DoubletFinder summary to: ", cfg$doubletfinder_summary_csv)
  message("Saved failed sample log to: ", cfg$failed_samples_csv)
  
  invisible(list(
    seurat_objects = seurat_objects,
    doubletfinder_summary = doubletfinder_summary_df,
    post_qc_summary = post_qc_summary_df,
    failed_samples = failed_samples_df
  ))
}

# ------------------------------------------------------------------------------
# 3. Run pipeline
# ------------------------------------------------------------------------------

results <- lapply(datasets, process_dataset)

message("\nAll DoubletFinder processing complete.")
