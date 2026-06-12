#!/usr/bin/env Rscript

# ==============================================================================
# QC filtering, immune receptor gene removal, SCTransform normalisation, and PCA
# for a named list of Seurat objects.
#
# Expected input:
#   - An RDS file containing either:
#       1) a named list of Seurat objects, or
#       2) a single Seurat object
#
# Main outputs:
#   - Filtered + SCTransformed Seurat object list as RDS
#   - Post-QC sample summary CSV
#   - Per-sample QC violin plots
#   - Per-sample variable feature plots
#   - Failed sample log, if any samples fail
#
# Notes:
#   - This script is intentionally generic for repository use.
#   - Edit the "User configuration" section before running.
#   - Designed for human scRNA-seq gene symbols by default.
# ==============================================================================


# -----------------------------#
# 0. Load packages
# -----------------------------#

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "ggplot2",
  "dplyr",
  "tibble",
  "sctransform"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Please install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(sctransform)
})


# -----------------------------#
# 1. User configuration
# -----------------------------#

# Edit this list for your project.
# You can add multiple datasets if you want to process several RDS files with
# different thresholds or output names.

datasets <- list(
  list(
    dataset_id = "example_dataset",
    
    # Input RDS containing a named list of Seurat objects.
    input_rds = "path/to/input/1_QC_seurat_objects.rds",
    
    # Output files/directories.
    output_rds = "path/to/output/2_SCTransform_normalised_seurat_objects.rds",
    qc_plot_dir = "path/to/output/qc_violin_plots",
    variable_feature_plot_dir = "path/to/output/variable_feature_plots",
    post_qc_summary_csv = "path/to/output/post_qc_summary.csv",
    failed_samples_csv = "path/to/output/failed_samples.csv",
    
    # Optional suffix added to plot filenames.
    plot_suffix = "",
    
    # QC thresholds.
    # Set any maximum/minimum to Inf/-Inf if you do not want to filter on it.
    qc_thresholds = list(
      percent_mt_max = 5,
      percent_ribo_max = 40,
      percent_hb_max = 0.1,
      nFeature_RNA_min = 200,
      nFeature_RNA_max = 3000,
      nCount_RNA_min = 200,
      nCount_RNA_max = 10000
    ),
    
    # QC gene patterns.
    # These defaults assume human gene symbols.
    qc_gene_patterns = list(
      mitochondrial = "^MT-",
      ribosomal = "^RPS|^RPL",
      haemoglobin = "^HB[AB]"
    ),
    
    # Immune receptor gene filtering.
    remove_immune_receptor_genes = TRUE,
    immune_receptor_gene_patterns = c(
      "^TRA", "^TRB", "^TRD", "^TRG",
      "^IGKV", "^IGKJ", "^IGKC",
      "^IGLV", "^IGLJ", "^IGLC",
      "^IGHV", "^IGHD", "^IGHJ",
      "^IGHM$", "^IGHD$", "^IGHG", "^IGHA", "^IGHE"
    ),
    
    # SCTransform settings.
    sctransform = list(
      vars_to_regress = c("percent.mt"),
      variable_features_n = 2000,
      ncells = 3000,
      return_only_var_genes = FALSE,
      seed_use = 1234,
      verbose = FALSE
    ),
    
    # PCA settings.
    run_pca = TRUE,
    max_pcs = 50
  )
)

# Increase this if SCTransform fails because of object size.
options(future.globals.maxSize = 2000 * 1024^2)


# -----------------------------#
# 2. Helper functions
# -----------------------------#

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


create_dir_if_missing <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}


safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}


get_counts_matrix <- function(seurat_obj, assay = "RNA") {
  # Seurat v5 uses layers; older Seurat versions use slots.
  counts <- tryCatch(
    GetAssayData(seurat_obj, assay = assay, layer = "counts"),
    error = function(e_layer) {
      tryCatch(
        GetAssayData(seurat_obj, assay = assay, slot = "counts"),
        error = function(e_slot) {
          stop(
            "Could not retrieve counts matrix for assay '", assay, "'.\n",
            "Layer error: ", e_layer$message, "\n",
            "Slot error: ", e_slot$message
          )
        }
      )
    }
  )
  
  counts
}


create_assay_from_counts <- function(counts) {
  # Prefer Assay5 when available, but remain compatible with older Seurat.
  if ("CreateAssay5Object" %in% getNamespaceExports("SeuratObject")) {
    SeuratObject::CreateAssay5Object(counts = counts)
  } else {
    Seurat::CreateAssayObject(counts = counts)
  }
}


calculate_feature_percentage <- function(seurat_obj, pattern, col_name, assay = "RNA") {
  counts <- get_counts_matrix(seurat_obj, assay = assay)
  matching_genes <- grep(pattern, rownames(counts), value = TRUE)
  
  total_counts <- Matrix::colSums(counts)
  
  if (length(matching_genes) == 0) {
    pct <- rep(0, ncol(counts))
  } else {
    feature_counts <- Matrix::colSums(counts[matching_genes, , drop = FALSE])
    pct <- feature_counts / total_counts * 100
    pct[is.na(pct) | is.infinite(pct)] <- 0
  }
  
  names(pct) <- colnames(counts)
  seurat_obj[[col_name]] <- pct
  
  seurat_obj
}


add_qc_metrics <- function(seurat_obj, qc_gene_patterns, assay = "RNA") {
  DefaultAssay(seurat_obj) <- assay
  
  seurat_obj <- calculate_feature_percentage(
    seurat_obj = seurat_obj,
    pattern = qc_gene_patterns$mitochondrial,
    col_name = "percent.mt",
    assay = assay
  )
  
  seurat_obj <- calculate_feature_percentage(
    seurat_obj = seurat_obj,
    pattern = qc_gene_patterns$ribosomal,
    col_name = "percent.ribo",
    assay = assay
  )
  
  seurat_obj <- calculate_feature_percentage(
    seurat_obj = seurat_obj,
    pattern = qc_gene_patterns$haemoglobin,
    col_name = "percent.hb",
    assay = assay
  )
  
  seurat_obj$log10_UMI <- log10(seurat_obj$nCount_RNA + 1)
  
  seurat_obj
}


qc_filter_object <- function(seurat_obj, qc_thresholds, sample_name) {
  required_columns <- c(
    "percent.mt", "percent.ribo", "percent.hb",
    "nFeature_RNA", "nCount_RNA"
  )
  
  missing_columns <- setdiff(required_columns, colnames(seurat_obj@meta.data))
  
  if (length(missing_columns) > 0) {
    stop(
      "Missing required QC metadata columns for sample '", sample_name, "': ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  cells_before <- ncol(seurat_obj)
  
  keep_cells <-
    seurat_obj$percent.mt < qc_thresholds$percent_mt_max &
    seurat_obj$percent.ribo < qc_thresholds$percent_ribo_max &
    seurat_obj$percent.hb < qc_thresholds$percent_hb_max &
    seurat_obj$nFeature_RNA > qc_thresholds$nFeature_RNA_min &
    seurat_obj$nFeature_RNA < qc_thresholds$nFeature_RNA_max &
    seurat_obj$nCount_RNA > qc_thresholds$nCount_RNA_min &
    seurat_obj$nCount_RNA < qc_thresholds$nCount_RNA_max
  
  keep_cells[is.na(keep_cells)] <- FALSE
  
  if (sum(keep_cells) == 0) {
    stop("No cells passed QC for sample '", sample_name, "'.")
  }
  
  seurat_obj <- seurat_obj[, keep_cells]
  
  filtering_summary <- tibble(
    sample = sample_name,
    cells_before_qc = cells_before,
    cells_after_qc = ncol(seurat_obj),
    cells_removed_qc = cells_before - ncol(seurat_obj),
    percent_removed_qc = round((cells_before - ncol(seurat_obj)) / cells_before * 100, 2)
  )
  
  list(
    seurat_obj = seurat_obj,
    filtering_summary = filtering_summary
  )
}


remove_immune_receptor_genes <- function(seurat_obj,
                                         gene_patterns,
                                         sample_name,
                                         assay = "RNA") {
  DefaultAssay(seurat_obj) <- assay
  
  counts <- get_counts_matrix(seurat_obj, assay = assay)
  
  receptor_genes <- unique(unlist(
    lapply(gene_patterns, function(pattern) {
      grep(pattern, rownames(counts), value = TRUE)
    })
  ))
  
  if (length(receptor_genes) == 0) {
    message("No immune receptor genes matched for sample: ", sample_name)
    seurat_obj$nCount_RNA_no_receptor <- seurat_obj$nCount_RNA
    seurat_obj$nFeature_RNA_no_receptor <- seurat_obj$nFeature_RNA
    
    return(list(
      seurat_obj = seurat_obj,
      n_removed_genes = 0
    ))
  }
  
  non_receptor_genes <- setdiff(rownames(counts), receptor_genes)
  
  if (length(non_receptor_genes) == 0) {
    stop("All genes matched immune receptor filters for sample '", sample_name, "'.")
  }
  
  filtered_counts <- counts[non_receptor_genes, , drop = FALSE]
  
  # Preserve original RNA QC metrics before replacing the assay.
  if ("nCount_RNA" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$nCount_RNA_pre_receptor_filter <- seurat_obj$nCount_RNA
  }
  
  if ("nFeature_RNA" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$nFeature_RNA_pre_receptor_filter <- seurat_obj$nFeature_RNA
  }
  
  seurat_obj[[assay]] <- create_assay_from_counts(filtered_counts)
  DefaultAssay(seurat_obj) <- assay
  
  seurat_obj$nCount_RNA_no_receptor <- Matrix::colSums(filtered_counts)
  seurat_obj$nFeature_RNA_no_receptor <- Matrix::colSums(filtered_counts > 0)
  
  list(
    seurat_obj = seurat_obj,
    n_removed_genes = length(receptor_genes)
  )
}


run_sctransform_and_pca <- function(seurat_obj,
                                    sct_config,
                                    run_pca = TRUE,
                                    max_pcs = 50,
                                    sample_name,
                                    assay = "RNA") {
  DefaultAssay(seurat_obj) <- assay
  
  vars_to_regress <- sct_config$vars_to_regress %||% NULL
  
  if (!is.null(vars_to_regress)) {
    vars_to_regress <- vars_to_regress[
      vars_to_regress %in% colnames(seurat_obj@meta.data)
    ]
    
    if (length(vars_to_regress) == 0) {
      vars_to_regress <- NULL
    }
  }
  
  sct_args <- list(
    object = seurat_obj,
    assay = assay,
    variable.features.n = sct_config$variable_features_n %||% 2000,
    ncells = min(sct_config$ncells %||% 3000, ncol(seurat_obj)),
    return.only.var.genes = sct_config$return_only_var_genes %||% FALSE,
    seed.use = sct_config$seed_use %||% 1234,
    verbose = sct_config$verbose %||% FALSE
  )
  
  if (!is.null(vars_to_regress)) {
    sct_args$vars.to.regress <- vars_to_regress
  }
  
  seurat_obj <- do.call(SCTransform, sct_args)
  
  variable_features <- VariableFeatures(seurat_obj)
  
  if (run_pca && length(variable_features) > 0) {
    npcs <- min(max_pcs, length(variable_features), ncol(seurat_obj) - 1)
    
    if (npcs >= 2) {
      seurat_obj <- RunPCA(
        seurat_obj,
        features = variable_features,
        npcs = npcs,
        verbose = FALSE
      )
    } else {
      warning("Skipping PCA for sample '", sample_name, "' because fewer than 2 PCs are possible.")
    }
  }
  
  seurat_obj
}


save_variable_feature_plot <- function(seurat_obj,
                                       sample_name,
                                       output_dir,
                                       plot_suffix = "") {
  create_dir_if_missing(output_dir)
  
  variable_features <- VariableFeatures(seurat_obj)
  
  if (length(variable_features) == 0) {
    warning("No variable features available for sample '", sample_name, "'.")
    return(invisible(NULL))
  }
  
  top10 <- head(variable_features, 10)
  
  plot1 <- VariableFeaturePlot(seurat_obj)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  
  output_file <- file.path(
    output_dir,
    paste0(safe_filename(sample_name), plot_suffix, "_variable_features.png")
  )
  
  ggsave(
    filename = output_file,
    plot = plot2,
    width = 10,
    height = 10,
    dpi = 300
  )
  
  invisible(output_file)
}


save_qc_violin_plot <- function(seurat_obj,
                                sample_name,
                                output_dir,
                                plot_suffix = "",
                                group_by = "sample") {
  create_dir_if_missing(output_dir)
  
  features_to_plot <- c(
    "percent.mt", "percent.ribo", "percent.hb",
    "nCount_RNA", "nFeature_RNA"
  )
  
  features_to_plot <- features_to_plot[
    features_to_plot %in% colnames(seurat_obj@meta.data)
  ]
  
  if (length(features_to_plot) == 0) {
    warning("No QC features available to plot for sample '", sample_name, "'.")
    return(invisible(NULL))
  }
  
  if (!group_by %in% colnames(seurat_obj@meta.data)) {
    group_by <- NULL
  }
  
  p_violin <- VlnPlot(
    seurat_obj,
    features = features_to_plot,
    group.by = group_by,
    ncol = 2,
    pt.size = 0
  )
  
  output_file <- file.path(
    output_dir,
    paste0(safe_filename(sample_name), plot_suffix, "_qc_violin.png")
  )
  
  ggsave(
    filename = output_file,
    plot = p_violin,
    width = 15,
    height = 10,
    dpi = 300
  )
  
  invisible(output_file)
}


make_post_qc_summary <- function(seurat_obj,
                                 sample_name,
                                 filtering_summary,
                                 n_receptor_genes_removed) {
  tibble(
    sample = sample_name,
    n_cells = ncol(seurat_obj),
    median_nCount_RNA = median(seurat_obj$nCount_RNA, na.rm = TRUE),
    median_nFeature_RNA = median(seurat_obj$nFeature_RNA, na.rm = TRUE),
    median_log10_UMI = median(log10(seurat_obj$nCount_RNA + 1), na.rm = TRUE),
    median_percent_mt = median(seurat_obj$percent.mt, na.rm = TRUE),
    median_percent_ribo = median(seurat_obj$percent.ribo, na.rm = TRUE),
    median_percent_hb = median(seurat_obj$percent.hb, na.rm = TRUE),
    n_variable_features = length(VariableFeatures(seurat_obj)),
    n_receptor_genes_removed = n_receptor_genes_removed
  ) %>%
    bind_cols(filtering_summary %>% select(-sample))
}


standardise_input_object <- function(input_object, dataset_id) {
  if (inherits(input_object, "Seurat")) {
    output <- list()
    output[[dataset_id]] <- input_object
    return(output)
  }
  
  if (!is.list(input_object)) {
    stop("Input RDS must contain either a Seurat object or a list of Seurat objects.")
  }
  
  if (is.null(names(input_object)) || any(names(input_object) == "")) {
    names(input_object) <- paste0("sample_", seq_along(input_object))
  }
  
  input_object
}


process_one_sample <- function(seurat_obj, sample_name, config) {
  message("\nProcessing sample: ", sample_name)
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Object is not a Seurat object.")
  }
  
  assay <- "RNA"
  
  if (!assay %in% Assays(seurat_obj)) {
    stop("RNA assay not found.")
  }
  
  DefaultAssay(seurat_obj) <- assay
  
  seurat_obj <- add_qc_metrics(
    seurat_obj = seurat_obj,
    qc_gene_patterns = config$qc_gene_patterns,
    assay = assay
  )
  
  qc_result <- qc_filter_object(
    seurat_obj = seurat_obj,
    qc_thresholds = config$qc_thresholds,
    sample_name = sample_name
  )
  
  seurat_obj <- qc_result$seurat_obj
  
  n_receptor_genes_removed <- 0
  
  if (isTRUE(config$remove_immune_receptor_genes)) {
    receptor_result <- remove_immune_receptor_genes(
      seurat_obj = seurat_obj,
      gene_patterns = config$immune_receptor_gene_patterns,
      sample_name = sample_name,
      assay = assay
    )
    
    seurat_obj <- receptor_result$seurat_obj
    n_receptor_genes_removed <- receptor_result$n_removed_genes
  }
  
  seurat_obj <- run_sctransform_and_pca(
    seurat_obj = seurat_obj,
    sct_config = config$sctransform,
    run_pca = isTRUE(config$run_pca),
    max_pcs = config$max_pcs %||% 50,
    sample_name = sample_name,
    assay = assay
  )
  
  save_variable_feature_plot(
    seurat_obj = seurat_obj,
    sample_name = sample_name,
    output_dir = config$variable_feature_plot_dir,
    plot_suffix = config$plot_suffix %||% ""
  )
  
  save_qc_violin_plot(
    seurat_obj = seurat_obj,
    sample_name = sample_name,
    output_dir = config$qc_plot_dir,
    plot_suffix = config$plot_suffix %||% "",
    group_by = "sample"
  )
  
  summary_row <- make_post_qc_summary(
    seurat_obj = seurat_obj,
    sample_name = sample_name,
    filtering_summary = qc_result$filtering_summary,
    n_receptor_genes_removed = n_receptor_genes_removed
  )
  
  list(
    seurat_obj = seurat_obj,
    summary_row = summary_row
  )
}


process_dataset <- function(config) {
  dataset_id <- config$dataset_id %||% "dataset"
  
  message("\n============================================================")
  message("Starting dataset: ", dataset_id)
  message("============================================================")
  
  output_dirs <- c(
    dirname(config$output_rds),
    dirname(config$post_qc_summary_csv),
    dirname(config$failed_samples_csv),
    config$qc_plot_dir,
    config$variable_feature_plot_dir
  )
  
  invisible(lapply(unique(output_dirs), create_dir_if_missing))
  
  input_object <- readRDS(config$input_rds)
  seurat_objects <- standardise_input_object(input_object, dataset_id)
  
  processed_objects <- list()
  summary_rows <- list()
  failed_samples <- list()
  
  for (sample_name in names(seurat_objects)) {
    result <- tryCatch(
      process_one_sample(
        seurat_obj = seurat_objects[[sample_name]],
        sample_name = sample_name,
        config = config
      ),
      error = function(e) {
        warning("Sample failed: ", sample_name, " | ", e$message)
        
        failed_samples[[sample_name]] <<- tibble(
          dataset_id = dataset_id,
          sample = sample_name,
          error_message = e$message
        )
        
        NULL
      }
    )
    
    if (!is.null(result)) {
      processed_objects[[sample_name]] <- result$seurat_obj
      summary_rows[[sample_name]] <- result$summary_row
    }
  }
  
  if (length(processed_objects) == 0) {
    stop("No samples were processed successfully for dataset: ", dataset_id)
  }
  
  saveRDS(processed_objects, file = config$output_rds)
  
  summary_df <- bind_rows(summary_rows) %>%
    mutate(dataset_id = dataset_id, .before = sample)
  
  write.csv(summary_df, config$post_qc_summary_csv, row.names = FALSE)
  
  if (length(failed_samples) > 0) {
    failed_df <- bind_rows(failed_samples)
    write.csv(failed_df, config$failed_samples_csv, row.names = FALSE)
  }
  
  message("\nFinished dataset: ", dataset_id)
  message("Processed samples: ", length(processed_objects))
  message("Failed samples: ", length(failed_samples))
  message("Saved RDS: ", config$output_rds)
  message("Saved summary: ", config$post_qc_summary_csv)
  
  invisible(list(
    processed_objects = processed_objects,
    summary = summary_df,
    failed_samples = failed_samples
  ))
}


# -----------------------------#
# 3. Run pipeline
# -----------------------------#

results <- lapply(datasets, process_dataset)

message("\nAll configured datasets processed.")
