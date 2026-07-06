#!/usr/bin/env Rscript

# ==============================================================================
# Create Seurat objects and QC plots from Cell Ranger sample outputs
# ==============================================================================
#
# Purpose:
#   This script searches one or more input directories for filtered 10x/Cell Ranger
#   count matrices, creates one Seurat object per sample, calculates basic QC
#   metrics, saves a list of Seurat objects, and writes QC summary outputs.
#
# Notes for repository use:
#   - Keep project-specific paths in the CONFIG section only.
#   - Do not commit large RDS files or raw count matrices to Git.
#   - This script assumes each sample has a Cell Ranger-style matrix folder,
#     for example: sample_filtered_feature_bc_matrix/
#
# ==============================================================================


# ----------------------------- #
# 0. Load packages
# ----------------------------- #

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(stringr)
  library(Matrix)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  # One or more directories to search for sample matrix folders.
  # Replace these with your own paths.
  input_roots = c(
    "path/to/cellranger_outputs"
  ),
  
  # Directory where QC plots and tables will be written.
  output_dir = "path/to/output/qc",
  
  # File path for the saved list of Seurat objects.
  output_rds = "path/to/output/seurat_objects_qc.rds",
  
  # File name for the per-sample QC summary.
  qc_summary_file = "sample_qc_summary.csv",
  
  # File name for any samples that fail during processing.
  failed_samples_file = "failed_samples.csv",
  
  # Pattern used to identify sample matrix directories.
  # Edit this if your folders are named differently.
  matrix_dir_pattern = "sample_filtered_feature_bc_matrix$",
  
  # Optional pattern used to exclude known empty/failed samples.
  # Set to NULL if no samples need to be excluded.
  exclude_pattern = NULL,
  
  # How many samples to include per QC violin plot.
  samples_per_plot = 16,
  
  # Which parent directory should be used as the sample name?
  # 1 = parent of the matrix folder
  # 2 = grandparent of the matrix folder
  # This matches basename(dirname(dirname(matrix_path))) from many Cell Ranger
  # multi outputs.
  sample_name_parent_level = 2,
  
  # Assay names expected from Read10X() when multiple assays are present.
  gene_expression_assay_name = "Gene Expression",
  antibody_capture_assay_name = "Antibody Capture",
  
  # RNA filtering during Seurat object creation.
  min_cells = 3,
  min_features = 0,
  
  # Gene-name patterns used for QC metrics.
  # These defaults support common human and mouse gene symbols.
  mitochondrial_pattern = "^MT-|^mt-",
  ribosomal_pattern = "^RPS|^RPL|^Rps|^Rpl",
  hemoglobin_pattern = "^HB[AB]|^HBA|^HBB|^Hba|^Hbb"
)


# ----------------------------- #
# 2. Helper functions
# ----------------------------- #

make_output_dirs <- function(config) {
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(config$output_rds), recursive = TRUE, showWarnings = FALSE)
}


find_matrix_dirs <- function(input_roots, matrix_dir_pattern, exclude_pattern = NULL) {
  matrix_dirs <- input_roots |>
    map(~ list.dirs(.x, recursive = TRUE, full.names = TRUE)) |>
    unlist(use.names = FALSE) |>
    str_subset(matrix_dir_pattern)
  
  if (!is.null(exclude_pattern)) {
    matrix_dirs <- matrix_dirs[!str_detect(matrix_dirs, exclude_pattern)]
  }
  
  unique(matrix_dirs)
}


get_parent_name <- function(path, parent_level = 2) {
  sample_path <- path
  
  for (i in seq_len(parent_level)) {
    sample_path <- dirname(sample_path)
  }
  
  basename(sample_path)
}


get_counts_layer <- function(seurat_obj, assay = "RNA") {
  # Seurat v5 uses layer = "counts"; older versions used slot = "counts".
  tryCatch(
    GetAssayData(seurat_obj, assay = assay, layer = "counts"),
    error = function(e) GetAssayData(seurat_obj, assay = assay, slot = "counts")
  )
}


create_assay_from_counts <- function(counts, min_cells = 0) {
  if ("CreateAssay5Object" %in% getNamespaceExports("Seurat")) {
    CreateAssay5Object(counts = counts, min.cells = min_cells)
  } else {
    CreateAssayObject(counts = counts, min.cells = min_cells)
  }
}


extract_10x_assays <- function(
    data,
    gene_expression_assay_name = "Gene Expression",
    antibody_capture_assay_name = "Antibody Capture"
) {
  if (is.list(data)) {
    if (!gene_expression_assay_name %in% names(data)) {
      stop(
        "Could not find gene expression assay named '",
        gene_expression_assay_name,
        "'. Available assays: ",
        paste(names(data), collapse = ", ")
      )
    }
    
    gene_counts <- data[[gene_expression_assay_name]]
    
    adt_counts <- NULL
    if (antibody_capture_assay_name %in% names(data)) {
      adt_counts <- data[[antibody_capture_assay_name]]
    }
  } else {
    gene_counts <- data
    adt_counts <- NULL
  }
  
  list(
    gene_counts = gene_counts,
    adt_counts = adt_counts
  )
}


add_percent_feature_metric <- function(seurat_obj, gene_pattern, metric_name) {
  counts <- get_counts_layer(seurat_obj, assay = "RNA")
  genes <- grep(gene_pattern, rownames(counts), value = TRUE)
  
  if (length(genes) == 0) {
    seurat_obj[[metric_name]] <- 0
  } else {
    total_counts <- seurat_obj$nCount_RNA
    total_counts[total_counts == 0] <- NA_real_
    
    seurat_obj[[metric_name]] <-
      Matrix::colSums(counts[genes, , drop = FALSE]) / total_counts * 100
    
    seurat_obj[[metric_name]][is.na(seurat_obj[[metric_name]])] <- 0
  }
  
  seurat_obj
}


process_sample <- function(data_dir, sample_name, config) {
  message("Processing: ", sample_name)
  
  raw_data <- Read10X(data.dir = data_dir)
  
  assays <- extract_10x_assays(
    data = raw_data,
    gene_expression_assay_name = config$gene_expression_assay_name,
    antibody_capture_assay_name = config$antibody_capture_assay_name
  )
  
  seurat_obj <- CreateSeuratObject(
    counts = assays$gene_counts,
    project = sample_name,
    min.cells = config$min_cells,
    min.features = config$min_features
  )
  
  seurat_obj$sample_id <- sample_name
  
  if (!is.null(assays$adt_counts)) {
    seurat_obj[["ADT"]] <- create_assay_from_counts(
      counts = assays$adt_counts,
      min_cells = 0
    )
  }
  
  seurat_obj$log10_UMI <- log10(seurat_obj$nCount_RNA + 1)
  
  seurat_obj <- add_percent_feature_metric(
    seurat_obj,
    gene_pattern = config$mitochondrial_pattern,
    metric_name = "percent.mt"
  )
  
  seurat_obj <- add_percent_feature_metric(
    seurat_obj,
    gene_pattern = config$ribosomal_pattern,
    metric_name = "percent.ribo"
  )
  
  seurat_obj <- add_percent_feature_metric(
    seurat_obj,
    gene_pattern = config$hemoglobin_pattern,
    metric_name = "percent.hb"
  )
  
  seurat_obj
}


make_cell_qc_metadata <- function(seurat_obj, sample_name) {
  seurat_obj@meta.data |>
    select(
      nCount_RNA,
      nFeature_RNA,
      log10_UMI,
      percent.mt,
      percent.ribo,
      percent.hb
    ) |>
    mutate(sample_id = sample_name)
}


make_sample_qc_summary <- function(seurat_obj, sample_name) {
  tibble(
    sample_id = sample_name,
    n_cells = ncol(seurat_obj),
    median_nCount_RNA = median(seurat_obj$nCount_RNA, na.rm = TRUE),
    median_nFeature_RNA = median(seurat_obj$nFeature_RNA, na.rm = TRUE),
    median_log10_UMI = median(seurat_obj$log10_UMI, na.rm = TRUE),
    median_percent_mt = median(seurat_obj$percent.mt, na.rm = TRUE),
    median_percent_ribo = median(seurat_obj$percent.ribo, na.rm = TRUE),
    median_percent_hb = median(seurat_obj$percent.hb, na.rm = TRUE)
  )
}


plot_qc_violins <- function(qc_metadata, output_dir, samples_per_plot = 16) {
  if (nrow(qc_metadata) == 0) {
    warning("No QC metadata available for plotting.")
    return(invisible(NULL))
  }
  
  metrics <- c(
    "nCount_RNA",
    "nFeature_RNA",
    "log10_UMI",
    "percent.mt",
    "percent.ribo",
    "percent.hb"
  )
  
  long_df <- qc_metadata |>
    mutate(sample_id = factor(sample_id, levels = unique(sample_id))) |>
    pivot_longer(
      cols = all_of(metrics),
      names_to = "metric",
      values_to = "value"
    )
  
  sample_groups <- split(
    levels(long_df$sample_id),
    ceiling(seq_along(levels(long_df$sample_id)) / samples_per_plot)
  )
  
  for (i in seq_along(sample_groups)) {
    group_samples <- sample_groups[[i]]
    
    df_subset <- long_df |>
      filter(sample_id %in% group_samples)
    
    p <- ggplot(df_subset, aes(x = sample_id, y = value, fill = sample_id)) +
      geom_violin(scale = "width", trim = TRUE) +
      facet_wrap(~ metric, scales = "free_y", ncol = 3) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)
      ) +
      labs(
        title = paste0("QC violin plots: sample group ", i),
        x = NULL,
        y = NULL
      )
    
    ggsave(
      filename = file.path(output_dir, paste0("QC_violin_samples_", i, ".png")),
      plot = p,
      width = 16,
      height = 10,
      dpi = 300
    )
  }
  
  invisible(NULL)
}


# ----------------------------- #
# 3. Run pipeline
# ----------------------------- #

make_output_dirs(config)

sample_paths <- find_matrix_dirs(
  input_roots = config$input_roots,
  matrix_dir_pattern = config$matrix_dir_pattern,
  exclude_pattern = config$exclude_pattern
)

if (length(sample_paths) == 0) {
  stop(
    "No matrix folders found. Check config$input_roots and ",
    "config$matrix_dir_pattern."
  )
}

sample_table <- tibble(
  matrix_dir = sample_paths,
  sample_id = map_chr(
    sample_paths,
    get_parent_name,
    parent_level = config$sample_name_parent_level
  )
) |>
  mutate(sample_id = make.unique(sample_id, sep = "_"))

message("Found ", nrow(sample_table), " sample matrix folder(s).")

seurat_objects <- list()
qc_metadata_list <- list()
qc_summary_list <- list()
failed_samples <- list()

for (i in seq_len(nrow(sample_table))) {
  sample_id <- sample_table$sample_id[[i]]
  matrix_dir <- sample_table$matrix_dir[[i]]
  
  result <- tryCatch(
    {
      obj <- process_sample(
        data_dir = matrix_dir,
        sample_name = sample_id,
        config = config
      )
      
      seurat_objects[[sample_id]] <- obj
      qc_metadata_list[[sample_id]] <- make_cell_qc_metadata(obj, sample_id)
      qc_summary_list[[sample_id]] <- make_sample_qc_summary(obj, sample_id)
      
      TRUE
    },
    error = function(e) {
      warning("Failed to process ", sample_id, ": ", conditionMessage(e))
      
      failed_samples[[sample_id]] <<- tibble(
        sample_id = sample_id,
        matrix_dir = matrix_dir,
        error = conditionMessage(e)
      )
      
      FALSE
    }
  )
}

qc_summary_df <- bind_rows(qc_summary_list)
qc_metadata_df <- bind_rows(qc_metadata_list)
failed_samples_df <- bind_rows(failed_samples)

write.csv(
  qc_summary_df,
  file.path(config$output_dir, config$qc_summary_file),
  row.names = FALSE
)

if (nrow(failed_samples_df) > 0) {
  write.csv(
    failed_samples_df,
    file.path(config$output_dir, config$failed_samples_file),
    row.names = FALSE
  )
}

plot_qc_violins(
  qc_metadata = qc_metadata_df,
  output_dir = config$output_dir,
  samples_per_plot = config$samples_per_plot
)

saveRDS(seurat_objects, config$output_rds)

message("Processing complete.")
message("Successfully processed samples: ", length(seurat_objects))
message("Failed samples: ", nrow(failed_samples_df))
message("QC summary: ", file.path(config$output_dir, config$qc_summary_file))
message("Seurat objects: ", config$output_rds)

sessionInfo()
