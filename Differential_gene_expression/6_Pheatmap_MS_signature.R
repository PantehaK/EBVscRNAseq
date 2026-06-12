#!/usr/bin/env Rscript

# ==============================================================================
# Heatmap of stimulated latent DE genes by sample
# ==============================================================================
#
# Purpose:
#   This script creates an unsupervised sample-by-gene heatmap using genes from
#   pseudobulk differential expression analysis of stimulated latent EBV-specific
#   T cells in MS versus Control.
#
#   It subsets to:
#     - latent_group == "Stimulated Latent"
#     - matched_cluster == target_cluster
#
#   It then averages expression per sample, z-scores each gene, transposes the
#   matrix so samples are rows and genes are columns, and plots a pheatmap with
#   cohort annotation.
#
# Expected input:
#   Activated EBV clustered/module-scored object, e.g.
#     15_activated_EBV_MS_stimLat_signature_scored.rds
#
# Optional input:
#   edgeR result CSV from stimulated latent pseudobulk analysis.
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
  library(stringr)
  library(readr)
  library(pheatmap)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  dataset_id = "activated_ebv_stimulated_latent_heatmap",
  
  input_rds = "path/to/input/15_activated_EBV_MS_stimLat_signature_scored.rds",
  
  output_dir = "path/to/output/Publication_data/EBV_activated/StimLat_heatmap",
  
  heatmap_pdf = "path/to/output/Publication_data/EBV_activated/StimLat_heatmap/Cluster2_StimulatedLatent_DEG_heatmap.pdf",
  heatmap_png = "path/to/output/Publication_data/EBV_activated/StimLat_heatmap/Cluster2_StimulatedLatent_DEG_heatmap.png",
  averaged_expression_csv = "path/to/output/Publication_data/EBV_activated/StimLat_heatmap/Cluster2_StimulatedLatent_average_expression_by_sample.csv",
  zscore_matrix_csv = "path/to/output/Publication_data/EBV_activated/StimLat_heatmap/Cluster2_StimulatedLatent_zscore_matrix_samples_by_genes.csv",
  gene_presence_csv = "path/to/output/Publication_data/EBV_activated/StimLat_heatmap/Cluster2_StimulatedLatent_gene_presence.csv",
  sample_annotation_csv = "path/to/output/Publication_data/EBV_activated/StimLat_heatmap/Cluster2_StimulatedLatent_sample_annotation.csv",
  session_info_file = "path/to/output/Publication_data/EBV_activated/StimLat_heatmap/sessionInfo_stimLat_heatmap.txt",
  
  # Metadata columns.
  batch_col = "batch",
  lifecycle_col = "lifecycle",
  sample_col = "sample",
  cohort_col = "cohort",
  cluster_col = "matched_cluster",
  
  # Latent group labels.
  latent_group_col = "latent_group",
  latent_label = "Latent",
  non_latent_label = "Non-Latent",
  baseline_latent_label = "Baseline Latent",
  stimulated_latent_label = "Stimulated Latent",
  baseline_batch_pattern = "^GEMEBV",
  stimulated_batch_pattern = "^EBVLCL",
  
  # Heatmap subset.
  target_cluster = "2",
  
  # Assay/expression settings.
  assay = "SCT",
  expression_slot = "scale.data",
  
  # Optional edgeR gene input.
  use_edgeR_de_csv = TRUE,
  edgeR_de_csv = "path/to/input/Age_adjusted_EdgeR_StimulatedLatent_MS_vs_Control_by_matched_cluster_ALL.csv",
  edgeR_cluster_col = "matched_cluster",
  edgeR_gene_col = "gene",
  edgeR_fdr_col = "p_val_adj",
  edgeR_logfc_col = "avg_log2FC",
  edgeR_fdr_threshold = 0.05,
  edgeR_abs_logfc_threshold = 0,
  
  # Fallback gene list, used if edgeR CSV is unavailable.
  fallback_gene_list = c(
    "FTH1", "JAK1", "ITM2B", "RPL13", "SMCHD1", "IL7R", "HLA-B", "PTMA",
    "HERC5", "TOMM7", "MT-CYB", "H3-3A", "HLA-F", "RPS10", "MT-ND4L",
    "HNRNPA1", "FYN", "PFN1", "RPS29", "MT-CO3", "RPS15", "RPS18",
    "SLFN5", "ZBTB20", "RPLP2", "MT-ND3", "ARPC2", "IFI6", "CD3D",
    "RPS27", "HLA-C", "CXCR4", "RPS27A", "CALM1", "CCSER2", "ANKRD12",
    "VIM", "RPS9", "CDC42SE2", "EEF1G", "RPS6", "PDE3B", "EEF1A1",
    "RPS16", "YWHAB", "RPL30", "SARAF", "SRSF7"
  ),
  
  heatmap_colours = c(
    "#011f4b", "#03396c", "#005b96", "#6497b1", "#b3cde0",
    "#ffefea", "#fbd9d3", "#ffb09c", "#fe5757", "#cb2424", "#900000"
  ),
  
  cohort_colours = c(
    "Control" = "#6CC3F4",
    "MS" = "#E65757"
  ),
  
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 9,
  fontsize_col = 10,
  pdf_width = 15,
  pdf_height = 8,
  png_width = 15,
  png_height = 8,
  png_dpi = 400,
  
  verbose = TRUE
)


# ----------------------------- #
# 2. Helper functions
# ----------------------------- #

create_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}


create_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}


setup_dirs <- function(cfg) {
  create_dir(cfg$output_dir)
  create_parent_dir(cfg$heatmap_pdf)
  create_parent_dir(cfg$heatmap_png)
  create_parent_dir(cfg$averaged_expression_csv)
  create_parent_dir(cfg$zscore_matrix_csv)
  create_parent_dir(cfg$gene_presence_csv)
  create_parent_dir(cfg$sample_annotation_csv)
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


check_required_metadata <- function(obj, cfg) {
  required_cols <- c(
    cfg$batch_col,
    cfg$lifecycle_col,
    cfg$sample_col,
    cfg$cohort_col,
    cfg$cluster_col
  )
  
  missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
  
  if (length(missing_cols) > 0) {
    stop(
      "Object metadata is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
}


add_latent_group_if_needed <- function(obj, cfg) {
  if (cfg$latent_group_col %in% colnames(obj@meta.data)) {
    return(obj)
  }
  
  md <- obj@meta.data
  
  obj[[cfg$latent_group_col]] <- case_when(
    !is.na(md[[cfg$batch_col]]) &
      grepl(cfg$baseline_batch_pattern, md[[cfg$batch_col]]) &
      md[[cfg$lifecycle_col]] == cfg$latent_label ~ cfg$baseline_latent_label,
    
    !is.na(md[[cfg$batch_col]]) &
      grepl(cfg$stimulated_batch_pattern, md[[cfg$batch_col]]) &
      md[[cfg$lifecycle_col]] == cfg$latent_label ~ cfg$stimulated_latent_label,
    
    TRUE ~ cfg$non_latent_label
  )
  
  obj[[cfg$latent_group_col]] <- factor(
    obj[[cfg$latent_group_col, drop = TRUE]],
    levels = c(
      cfg$non_latent_label,
      cfg$baseline_latent_label,
      cfg$stimulated_latent_label
    )
  )
  
  obj
}


subset_heatmap_cells <- function(obj, cfg) {
  md <- obj@meta.data
  
  keep_cells <- rownames(md)[
    md[[cfg$latent_group_col]] == cfg$stimulated_latent_label &
      as.character(md[[cfg$cluster_col]]) == as.character(cfg$target_cluster)
  ]
  
  if (length(keep_cells) == 0) {
    stop(
      "No cells found for ",
      cfg$stimulated_latent_label,
      " and ",
      cfg$cluster_col,
      " == ",
      cfg$target_cluster
    )
  }
  
  subset(obj, cells = keep_cells)
}


read_gene_list <- function(cfg) {
  if (isTRUE(cfg$use_edgeR_de_csv) && file.exists(cfg$edgeR_de_csv)) {
    de <- read_csv(cfg$edgeR_de_csv, show_col_types = FALSE)
    
    required_cols <- c(
      cfg$edgeR_gene_col,
      cfg$edgeR_fdr_col,
      cfg$edgeR_logfc_col
    )
    
    missing_cols <- setdiff(required_cols, colnames(de))
    
    if (length(missing_cols) == 0) {
      if (cfg$edgeR_cluster_col %in% colnames(de)) {
        de <- de |>
          filter(as.character(.data[[cfg$edgeR_cluster_col]]) == as.character(cfg$target_cluster))
      }
      
      genes <- de |>
        filter(
          !is.na(.data[[cfg$edgeR_fdr_col]]),
          .data[[cfg$edgeR_fdr_col]] < cfg$edgeR_fdr_threshold,
          !is.na(.data[[cfg$edgeR_logfc_col]]),
          abs(.data[[cfg$edgeR_logfc_col]]) > cfg$edgeR_abs_logfc_threshold
        ) |>
        arrange(.data[[cfg$edgeR_fdr_col]]) |>
        pull(.data[[cfg$edgeR_gene_col]]) |>
        unique()
      
      if (length(genes) >= 2) {
        return(genes)
      }
      
      warning("edgeR CSV found but fewer than two genes passed thresholds. Using fallback gene list.")
    } else {
      warning(
        "edgeR CSV is missing required column(s): ",
        paste(missing_cols, collapse = ", "),
        ". Using fallback gene list."
      )
    }
  } else {
    warning("edgeR CSV not found. Using fallback gene list.")
  }
  
  unique(cfg$fallback_gene_list)
}


get_genes_present <- function(obj, genes, cfg) {
  if (!cfg$assay %in% Assays(obj)) {
    stop("Assay not found in object: ", cfg$assay)
  }
  
  genes_present <- intersect(genes, rownames(obj[[cfg$assay]]))
  genes_missing <- setdiff(genes, genes_present)
  
  gene_presence <- tibble(
    gene = genes,
    present = genes %in% genes_present
  )
  
  write.csv(
    gene_presence,
    cfg$gene_presence_csv,
    row.names = FALSE
  )
  
  if (length(genes_present) < 2) {
    stop("Too few genes from gene list were found in the object after subsetting.")
  }
  
  if (length(genes_missing) > 0) {
    message(
      "Missing genes: ",
      paste(genes_missing, collapse = ", ")
    )
  }
  
  genes_present
}


make_average_expression_matrix <- function(obj, genes_present, cfg) {
  DefaultAssay(obj) <- cfg$assay
  
  agg <- AverageExpression(
    object = obj,
    features = genes_present,
    group.by = cfg$sample_col,
    assays = cfg$assay,
    slot = cfg$expression_slot,
    return.seurat = FALSE
  )
  
  mat <- agg[[cfg$assay]]
  
  if (is.null(mat) || nrow(mat) < 2 || ncol(mat) < 2) {
    stop("Average expression matrix is too small for heatmap clustering.")
  }
  
  mat
}


zscore_rows <- function(mat) {
  mat_z <- t(scale(t(as.matrix(mat))))
  mat_z[is.na(mat_z)] <- 0
  mat_z
}


make_sample_annotation <- function(obj, sample_names, cfg) {
  meta_samp <- obj@meta.data |>
    select(
      sample = .data[[cfg$sample_col]],
      cohort = .data[[cfg$cohort_col]]
    ) |>
    distinct(sample, .keep_all = TRUE)
  
  anno_row <- data.frame(
    cohort = meta_samp$cohort
  )
  
  rownames(anno_row) <- meta_samp$sample
  
  anno_row <- anno_row[sample_names, , drop = FALSE]
  
  if (any(is.na(anno_row$cohort))) {
    warning("Some samples are missing cohort annotation.")
  }
  
  anno_row
}


save_heatmap <- function(mat_z_t, anno_row, cfg) {
  anno_colors <- list(
    cohort = cfg$cohort_colours
  )
  
  heatmap_cols <- colorRampPalette(cfg$heatmap_colours)(200)
  
  pheatmap(
    mat_z_t,
    cluster_rows = cfg$cluster_rows,
    cluster_cols = cfg$cluster_cols,
    annotation_row = anno_row,
    annotation_colors = anno_colors,
    color = heatmap_cols,
    fontsize_row = cfg$fontsize_row,
    fontsize_col = cfg$fontsize_col,
    border_color = NA,
    filename = cfg$heatmap_pdf,
    width = cfg$pdf_width,
    height = cfg$pdf_height
  )
  
  pheatmap(
    mat_z_t,
    cluster_rows = cfg$cluster_rows,
    cluster_cols = cfg$cluster_cols,
    annotation_row = anno_row,
    annotation_colors = anno_colors,
    color = heatmap_cols,
    fontsize_row = cfg$fontsize_row,
    fontsize_col = cfg$fontsize_col,
    border_color = NA,
    filename = cfg$heatmap_png,
    width = cfg$png_width,
    height = cfg$png_height
  )
}


# ----------------------------- #
# 3. Run pipeline
# ----------------------------- #

setup_dirs(config)

message("Loading activated EBV object...")
obj <- load_seurat_object(config$input_rds)

message("Checking required metadata...")
check_required_metadata(obj, config)

message("Defining latent_group if needed...")
obj <- add_latent_group_if_needed(obj, config)

message("Subsetting stimulated latent cluster ", config$target_cluster, "...")
obj_sub <- subset_heatmap_cells(obj, config)

message("Reading heatmap gene list...")
genes <- read_gene_list(config)

message("Checking genes present in object...")
genes_present <- get_genes_present(obj_sub, genes, config)

message("Calculating average expression per sample...")
mat <- make_average_expression_matrix(obj_sub, genes_present, config)

write.csv(
  as.data.frame(mat) |>
    rownames_to_column("gene"),
  config$averaged_expression_csv,
  row.names = FALSE
)

message("Z-scoring genes and transposing matrix...")
mat_z <- zscore_rows(mat)
mat_z_t <- t(mat_z)

write.csv(
  as.data.frame(mat_z_t) |>
    rownames_to_column("sample"),
  config$zscore_matrix_csv,
  row.names = FALSE
)

message("Creating sample annotation...")
anno_row <- make_sample_annotation(
  obj = obj_sub,
  sample_names = rownames(mat_z_t),
  cfg = config
)

write.csv(
  anno_row |>
    rownames_to_column("sample"),
  config$sample_annotation_csv,
  row.names = FALSE
)

message("Saving heatmap...")
save_heatmap(
  mat_z_t = mat_z_t,
  anno_row = anno_row,
  cfg = config
)

writeLines(
  capture.output(sessionInfo()),
  config$session_info_file
)

message("\nStimulated latent DEG heatmap complete.")
message("Target matched_cluster: ", config$target_cluster)
message("Saved PDF heatmap to: ", config$heatmap_pdf)
message("Saved PNG heatmap to: ", config$heatmap_png)
message("Saved averaged expression matrix to: ", config$averaged_expression_csv)
message("Saved z-score matrix to: ", config$zscore_matrix_csv)
message("Saved gene presence table to: ", config$gene_presence_csv)
message("Saved sample annotation to: ", config$sample_annotation_csv)
message("Saved session info to: ", config$session_info_file)