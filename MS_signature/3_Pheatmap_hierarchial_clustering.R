#!/usr/bin/env Rscript

# ==============================================================================
# Unsupervised sample x gene heatmap for stimulated latent EBV-specific T cells
# ==============================================================================
#
# Purpose:
#   This script creates an unsupervised clustered heatmap of selected genes across
#   samples in stimulated latent EBV-specific T cells from an activated/baseline
#   integrated Seurat object.
#
#   The default analysis:
#     - Loads 3_activated_EBV_module_scored_reannotated.rds
#     - Defines latent_group from batch and lifecycle
#     - Subsets to Stimulated Latent cells
#     - Subsets to new_cluster == "2"
#     - Averages expression per sample
#     - Z-scores each gene across samples
#     - Plots samples as rows and genes as columns using pheatmap
#
# Notes:
#   - Paths are intentionally generic for GitHub.
#   - Update the config block before running.
#   - The donor/sample column defaults to `sample`; change to `id` if needed.
#   - By default, expression is taken from the SCT assay scale.data layer/slot.
#   - If genes are not present in scale.data, the script can optionally run
#     ScaleData() for the requested genes in the subset.
#
# ==============================================================================


# ----------------------------- #
# 0. Load packages
# ----------------------------- #

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
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
  
  # Input RDS
  input_rds = "path/to/Publication_RDS/Activated_EBV_integrated_seurat.rds",
  
  # Output directory
  output_dir = "path/to/Publication_data/EBV_activated/heatmaps",
  
  # Metadata columns
  donor_col = "sample",
  cohort_col = "cohort",
  cluster_col = "new_cluster",
  batch_col = "batch",
  lifecycle_col = "lifecycle",
  
  # Latent group definitions
  baseline_batch_pattern = "^GEMEBV",
  stimulated_batch_pattern = "^EBVLCL",
  latent_lifecycle_label = "Latent",
  baseline_label = "Baseline Latent",
  stimulated_label = "Stimulated Latent",
  non_latent_label = "Non-Latent",
  
  # Analysis subset
  target_latent_group = "Stimulated Latent",
  target_cluster = "2",
  cohort_order = c("Control", "MS"),
  
  # Expression settings
  assay = "SCT",
  expression_layer = "scale.data",
  
  # If TRUE, run ScaleData() on missing heatmap genes when expression_layer is scale.data.
  scale_missing_genes = TRUE,
  
  # Genes to plot
  genes = c(
    "FTH1", "JAK1", "ITM2B", "RPL13", "SMCHD1", "IL7R", "HLA-B", "PTMA",
    "HERC5", "TOMM7", "MT-CYB", "H3-3A", "HLA-F", "RPS10", "MT-ND4L",
    "HNRNPA1", "FYN", "PFN1", "RPS29", "MT-CO3", "RPS15", "RPS18",
    "SLFN5", "ZBTB20", "RPLP2", "MT-ND3", "ARPC2", "IFI6", "CD3D",
    "RPS27", "HLA-C", "CXCR4", "RPS27A", "CALM1", "CCSER2", "ANKRD12",
    "VIM", "RPS9", "CDC42SE2", "EEF1G", "RPS6", "PDE3B", "EEF1A1",
    "RPS16", "YWHAB", "RPL30", "SARAF", "SRSF7"
  ),
  
  # Heatmap colours
  heatmap_colours = c(
    "#011f4b", "#03396c", "#005b96", "#6497b1", "#b3cde0",
    "#ffefea", "#fbd9d3", "#ffb09c", "#fe5757", "#cb2424", "#900000"
  ),
  
  annotation_colours = list(
    cohort = c(
      "Control" = "#6CC3F4",
      "MS" = "#E65757"
    )
  ),
  
  # Heatmap settings
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 9,
  fontsize_col = 10,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = NA,
  
  pdf_width = 15,
  pdf_height = 8,
  png_width = 15,
  png_height = 8,
  png_dpi = 400,
  
  # Output files
  filtered_metadata_csv = "StimulatedLatent_cluster2_heatmap_filtered_metadata.csv",
  genes_requested_csv = "StimulatedLatent_cluster2_heatmap_requested_genes.csv",
  sample_average_expression_csv = "StimulatedLatent_cluster2_heatmap_sample_average_expression.csv",
  sample_gene_zscore_matrix_csv = "StimulatedLatent_cluster2_heatmap_sample_gene_zscore_matrix.csv",
  heatmap_pdf = "StimulatedLatent_cluster2_gene_by_sample_heatmap.pdf",
  heatmap_png = "StimulatedLatent_cluster2_gene_by_sample_heatmap.png",
  session_info_file = "sessionInfo_StimulatedLatent_cluster2_gene_by_sample_heatmap.txt"
)


# ----------------------------- #
# 2. Helper functions
# ----------------------------- #

create_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}


load_seurat_object <- function(path) {
  
  if (!file.exists(path)) {
    stop("Input RDS file does not exist: ", path)
  }
  
  obj <- readRDS(path)
  
  if (!inherits(obj, "Seurat")) {
    stop("Input RDS does not contain a Seurat object: ", path)
  }
  
  obj
}


check_required_columns <- function(df, required_cols, object_name = "dataframe") {
  
  missing_cols <- setdiff(required_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ",
      object_name,
      ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  invisible(TRUE)
}


get_assay_data_compat <- function(object, assay, layer_name) {
  
  # Compatible with Seurat v5 and Seurat v4.
  tryCatch(
    {
      GetAssayData(object, assay = assay, layer = layer_name)
    },
    error = function(e) {
      GetAssayData(object, assay = assay, slot = layer_name)
    }
  )
}


clean_chr <- function(x) {
  
  x <- str_trim(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL", "null", "None", "none")] <- NA_character_
  x
}


clean_cohort <- function(x, cohort_order = c("Control", "MS")) {
  
  x_raw <- as.character(x)
  x_clean <- str_to_lower(str_trim(x_raw))
  
  out <- case_when(
    x_clean %in% c(
      "control", "controls", "healthy control", "healthy controls",
      "ctrl", "hc", "nms", "non-ms", "nonms", "non_ms"
    ) ~ "Control",
    
    x_clean %in% c(
      "ms", "rrms", "pwms", "multiple sclerosis"
    ) ~ "MS",
    
    TRUE ~ x_raw
  )
  
  factor(out, levels = cohort_order)
}


safe_zscore_rows <- function(mat) {
  
  mat_z <- t(scale(t(as.matrix(mat))))
  mat_z[is.na(mat_z)] <- 0
  mat_z
}


safe_filename <- function(x) {
  
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+$", "") %>%
    str_replace_all("^_+", "")
}


# ----------------------------- #
# 3. Load object and prepare metadata
# ----------------------------- #

create_dir(config$output_dir)

message("Loading activated EBV object:")
message(config$input_rds)

obj <- load_seurat_object(config$input_rds)

required_cols <- c(
  config$donor_col,
  config$cohort_col,
  config$cluster_col,
  config$batch_col,
  config$lifecycle_col
)

check_required_columns(
  df = obj@meta.data,
  required_cols = required_cols,
  object_name = "activated EBV metadata"
)

if (!config$assay %in% Assays(obj)) {
  stop(
    "Assay not found: ",
    config$assay,
    ". Available assays: ",
    paste(Assays(obj), collapse = ", ")
  )
}

DefaultAssay(obj) <- config$assay

meta_clean <- obj@meta.data %>%
  as.data.frame() %>%
  mutate(
    donor_id = clean_chr(.data[[config$donor_col]]),
    cohort_clean = clean_cohort(
      .data[[config$cohort_col]],
      cohort_order = config$cohort_order
    ),
    cluster_clean = clean_chr(.data[[config$cluster_col]]),
    batch_clean = clean_chr(.data[[config$batch_col]]),
    lifecycle_clean = clean_chr(.data[[config$lifecycle_col]]),
    
    latent_group = case_when(
      !is.na(batch_clean) &
        grepl(config$baseline_batch_pattern, batch_clean) &
        lifecycle_clean == config$latent_lifecycle_label ~ config$baseline_label,
      
      !is.na(batch_clean) &
        grepl(config$stimulated_batch_pattern, batch_clean) &
        lifecycle_clean == config$latent_lifecycle_label ~ config$stimulated_label,
      
      TRUE ~ config$non_latent_label
    ),
    
    latent_group = factor(
      latent_group,
      levels = c(
        config$non_latent_label,
        config$baseline_label,
        config$stimulated_label
      )
    )
  )

rownames(meta_clean) <- colnames(obj)
obj@meta.data <- meta_clean

message("Latent group counts:")
print(table(obj$latent_group, useNA = "ifany"))

message("Cluster counts:")
print(table(obj$cluster_clean, useNA = "ifany"))


# ----------------------------- #
# 4. Subset to target cells
# ----------------------------- #

cells_keep <- rownames(obj@meta.data)[
  obj$latent_group == config$target_latent_group &
    obj$cluster_clean == config$target_cluster &
    obj$cohort_clean %in% config$cohort_order &
    !is.na(obj$donor_id) &
    obj$donor_id != ""
]

if (length(cells_keep) == 0) {
  stop("No cells matched the target latent group / cluster / cohort filters.")
}

obj_sub <- subset(
  obj,
  cells = cells_keep
)

DefaultAssay(obj_sub) <- config$assay

filtered_metadata_path <- file.path(
  config$output_dir,
  config$filtered_metadata_csv
)

write_csv(
  obj_sub@meta.data %>%
    as.data.frame() %>%
    rownames_to_column("cell_barcode"),
  filtered_metadata_path
)

message("Cells retained for heatmap:")
print(ncol(obj_sub))

message("Cells by cohort:")
print(table(obj_sub$cohort_clean, useNA = "ifany"))

message("Cells by sample:")
print(table(obj_sub$donor_id, useNA = "ifany"))


# ----------------------------- #
# 5. Check genes and expression data
# ----------------------------- #

available_genes <- rownames(obj_sub[[config$assay]])

genes_present_in_assay <- intersect(config$genes, available_genes)
genes_missing_from_assay <- setdiff(config$genes, available_genes)

if (length(genes_present_in_assay) < 2) {
  stop("Too few requested genes were found in assay: ", config$assay)
}

expr_mat <- get_assay_data_compat(
  object = obj_sub,
  assay = config$assay,
  layer_name = config$expression_layer
)

genes_present_in_expr <- intersect(genes_present_in_assay, rownames(expr_mat))
genes_missing_from_expr <- setdiff(genes_present_in_assay, rownames(expr_mat))

if (
  length(genes_missing_from_expr) > 0 &&
  config$expression_layer == "scale.data" &&
  isTRUE(config$scale_missing_genes)
) {
  
  message(
    "Scaling missing genes in subset because they are not present in scale.data: ",
    paste(genes_missing_from_expr, collapse = ", ")
  )
  
  obj_sub <- ScaleData(
    obj_sub,
    assay = config$assay,
    features = genes_present_in_assay,
    verbose = FALSE
  )
  
  expr_mat <- get_assay_data_compat(
    object = obj_sub,
    assay = config$assay,
    layer_name = config$expression_layer
  )
  
  genes_present_in_expr <- intersect(genes_present_in_assay, rownames(expr_mat))
  genes_missing_from_expr <- setdiff(genes_present_in_assay, rownames(expr_mat))
}

if (length(genes_present_in_expr) < 2) {
  stop(
    "Too few requested genes were available in expression layer/slot: ",
    config$expression_layer
  )
}

genes_requested_path <- file.path(
  config$output_dir,
  config$genes_requested_csv
)

write_csv(
  tibble(
    gene = config$genes,
    present_in_assay = config$genes %in% available_genes,
    present_in_expression_layer = config$genes %in% rownames(expr_mat)
  ),
  genes_requested_path
)

if (length(genes_missing_from_assay) > 0) {
  warning(
    "These requested genes were not found in the assay and will be omitted: ",
    paste(genes_missing_from_assay, collapse = ", ")
  )
}

if (length(genes_missing_from_expr) > 0) {
  warning(
    "These requested genes were not found in the expression layer/slot and will be omitted: ",
    paste(genes_missing_from_expr, collapse = ", ")
  )
}


# ----------------------------- #
# 6. Average expression per sample
# ----------------------------- #
#
# Expression matrix is genes x cells. We average across cells within donor_id.
#
# ----------------------------- #

meta_sub <- obj_sub@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell_barcode")

sample_levels <- sort(unique(meta_sub$donor_id))
sample_factor <- factor(meta_sub$donor_id, levels = sample_levels)

aggregation_matrix <- Matrix::sparseMatrix(
  i = seq_along(sample_factor),
  j = as.integer(sample_factor),
  x = 1,
  dims = c(length(sample_factor), nlevels(sample_factor)),
  dimnames = list(meta_sub$cell_barcode, sample_levels)
)

cell_counts_per_sample <- Matrix::colSums(aggregation_matrix)

expr_use <- expr_mat[genes_present_in_expr, meta_sub$cell_barcode, drop = FALSE]

sample_sum <- expr_use %*% aggregation_matrix

sample_average <- sweep(
  sample_sum,
  MARGIN = 2,
  STATS = cell_counts_per_sample,
  FUN = "/"
)

sample_average <- as.matrix(sample_average)

# Preserve requested gene order.
gene_order_present <- config$genes[config$genes %in% rownames(sample_average)]
sample_average <- sample_average[gene_order_present, , drop = FALSE]

sample_average_expression_path <- file.path(
  config$output_dir,
  config$sample_average_expression_csv
)

write_csv(
  as.data.frame(sample_average) %>%
    rownames_to_column("gene"),
  sample_average_expression_path
)


# ----------------------------- #
# 7. Z-score per gene and transpose
# ----------------------------- #

mat_z <- safe_zscore_rows(sample_average)

# Samples as rows, genes as columns.
mat_z_t <- t(mat_z)

sample_gene_zscore_matrix_path <- file.path(
  config$output_dir,
  config$sample_gene_zscore_matrix_csv
)

write_csv(
  as.data.frame(mat_z_t) %>%
    rownames_to_column("sample"),
  sample_gene_zscore_matrix_path
)


# ----------------------------- #
# 8. Row annotation
# ----------------------------- #

meta_samp <- meta_sub %>%
  distinct(donor_id, cohort_clean) %>%
  arrange(match(donor_id, rownames(mat_z_t)))

if (anyDuplicated(meta_samp$donor_id) > 0) {
  stop("Each donor/sample should have only one cohort annotation.")
}

anno_row <- data.frame(
  cohort = as.character(meta_samp$cohort_clean)
)

rownames(anno_row) <- meta_samp$donor_id
anno_row <- anno_row[rownames(mat_z_t), , drop = FALSE]

if (!identical(rownames(anno_row), rownames(mat_z_t))) {
  stop("Heatmap row annotation does not align with z-score matrix.")
}


# ----------------------------- #
# 9. Plot heatmap
# ----------------------------- #

heatmap_pdf_path <- file.path(
  config$output_dir,
  config$heatmap_pdf
)

heatmap_png_path <- file.path(
  config$output_dir,
  config$heatmap_png
)

heatmap_title <- paste0(
  config$target_latent_group,
  ", cluster ",
  config$target_cluster,
  ": selected gene expression by sample"
)

pheatmap_args <- list(
  mat = mat_z_t,
  cluster_rows = config$cluster_rows,
  cluster_cols = config$cluster_cols,
  annotation_row = anno_row,
  annotation_colors = config$annotation_colours,
  color = colorRampPalette(config$heatmap_colours)(200),
  fontsize_row = config$fontsize_row,
  fontsize_col = config$fontsize_col,
  show_rownames = config$show_rownames,
  show_colnames = config$show_colnames,
  border_color = config$border_color,
  main = heatmap_title
)

pdf(
  heatmap_pdf_path,
  width = config$pdf_width,
  height = config$pdf_height,
  useDingbats = FALSE
)

do.call(pheatmap::pheatmap, pheatmap_args)

dev.off()

png(
  filename = heatmap_png_path,
  width = config$png_width,
  height = config$png_height,
  units = "in",
  res = config$png_dpi
)

do.call(pheatmap::pheatmap, pheatmap_args)

dev.off()


# ----------------------------- #
# 10. Save session information
# ----------------------------- #

session_info_path <- file.path(
  config$output_dir,
  config$session_info_file
)

writeLines(
  c(
    "Configuration:",
    capture.output(str(config)),
    "",
    "Latent group counts:",
    capture.output(print(table(obj$latent_group, useNA = "ifany"))),
    "",
    "Cells retained for heatmap:",
    as.character(ncol(obj_sub)),
    "",
    "Cells by cohort:",
    capture.output(print(table(obj_sub$cohort_clean, useNA = "ifany"))),
    "",
    "Cells by sample:",
    capture.output(print(table(obj_sub$donor_id, useNA = "ifany"))),
    "",
    "Genes present in assay:",
    paste(genes_present_in_assay, collapse = ", "),
    "",
    "Genes present in expression layer:",
    paste(genes_present_in_expr, collapse = ", "),
    "",
    "Genes missing from assay:",
    paste(genes_missing_from_assay, collapse = ", "),
    "",
    "Genes missing from expression layer:",
    paste(genes_missing_from_expr, collapse = ", "),
    "",
    "Output files:",
    filtered_metadata_path,
    genes_requested_path,
    sample_average_expression_path,
    sample_gene_zscore_matrix_path,
    heatmap_pdf_path,
    heatmap_png_path,
    "",
    "Session information:",
    capture.output(sessionInfo())
  ),
  session_info_path
)


# ----------------------------- #
# 11. Completion messages
# ----------------------------- #

message("\nStimulated latent gene-by-sample heatmap workflow complete.")
message("Filtered metadata: ", filtered_metadata_path)
message("Requested genes: ", genes_requested_path)
message("Sample average expression: ", sample_average_expression_path)
message("Sample gene z-score matrix: ", sample_gene_zscore_matrix_path)
message("Heatmap PDF: ", heatmap_pdf_path)
message("Heatmap PNG: ", heatmap_png_path)
message("Session info: ", session_info_path)
