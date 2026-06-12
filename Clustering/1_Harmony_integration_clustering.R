#!/usr/bin/env Rscript

# ==============================================================================
# Harmony integration, clustering, and annotation of EBV CD8+ T cell objects
# ==============================================================================
#
# Purpose:
#   This script performs Harmony integration and clustering for publication-ready
#   EBV CD8+ T cell Seurat objects. It supports:
#
#     1. baseline EBV-enriched CD8+ T cells,
#     2. activated EBV/LCL-responsive paired CD8+ T cells.
#
#   It performs:
#     - run/chip annotation,
#     - SCTransform with cell cycle and mitochondrial regression,
#     - Harmony integration by single-cell run,
#     - UMAP and clustering,
#     - cluster marker detection,
#     - contaminant cluster removal,
#     - illegitimate tetramer/sample filtering for sticky multimers,
#     - final reclustering,
#     - cluster renumbering and cell type annotation,
#     - canonical marker feature plots,
#     - activated CITE-seq feature plots,
#     - cytotoxicity module scoring for activated data,
#     - final publication RDS export.
#
# Expected inputs:
#   Baseline object:
#     13_EBV_baseline_TCR_publication.rds
#
#   Activated paired object:
#     13_EBV_LCL_activated_TCR_publication.rds
#
# Notes:
#   - Keep project-specific paths in the CONFIG section only.
#   - The baseline illegitimate tetramer barcode file can be reused for the
#     activated object so the same cells are excluded consistently.
#   - This script does not perform TCR QC. It assumes that has already been done.
#
# ==============================================================================


# ----------------------------- #
# 0. Load packages
# ----------------------------- #

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(ggplot2)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

run_map <- list(
  "1"  = c("GEMA", "GEMB"),
  "2"  = c("GEM1", "GEM2"),
  "3"  = c("GEM3", "GEM4"),
  "4"  = c("GEM5", "GEM6"),
  "5"  = c("GEM7", "GEM8"),
  "6"  = c("GEM9", "GEM10"),
  "7"  = c("GEM11"),
  "8"  = c("GEM13"),
  "9"  = c("GEMEBV1"),
  "10" = c("GEMEBV2", "GEMEBV3"),
  "11" = c("GEMEBV4", "GEMEBV5"),
  "12" = c("EBVLCL1DEX", "EBVLCL2DEX", "EBVLCL1CD8", "EBVLCL2CD8", "EBVLCL1CD4", "EBVLCL2CD4")
)

baseline_marker_genes <- c(
  "BACH2", "LTB", "CCR7", "CD55", "GZMK", "TC2N", "CXCR3", "EPHA4",
  "GZMH", "NKG7", "HLA-DRB1", "HLA-DPB1", "GNLY", "GZMB", "FGFBP2",
  "PRF1", "TPM4", "CD69", "ACTB", "GBP1", "KLRB1", "CEBPD", "CCR6",
  "SLC4A10", "CCL4", "CCL4L2", "IFNG", "TNF"
)

activated_marker_genes <- c(
  "BACH2", "LTB", "CCR7", "CD55", "GZMK", "TC2N", "CXCR3", "EPHA4",
  "FLNA", "MYADM", "FOSB", "IL2RB", "GZMH", "NKG7", "HLA-DRB1",
  "HLA-DPB1", "ISG15", "IFIT1", "STAT1", "GNLY", "GZMB", "FGFBP2",
  "PRF1", "TPM4", "CD69", "ACTB", "GBP1", "KLRB1", "CEBPD", "CCR6",
  "SLC4A10", "MIR155HG", "TNFRSF9", "FABP5", "CD82", "RGS1", "FAS",
  "CD70", "IL32", "CCL4", "CCL4L2", "IFNG", "TNF", "IL2RA"
)

canonical_feature_genes <- c(
  "CCR7", "IL7R", "TCF7", "CD28", "GZMK", "CD69", "IL2RA", "GZMB"
)

cytotoxicity_genes <- c(
  "GZMB", "GZMH", "GZMK", "GZMA", "GZMM",
  "PRF1", "GNLY", "NKG7", "FGFBP2", "CTSW",
  "KLRD1", "KLRB1", "KLRK1",
  "TBX21", "EOMES", "ZEB2", "ID2",
  "IFNG", "IFITM1", "IFITM2", "IFITM3",
  "CX3CR1", "S100A4", "S100A10"
)

allowed_tetramer_sample_map <- tibble::tribble(
  ~tet_core, ~sample,
  "LPRR", "C246",
  "LPRR", "C287",
  "RPQG", "C246",
  "RPQG", "C287",
  "ELRS", "MS359",
  "ELRS", "MS394",
  "QEIR", "MS173",
  "LPFE", "C339",
  "LPFE", "MS069",
  "RYGF", "MS045",
  "RYGF", "MS136",
  "LEKA", "MS358",
  "LEKA", "MS173"
)

datasets <- list(
  list(
    dataset_id = "baseline_ebv",
    input_rds = "path/to/input/13_EBV_baseline_TCR_publication.rds",
    output_rds = "path/to/output/publication_RDS/14_baseline_EBV_clustered.rds",
    output_filtered_rds = "path/to/output/publication_RDS/14_baseline_EBV_filtered_pre_annotation.rds",
    
    plot_dir = "path/to/output/Publication_data/EBV_baseline/Clustering",
    marker_csv = "path/to/output/Publication_data/EBV_baseline/Clustering/top30_cluster_markers_post_subset.csv",
    avg_expression_csv = "path/to/output/Publication_data/EBV_baseline/Clustering/avg_expression_markers_by_cluster.csv",
    illegitimate_tetramer_csv = "path/to/output/Publication_data/EBV_baseline/Clustering/flagged_illegitimate_tetramer_cells.csv",
    illegitimate_tetramer_barcodes_txt = "path/to/output/Publication_data/EBV_baseline/Clustering/flagged_illegitimate_tetramer_barcodes.txt",
    session_info_file = "path/to/output/Publication_data/EBV_baseline/Clustering/sessionInfo_baseline_clustering.txt",
    
    run_map = run_map,
    include_activated_run = FALSE,
    
    first_pass_resolution = 0.4,
    final_resolution = 0.5,
    dims = 1:30,
    contaminant_clusters_first_pass = c("5"),
    
    perform_illegitimate_tetramer_filter = TRUE,
    allowed_tetramer_sample_map = allowed_tetramer_sample_map,
    remove_barcodes_file = NULL,
    
    cluster_renumber_map = c(
      "2" = "0",
      "0" = "1",
      "1" = "2",
      "3" = "3",
      "4" = "4",
      "6" = "5",
      "5" = "6"
    ),
    
    cluster_label_col = "new_cluster",
    cluster_cols = c(
      "#E65757", "#6CC3F4", "#B560DD", "#9BE599",
      "#DD9560", "#c90076", "#f467ba"
    ),
    
    celltype_map = c(
      "0" = "Naive/early TCM",
      "1" = "TEM",
      "2" = "Late TEM",
      "3" = "CTL",
      "4" = "CD69+ TCM",
      "5" = "CD69+ TCM",
      "6" = "CD69+ TEM"
    ),
    
    marker_genes = baseline_marker_genes,
    feature_genes = canonical_feature_genes,
    
    make_adt_activation_plots = FALSE,
    make_cytotoxicity_score = FALSE,
    make_latent_group_plot = FALSE
  ),
  
  list(
    dataset_id = "activated_ebv",
    input_rds = "path/to/input/13_EBV_LCL_activated_TCR_publication.rds",
    output_rds = "path/to/output/publication_RDS/14_activated_EBV_clustered_module_scored.rds",
    output_filtered_rds = "path/to/output/publication_RDS/14_activated_EBV_filtered_pre_annotation.rds",
    
    plot_dir = "path/to/output/Publication_data/EBV_activated/Clustering",
    marker_csv = "path/to/output/Publication_data/EBV_activated/Clustering/top30_cluster_markers_post_subset.csv",
    avg_expression_csv = "path/to/output/Publication_data/EBV_activated/Clustering/avg_expression_markers_by_cluster.csv",
    illegitimate_tetramer_csv = "path/to/output/Publication_data/EBV_activated/Clustering/flagged_illegitimate_tetramer_cells.csv",
    illegitimate_tetramer_barcodes_txt = "path/to/output/Publication_data/EBV_activated/Clustering/flagged_illegitimate_tetramer_barcodes.txt",
    session_info_file = "path/to/output/Publication_data/EBV_activated/Clustering/sessionInfo_activated_clustering.txt",
    
    run_map = run_map,
    include_activated_run = TRUE,
    
    first_pass_resolution = 0.5,
    final_resolution = 0.5,
    dims = 1:30,
    contaminant_clusters_first_pass = c("8"),
    
    perform_illegitimate_tetramer_filter = FALSE,
    allowed_tetramer_sample_map = allowed_tetramer_sample_map,
    
    # Use the baseline barcode file here if you want the same sticky tetramer
    # cells removed from the activated paired object.
    remove_barcodes_file = "path/to/output/Publication_data/EBV_baseline/Clustering/flagged_illegitimate_tetramer_barcodes.txt",
    
    cluster_renumber_map = c(
      "0" = "0",
      "1" = "2",
      "2" = "1",
      "3" = "7",
      "4" = "3",
      "5" = "5",
      "6" = "6",
      "7" = "4",
      "8" = "8"
    ),
    
    cluster_label_col = "matched_cluster",
    cluster_cols = c(
      "#6CC3F4", "#B560DD", "#E65757", "#9BE599", "#DD9560",
      "#f467ba", "#c90076", "#72BEB7", "#49d10b"
    ),
    
    celltype_map = c(
      "0" = "TEM",
      "1" = "T Late EM",
      "2" = "Naive/early TCM",
      "3" = "CTL",
      "4" = "T activated CM",
      "5" = "T activated EM",
      "6" = "T CD69+ CM",
      "7" = "T activated EM",
      "8" = "TEM"
    ),
    
    marker_genes = activated_marker_genes,
    feature_genes = canonical_feature_genes,
    
    make_adt_activation_plots = TRUE,
    adt_activation_features = c("C0146-anti-human-CD69", "C0355-anti-human-CD137"),
    
    make_cytotoxicity_score = TRUE,
    cytotoxicity_genes = cytotoxicity_genes,
    
    make_latent_group_plot = TRUE
  )
)

global_config <- list(
  assay_rna = "RNA",
  assay_sct = "SCT",
  assay_adt = "ADT",
  
  vars_to_regress = c("S.Score", "G2M.Score", "percent.mt"),
  return_only_var_genes = TRUE,
  
  reduction_harmony = "harmony",
  reduction_umap = "umap.harmony",
  graph_name = "harmony_snn",
  cluster_name = "harmony_clusters",
  
  marker_min_log2fc_initial = 0.5,
  marker_min_log2fc_top = 1,
  top_n_markers = 30,
  marker_exclude_pattern = "^RPS|^RPL|^MT-",
  
  feature_plot_cols = c("#ffe1e1", "#c72c2c"),
  
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

setup_dataset_dirs <- function(cfg) {
  create_dir(cfg$plot_dir)
  create_parent_dir(cfg$output_rds)
  create_parent_dir(cfg$output_filtered_rds)
  create_parent_dir(cfg$marker_csv)
  create_parent_dir(cfg$avg_expression_csv)
  create_parent_dir(cfg$illegitimate_tetramer_csv)
  create_parent_dir(cfg$illegitimate_tetramer_barcodes_txt)
  create_parent_dir(cfg$session_info_file)
}

load_seurat_object <- function(path) {
  if (!file.exists(path)) {
    stop("Input RDS does not exist: ", path)
  }
  
  obj <- readRDS(path)
  
  if (!inherits(obj, "Seurat")) {
    stop("Input RDS must contain a Seurat object.")
  }
  
  obj
}

join_layers_if_possible <- function(assay_obj) {
  tryCatch(
    JoinLayers(assay_obj),
    error = function(e) assay_obj
  )
}

assign_run_from_batch <- function(obj, cfg) {
  batch_values <- as.character(obj$batch)
  run_values <- rep(NA_character_, length(batch_values))
  
  run_map_use <- cfg$run_map
  
  if (!isTRUE(cfg$include_activated_run)) {
    run_map_use <- run_map_use[names(run_map_use) != "12"]
  }
  
  for (run_id in names(run_map_use)) {
    run_values[batch_values %in% run_map_use[[run_id]]] <- run_id
  }
  
  obj$run <- factor(run_values)
  
  if (any(is.na(obj$run))) {
    warning(
      "Some cells have NA run assignments. Check batch names: ",
      paste(unique(batch_values[is.na(obj$run)]), collapse = ", ")
    )
  }
  
  obj
}

prepare_rna_for_sct <- function(obj, cfg, global_config) {
  DefaultAssay(obj) <- global_config$assay_rna
  
  obj[[global_config$assay_rna]] <- join_layers_if_possible(
    obj[[global_config$assay_rna]]
  )
  
  obj[[global_config$assay_rna]] <- split(
    obj[[global_config$assay_rna]],
    f = obj$run
  )
  
  obj
}

run_harmony_clustering <- function(obj, cfg, global_config, run_label = "final") {
  message("Running SCTransform: ", cfg$dataset_id, " / ", run_label)
  
  obj <- assign_run_from_batch(obj, cfg)
  obj <- prepare_rna_for_sct(obj, cfg, global_config)
  
  obj <- SCTransform(
    object = obj,
    assay = global_config$assay_rna,
    new.assay.name = global_config$assay_sct,
    vars.to.regress = global_config$vars_to_regress,
    return.only.var.genes = global_config$return_only_var_genes,
    verbose = global_config$verbose
  )
  
  DefaultAssay(obj) <- global_config$assay_sct
  
  obj[[global_config$assay_rna]] <- join_layers_if_possible(
    obj[[global_config$assay_rna]]
  )
  
  obj <- RunPCA(
    object = obj,
    assay = global_config$assay_sct,
    verbose = global_config$verbose
  )
  
  elbow_plot <- ElbowPlot(obj)
  
  ggsave(
    filename = file.path(cfg$plot_dir, paste0("elbow_", run_label, ".jpeg")),
    plot = elbow_plot,
    height = 5,
    width = 20,
    dpi = 300
  )
  
  obj <- RunHarmony(
    object = obj,
    group.by.vars = "run",
    assay.use = global_config$assay_sct,
    verbose = global_config$verbose
  )
  
  obj <- FindNeighbors(
    object = obj,
    reduction = global_config$reduction_harmony,
    dims = cfg$dims,
    graph.name = global_config$graph_name,
    verbose = global_config$verbose
  )
  
  obj <- FindClusters(
    object = obj,
    resolution = ifelse(run_label == "first_pass", cfg$first_pass_resolution, cfg$final_resolution),
    graph.name = global_config$graph_name,
    cluster.name = global_config$cluster_name,
    verbose = global_config$verbose
  )
  
  obj <- RunUMAP(
    object = obj,
    reduction = global_config$reduction_harmony,
    dims = cfg$dims,
    reduction.name = global_config$reduction_umap,
    verbose = global_config$verbose
  )
  
  obj
}

save_cluster_markers <- function(obj, cfg, global_config, out_csv = cfg$marker_csv) {
  message("Finding cluster markers for: ", cfg$dataset_id)
  
  DefaultAssay(obj) <- global_config$assay_sct
  
  obj_prep <- PrepSCTFindMarkers(
    object = obj,
    assay = global_config$assay_sct
  )
  
  cluster_markers <- FindAllMarkers(
    object = obj_prep,
    only.pos = TRUE,
    assay = global_config$assay_sct,
    slot = "data"
  )
  
  markers_filtered <- cluster_markers |>
    group_by(cluster) |>
    filter(avg_log2FC > global_config$marker_min_log2fc_initial) |>
    ungroup()
  
  excluded_genes <- grep(
    global_config$marker_exclude_pattern,
    markers_filtered$gene,
    value = TRUE
  )
  
  markers_filtered <- markers_filtered |>
    filter(!gene %in% excluded_genes)
  
  top_markers <- markers_filtered |>
    group_by(cluster) |>
    filter(avg_log2FC > global_config$marker_min_log2fc_top) |>
    slice_head(n = global_config$top_n_markers) |>
    ungroup()
  
  write.csv(top_markers, out_csv, row.names = FALSE)
  
  invisible(top_markers)
}

remove_contaminant_clusters <- function(obj, cfg, global_config) {
  if (length(cfg$contaminant_clusters_first_pass) == 0) {
    return(obj)
  }
  
  subset(
    obj,
    subset = !(harmony_clusters %in% cfg$contaminant_clusters_first_pass)
  )
}

extract_tetramer_core <- function(x) {
  x |>
    as.character() |>
    str_remove("\\*$") |>
    str_extract("[A-Z]{4}$")
}

flag_illegitimate_tetramer_cells <- function(obj, cfg) {
  meta <- obj@meta.data |>
    rownames_to_column("cell_barcode") |>
    mutate(tet_core = extract_tetramer_core(tetramer))
  
  tetramers_to_check <- unique(cfg$allowed_tetramer_sample_map$tet_core)
  
  bad_cells <- meta |>
    filter(tet_core %in% tetramers_to_check) |>
    anti_join(
      cfg$allowed_tetramer_sample_map,
      by = c("tet_core", "sample")
    ) |>
    select(cell_barcode, sample, tetramer, tet_core)
  
  write_csv(
    bad_cells,
    cfg$illegitimate_tetramer_csv
  )
  
  write_lines(
    bad_cells$cell_barcode,
    cfg$illegitimate_tetramer_barcodes_txt
  )
  
  bad_cells$cell_barcode
}

read_barcodes_to_remove <- function(cfg) {
  if (is.null(cfg$remove_barcodes_file) || is.na(cfg$remove_barcodes_file)) {
    return(character())
  }
  
  if (!file.exists(cfg$remove_barcodes_file)) {
    warning("Barcode removal file not found: ", cfg$remove_barcodes_file)
    return(character())
  }
  
  read_lines(cfg$remove_barcodes_file)
}

remove_barcodes <- function(obj, barcodes) {
  barcodes <- intersect(barcodes, colnames(obj))
  
  if (length(barcodes) == 0) {
    return(obj)
  }
  
  subset(
    obj,
    cells = setdiff(colnames(obj), barcodes)
  )
}

apply_cluster_renumbering <- function(obj, cfg, global_config) {
  old_clusters <- as.character(obj[[global_config$cluster_name, drop = TRUE]])
  new_clusters <- unname(cfg$cluster_renumber_map[old_clusters])
  
  obj[[cfg$cluster_label_col]] <- factor(
    new_clusters,
    levels = sort(unique(na.omit(new_clusters)))
  )
  
  obj
}

save_cluster_umap <- function(obj, cfg, global_config) {
  p <- DimPlot(
    object = obj,
    reduction = global_config$reduction_umap,
    group.by = cfg$cluster_label_col,
    label = TRUE,
    pt.size = 0.3,
    alpha = 0.5,
    label.box = TRUE,
    stroke.size = 0.4,
    raster = FALSE,
    repel = TRUE,
    cols = cfg$cluster_cols,
    label.size = 5
  )
  
  ggsave(
    plot = p,
    height = 8,
    width = 9,
    dpi = 500,
    filename = file.path(cfg$plot_dir, "umap_clusters.png"),
    bg = "transparent"
  )
}

save_average_marker_expression <- function(obj, cfg, global_config) {
  DefaultAssay(obj) <- global_config$assay_sct
  
  genes_present <- intersect(cfg$marker_genes, rownames(obj))
  genes_missing <- setdiff(cfg$marker_genes, genes_present)
  
  if (length(genes_missing) > 0) {
    message(
      "Missing genes in ",
      cfg$dataset_id,
      ": ",
      paste(genes_missing, collapse = ", ")
    )
  }
  
  avg <- AverageExpression(
    object = obj,
    assays = global_config$assay_sct,
    features = genes_present,
    group.by = cfg$cluster_label_col,
    slot = "data",
    verbose = FALSE
  )[[global_config$assay_sct]]
  
  avg_df <- as.data.frame(avg)
  avg_df$gene <- rownames(avg_df)
  
  write.csv(
    avg_df,
    cfg$avg_expression_csv,
    row.names = FALSE
  )
  
  invisible(avg_df)
}

save_feature_plots <- function(obj, cfg, global_config) {
  DefaultAssay(obj) <- global_config$assay_sct
  
  features_present <- intersect(cfg$feature_genes, rownames(obj))
  
  p <- FeaturePlot(
    object = obj,
    features = features_present,
    reduction = global_config$reduction_umap,
    ncol = 4,
    cols = global_config$feature_plot_cols,
    pt.size = 0.3
  )
  
  ggsave(
    plot = p,
    height = 8,
    width = 20,
    dpi = 300,
    filename = file.path(cfg$plot_dir, "canonical_subset_markers_umap.png"),
    bg = "transparent"
  )
}

assign_celltypes <- function(obj, cfg) {
  clusters <- as.character(obj[[cfg$cluster_label_col, drop = TRUE]])
  obj$celltype <- unname(cfg$celltype_map[clusters])
  obj
}

save_adt_activation_plots <- function(obj, cfg, global_config) {
  if (!isTRUE(cfg$make_adt_activation_plots)) {
    return(invisible(NULL))
  }
  
  if (!global_config$assay_adt %in% Assays(obj)) {
    warning("ADT assay not found. Skipping ADT activation plots.")
    return(invisible(NULL))
  }
  
  DefaultAssay(obj) <- global_config$assay_adt
  
  obj <- NormalizeData(
    object = obj,
    normalization.method = "CLR",
    verbose = global_config$verbose
  )
  
  features_present <- intersect(cfg$adt_activation_features, rownames(obj))
  
  if (length(features_present) == 0) {
    warning("No configured ADT activation features found.")
    return(invisible(NULL))
  }
  
  p <- FeaturePlot(
    object = obj,
    features = features_present,
    reduction = global_config$reduction_umap,
    ncol = 2,
    cols = global_config$feature_plot_cols,
    pt.size = 0.4
  )
  
  ggsave(
    plot = p,
    height = 8,
    width = 15,
    dpi = 300,
    filename = file.path(cfg$plot_dir, "CITEseq_activation_markers_umap.png"),
    bg = "transparent"
  )
  
  invisible(NULL)
}

add_cytotoxicity_score <- function(obj, cfg, global_config) {
  if (!isTRUE(cfg$make_cytotoxicity_score)) {
    return(obj)
  }
  
  DefaultAssay(obj) <- global_config$assay_sct
  
  genes_present <- intersect(cfg$cytotoxicity_genes, rownames(obj))
  
  obj <- AddModuleScore(
    object = obj,
    assay = global_config$assay_sct,
    features = list(Cytotoxicity = genes_present),
    name = "Cytotoxicity"
  )
  
  score_cols <- grep("^Cytotoxicity", colnames(obj@meta.data), value = TRUE)
  newest_score_col <- tail(score_cols, 1)
  
  obj$Cytotoxicity_score <- obj[[newest_score_col, drop = TRUE]]
  
  obj
}

save_latent_group_umap <- function(obj, cfg, global_config) {
  if (!isTRUE(cfg$make_latent_group_plot)) {
    return(invisible(NULL))
  }
  
  md <- obj@meta.data
  
  obj$latent_group <- case_when(
    !is.na(md$batch) &
      grepl("^GEMEBV", md$batch) &
      md$lifecycle == "Latent" ~ "Baseline Latent",
    
    !is.na(md$batch) &
      grepl("^EBVLCL", md$batch) &
      md$lifecycle == "Latent" ~ "Stimulated Latent",
    
    TRUE ~ "Non-Latent"
  )
  
  obj$latent_group <- factor(
    obj$latent_group,
    levels = c("Non-Latent", "Baseline Latent", "Stimulated Latent")
  )
  
  p <- DimPlot(
    object = obj,
    reduction = global_config$reduction_umap,
    group.by = "latent_group",
    pt.size = 0.5,
    raster = FALSE
  ) +
    scale_color_manual(
      values = c(
        "Non-Latent" = "grey85",
        "Baseline Latent" = "#F467BA",
        "Stimulated Latent" = "#B560DD"
      )
    ) +
    theme_classic() +
    theme(
      legend.title = element_blank(),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  ggsave(
    plot = p,
    height = 8,
    width = 9,
    dpi = 300,
    filename = file.path(cfg$plot_dir, "umap_latent_stim_vs_unstim.png"),
    bg = "transparent"
  )
  
  invisible(NULL)
}


# ----------------------------- #
# 3. Dataset-specific workflow
# ----------------------------- #

process_dataset <- function(cfg, global_config) {
  setup_dataset_dirs(cfg)
  
  message("\n==============================")
  message("Processing dataset: ", cfg$dataset_id)
  message("==============================")
  
  obj <- load_seurat_object(cfg$input_rds)
  
  # Optional pre-removal of known bad barcodes from another object.
  pre_remove_barcodes <- read_barcodes_to_remove(cfg)
  obj <- remove_barcodes(obj, pre_remove_barcodes)
  
  message("First-pass Harmony integration...")
  first_pass <- run_harmony_clustering(
    obj = obj,
    cfg = cfg,
    global_config = global_config,
    run_label = "first_pass"
  )
  
  first_pass_marker_csv <- file.path(
    cfg$plot_dir,
    "top30_cluster_markers_first_pass.csv"
  )
  
  save_cluster_markers(
    obj = first_pass,
    cfg = cfg,
    global_config = global_config,
    out_csv = first_pass_marker_csv
  )
  
  message("Removing contaminant clusters if configured...")
  obj_clean <- remove_contaminant_clusters(
    obj = first_pass,
    cfg = cfg,
    global_config = global_config
  )
  
  # For baseline EBV, flag sticky/illegitimate tetramer calls and remove them
  # before final re-integration.
  if (isTRUE(cfg$perform_illegitimate_tetramer_filter)) {
    message("Flagging illegitimate tetramer/sample combinations...")
    
    bad_tetramer_barcodes <- flag_illegitimate_tetramer_cells(
      obj = obj_clean,
      cfg = cfg
    )
    
    obj_clean <- remove_barcodes(
      obj = obj_clean,
      barcodes = bad_tetramer_barcodes
    )
  }
  
  saveRDS(
    obj_clean,
    cfg$output_filtered_rds
  )
  
  message("Final Harmony integration after filtering...")
  merged <- run_harmony_clustering(
    obj = obj_clean,
    cfg = cfg,
    global_config = global_config,
    run_label = "final"
  )
  
  message("Applying cluster renumbering and annotations...")
  merged <- apply_cluster_renumbering(
    obj = merged,
    cfg = cfg,
    global_config = global_config
  )
  
  save_cluster_umap(
    obj = merged,
    cfg = cfg,
    global_config = global_config
  )
  
  save_cluster_markers(
    obj = merged,
    cfg = cfg,
    global_config = global_config,
    out_csv = cfg$marker_csv
  )
  
  save_average_marker_expression(
    obj = merged,
    cfg = cfg,
    global_config = global_config
  )
  
  save_feature_plots(
    obj = merged,
    cfg = cfg,
    global_config = global_config
  )
  
  merged <- assign_celltypes(
    obj = merged,
    cfg = cfg
  )
  
  save_adt_activation_plots(
    obj = merged,
    cfg = cfg,
    global_config = global_config
  )
  
  merged <- add_cytotoxicity_score(
    obj = merged,
    cfg = cfg,
    global_config = global_config
  )
  
  save_latent_group_umap(
    obj = merged,
    cfg = cfg,
    global_config = global_config
  )
  
  saveRDS(
    merged,
    cfg$output_rds
  )
  
  writeLines(
    capture.output(sessionInfo()),
    cfg$session_info_file
  )
  
  message("Saved clustered object to: ", cfg$output_rds)
  message("Saved plots to: ", cfg$plot_dir)
  message("Saved marker CSV to: ", cfg$marker_csv)
  message("Saved session info to: ", cfg$session_info_file)
  
  invisible(merged)
}


# ----------------------------- #
# 4. Run pipeline
# ----------------------------- #

results <- lapply(
  datasets,
  process_dataset,
  global_config = global_config
)

message("\nAll Harmony clustering workflows complete.")