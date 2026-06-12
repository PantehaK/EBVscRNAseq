#!/usr/bin/env Rscript

# ==============================================================================
# hdWGCNA module analysis in stimulated latent EBV-specific T cells
# ==============================================================================
#
# Purpose:
#   This script uses hdWGCNA to identify weighted gene co-expression modules in
#   stimulated latent EBV-specific T cells from a selected matched cluster.
#
#   It then tests whether MS-associated differentially expressed genes from the
#   previous pseudobulk edgeR analysis are enriched within the identified modules.
#
#   It saves:
#     1. stimulated latent cluster-specific hdWGCNA Seurat object,
#     2. soft-power plots,
#     3. WGCNA dendrogram,
#     4. module eigengene sample-level table,
#     5. module eigengene MS vs Control statistics,
#     6. gene module/kME/degree table,
#     7. edgeR DEG overlap table by module,
#     8. shared DEG/module gene list,
#     9. session information.
#
# Expected input:
#   Activated paired EBV CD8+ object, for example:
#     15_activated_EBV_MS_stimLat_signature_scored.rds
#
# Optional input:
#   edgeR result from previous stimulated latent cluster pseudobulk analysis:
#     Age_adjusted_EdgeR_StimulatedLatent_MS_vs_Control_by_matched_cluster_ALL.csv
#
# Notes:
#   - The default target cluster is matched_cluster == "2".
#   - The default WGCNA group is stimulated latent cells only.
#   - DE genes are pulled from the edgeR result for the same matched cluster.
#   - If the edgeR CSV is unavailable, a fallback hard-coded gene set is used.
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
  library(ggplot2)
  library(patchwork)
  library(hdWGCNA)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  dataset_id = "activated_ebv_stimulated_latent_hdWGCNA",
  
  input_rds = "path/to/input/15_activated_EBV_MS_stimLat_signature_scored.rds",
  
  output_dir = "path/to/output/Publication_data/EBV_activated/hdWGCNA_stimLat_cluster2",
  
  output_wgcna_rds = "path/to/output/Publication_data/EBV_activated/hdWGCNA_stimLat_cluster2/24_stimLat_cluster2_hdWGCNA_object.rds",
  soft_power_plot_png = "path/to/output/Publication_data/EBV_activated/hdWGCNA_stimLat_cluster2/soft_power_plots.png",
  dendrogram_png = "path/to/output/Publication_data/EBV_activated/hdWGCNA_stimLat_cluster2/dendrogram_stimLat_cluster2.png",
  module_eigengene_cell_csv = "path/to/output/Publication_data/EBV_activated/hdWGCNA_stimLat_cluster2/module_eigengenes_cell_level.csv",
  module_eigengene_sample_csv = "path/to/output/Publication_data/EBV_activated/hdWGCNA_stimLat_cluster2/module_eigengenes_sample_level.csv",
  module_eigengene_stats_csv = "path/to/output/Publication_data/EBV_activated/hdWGCNA_stimLat_cluster2/module_eigengene_MS_vs_Control_stats.csv",
  gene_module_table_csv = "path/to/output/Publication_data/EBV_activated/hdWGCNA_stimLat_cluster2/stimLat_cluster2_gene_modules_kME_degree_and_stats.csv",
  edgeR_gene_set_csv = "path/to/output/Publication_data/EBV_activated/hdWGCNA_stimLat_cluster2/edgeR_gene_set_used_for_overlap.csv",
  module_overlap_csv = "path/to/output/Publication_data/EBV_activated/hdWGCNA_stimLat_cluster2/edgeR_DEG_overlap_by_module.csv",
  shared_genes_csv = "path/to/output/Publication_data/EBV_activated/hdWGCNA_stimLat_cluster2/shared_edgeR_DEGs_and_WGCNA_modules.csv",
  session_info_file = "path/to/output/Publication_data/EBV_activated/hdWGCNA_stimLat_cluster2/sessionInfo_hdWGCNA_stimLat_cluster2.txt",
  
  # Metadata columns.
  batch_col = "batch",
  lifecycle_col = "lifecycle",
  cohort_col = "cohort",
  sample_col = "sample",
  cluster_col = "matched_cluster",
  
  # Latent group definitions.
  latent_label = "Latent",
  non_latent_label = "Non-Latent",
  baseline_latent_label = "Baseline Latent",
  stimulated_latent_label = "Stimulated Latent",
  baseline_batch_pattern = "^GEMEBV",
  stimulated_batch_pattern = "^EBVLCL",
  
  # WGCNA subset.
  target_cluster = "2",
  cohorts_keep = c("MS", "Control"),
  
  # Assay and preprocessing.
  assay = "RNA",
  n_variable_features = 5000,
  npcs = 30,
  
  # hdWGCNA settings.
  wgcna_name = "stimLat_cluster2",
  gene_select = "fraction",
  fraction = 0.05,
  
  metacell_group_by = "sample",
  metacell_reduction = "pca",
  metacell_k = 3,
  metacell_max_shared = 10,
  metacell_min_cells = 3,
  metacell_ident_group = "sample",
  
  set_dat_expr_group_by = "sample",
  set_dat_expr_assay = "RNA",
  set_dat_expr_layer = "data",
  
  network_type = "signed",
  soft_power = 9,
  tom_name = "stimLat_cluster2_TOM",
  
  module_eigengene_group_by_vars = "sample",
  
  # edgeR DEG input from previous script.
  use_edgeR_de_csv = TRUE,
  edgeR_de_csv = "path/to/input/Age_adjusted_EdgeR_StimulatedLatent_MS_vs_Control_by_matched_cluster_ALL.csv",
  edgeR_cluster_col = "matched_cluster",
  edgeR_gene_col = "gene",
  edgeR_fdr_col = "p_val_adj",
  edgeR_logfc_col = "avg_log2FC",
  edgeR_fdr_threshold = 0.05,
  edgeR_abs_logfc_threshold = 0,
  
  # Fallback gene set if edgeR CSV does not exist.
  fallback_gene_set = c(
    "FTH1", "JAK1", "ITM2B", "RPL13", "SMCHD1", "IL7R", "HLA-B", "PTMA",
    "HERC5", "TOMM7", "MT-CYB", "H3-3A", "HLA-F", "RPS10", "MT-ND4L",
    "HNRNPA1", "FYN", "PFN1", "RPS29", "MT-CO3", "RPS15", "RPS18",
    "SLFN5", "ZBTB20", "RPLP2", "MT-ND3", "ARPC2", "IFI6", "CD3D",
    "RPS27", "HLA-C", "CXCR4", "RPS27A", "CALM1", "CCSER2", "ANKRD12",
    "VIM", "RPS9", "CDC42SE2", "EEF1G", "RPS6", "PDE3B", "EEF1A1",
    "RPS16", "YWHAB", "RPL30", "SARAF", "SRSF7"
  ),
  
  # Optional module of special interest.
  module_of_interest = "turquoise",
  
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
  create_parent_dir(cfg$output_wgcna_rds)
  create_parent_dir(cfg$soft_power_plot_png)
  create_parent_dir(cfg$dendrogram_png)
  create_parent_dir(cfg$module_eigengene_cell_csv)
  create_parent_dir(cfg$module_eigengene_sample_csv)
  create_parent_dir(cfg$module_eigengene_stats_csv)
  create_parent_dir(cfg$gene_module_table_csv)
  create_parent_dir(cfg$edgeR_gene_set_csv)
  create_parent_dir(cfg$module_overlap_csv)
  create_parent_dir(cfg$shared_genes_csv)
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


join_layers_if_possible <- function(obj, assay) {
  obj <- tryCatch(
    {
      JoinLayers(obj, assay = assay)
    },
    error = function(e) {
      obj
    }
  )
  
  obj
}


check_required_metadata <- function(obj, cfg) {
  required_cols <- c(
    cfg$batch_col,
    cfg$lifecycle_col,
    cfg$cohort_col,
    cfg$sample_col,
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


add_latent_group <- function(obj, cfg) {
  md <- obj@meta.data
  
  obj$latent_group <- case_when(
    !is.na(md[[cfg$batch_col]]) &
      grepl(cfg$baseline_batch_pattern, md[[cfg$batch_col]]) &
      md[[cfg$lifecycle_col]] == cfg$latent_label ~ cfg$baseline_latent_label,
    
    !is.na(md[[cfg$batch_col]]) &
      grepl(cfg$stimulated_batch_pattern, md[[cfg$batch_col]]) &
      md[[cfg$lifecycle_col]] == cfg$latent_label ~ cfg$stimulated_latent_label,
    
    TRUE ~ cfg$non_latent_label
  )
  
  obj$latent_group <- factor(
    obj$latent_group,
    levels = c(
      cfg$non_latent_label,
      cfg$baseline_latent_label,
      cfg$stimulated_latent_label
    )
  )
  
  obj
}


subset_stimulated_latent_cluster <- function(obj, cfg) {
  md <- obj@meta.data
  
  keep_cells <- rownames(md)[
    md$latent_group == cfg$stimulated_latent_label &
      md[[cfg$cohort_col]] %in% cfg$cohorts_keep &
      as.character(md[[cfg$cluster_col]]) == as.character(cfg$target_cluster)
  ]
  
  if (length(keep_cells) == 0) {
    stop(
      "No cells found for stimulated latent cluster ",
      cfg$target_cluster
    )
  }
  
  subset(obj, cells = keep_cells)
}


preprocess_for_hdWGCNA <- function(obj, cfg) {
  DefaultAssay(obj) <- cfg$assay
  
  obj <- join_layers_if_possible(obj, assay = cfg$assay)
  
  obj <- NormalizeData(
    object = obj,
    assay = cfg$assay,
    verbose = cfg$verbose
  )
  
  obj <- FindVariableFeatures(
    object = obj,
    assay = cfg$assay,
    nfeatures = cfg$n_variable_features,
    verbose = cfg$verbose
  )
  
  obj <- ScaleData(
    object = obj,
    assay = cfg$assay,
    features = VariableFeatures(obj),
    verbose = FALSE
  )
  
  obj <- RunPCA(
    object = obj,
    assay = cfg$assay,
    npcs = cfg$npcs,
    verbose = FALSE
  )
  
  obj
}


setup_hdWGCNA <- function(obj, cfg) {
  obj <- SetupForWGCNA(
    seurat_obj = obj,
    gene_select = cfg$gene_select,
    fraction = cfg$fraction,
    wgcna_name = cfg$wgcna_name
  )
  
  obj <- MetacellsByGroups(
    seurat_obj = obj,
    group.by = cfg$metacell_group_by,
    reduction = cfg$metacell_reduction,
    k = cfg$metacell_k,
    max_shared = cfg$metacell_max_shared,
    min_cells = cfg$metacell_min_cells,
    ident.group = cfg$metacell_ident_group
  )
  
  obj <- NormalizeMetacells(obj)
  
  mc <- GetMetacellObject(obj)
  metacell_groups <- sort(unique(mc[[cfg$set_dat_expr_group_by, drop = TRUE]]))
  
  obj <- SetDatExpr(
    seurat_obj = obj,
    group_name = metacell_groups,
    group.by = cfg$set_dat_expr_group_by,
    assay = cfg$set_dat_expr_assay,
    layer = cfg$set_dat_expr_layer
  )
  
  obj
}


run_network_construction <- function(obj, cfg) {
  obj <- TestSoftPowers(
    seurat_obj = obj,
    networkType = cfg$network_type
  )
  
  soft_power_plots <- PlotSoftPowers(obj)
  
  p_soft <- wrap_plots(soft_power_plots, ncol = 2)
  
  ggsave(
    filename = cfg$soft_power_plot_png,
    plot = p_soft,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  obj <- ConstructNetwork(
    seurat_obj = obj,
    soft_power = cfg$soft_power,
    tom_name = cfg$tom_name
  )
  
  png(
    filename = cfg$dendrogram_png,
    width = 6,
    height = 4,
    units = "in",
    res = 400
  )
  
  PlotDendrogram(obj)
  
  dev.off()
  
  obj
}


run_module_eigengenes_and_connectivity <- function(obj, cfg) {
  obj <- ModuleEigengenes(
    seurat_obj = obj,
    group.by.vars = cfg$module_eigengene_group_by_vars,
    wgcna_name = cfg$wgcna_name
  )
  
  obj <- ModuleConnectivity(
    seurat_obj = obj,
    group_name = cfg$tom_name,
    wgcna_name = cfg$wgcna_name
  )
  
  obj
}


get_module_eigengenes_safe <- function(obj, cfg) {
  MEs <- GetMEs(
    seurat_obj = obj,
    wgcna_name = cfg$wgcna_name,
    harmonized = TRUE
  )
  
  MEs <- MEs[colnames(obj), , drop = FALSE]
  
  MEs
}


summarise_module_eigengenes <- function(obj, MEs, cfg) {
  meta <- obj@meta.data[colnames(obj), c(cfg$sample_col, cfg$cohort_col), drop = FALSE]
  
  colnames(meta) <- c("sample", "cohort")
  
  ME_cell <- as.data.frame(MEs) |>
    rownames_to_column("cell_barcode") |>
    mutate(
      sample = meta$sample,
      cohort = meta$cohort
    )
  
  ME_sample <- ME_cell |>
    select(-cell_barcode) |>
    group_by(sample, cohort) |>
    summarise(
      across(where(is.numeric), mean, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(
    cell = ME_cell,
    sample = ME_sample
  )
}


wilcox_module_eigengenes <- function(ME_sample) {
  module_cols <- setdiff(colnames(ME_sample), c("sample", "cohort"))
  
  bind_rows(lapply(module_cols, function(module_name) {
    dat <- ME_sample |>
      select(sample, cohort, value = all_of(module_name)) |>
      filter(!is.na(value), !is.na(cohort))
    
    n_control <- sum(dat$cohort == "Control")
    n_ms <- sum(dat$cohort == "MS")
    
    if (n_control < 2 || n_ms < 2) {
      return(tibble(
        module = module_name,
        n_control = n_control,
        n_ms = n_ms,
        p_value = NA_real_,
        median_control = median(dat$value[dat$cohort == "Control"], na.rm = TRUE),
        median_ms = median(dat$value[dat$cohort == "MS"], na.rm = TRUE),
        median_difference_ms_minus_control = NA_real_
      ))
    }
    
    test <- wilcox.test(value ~ cohort, data = dat, exact = FALSE)
    
    median_control <- median(dat$value[dat$cohort == "Control"], na.rm = TRUE)
    median_ms <- median(dat$value[dat$cohort == "MS"], na.rm = TRUE)
    
    tibble(
      module = module_name,
      n_control = n_control,
      n_ms = n_ms,
      p_value = test$p.value,
      median_control = median_control,
      median_ms = median_ms,
      median_difference_ms_minus_control = median_ms - median_control
    )
  })) |>
    mutate(
      p_adj = p.adjust(p_value, method = "BH")
    ) |>
    arrange(p_adj)
}


standardise_module_table <- function(modules_df) {
  if (!"gene_name" %in% colnames(modules_df)) {
    possible_gene_cols <- intersect(
      c("gene", "Gene", "name", "gene_id"),
      colnames(modules_df)
    )
    
    if (length(possible_gene_cols) > 0) {
      modules_df <- modules_df |>
        rename(gene_name = all_of(possible_gene_cols[1]))
    } else {
      modules_df$gene_name <- rownames(modules_df)
    }
  }
  
  if (!"module" %in% colnames(modules_df)) {
    possible_module_cols <- intersect(
      c("module_df", "color", "module_color"),
      colnames(modules_df)
    )
    
    if (length(possible_module_cols) > 0) {
      modules_df <- modules_df |>
        rename(module = all_of(possible_module_cols[1]))
    } else {
      stop("Could not identify module column in hdWGCNA module table.")
    }
  }
  
  modules_df
}


extract_gene_module_table <- function(obj, module_stats, cfg) {
  modules_df <- GetModules(
    seurat_obj = obj,
    wgcna_name = cfg$wgcna_name
  ) |>
    as_tibble() |>
    standardise_module_table()
  
  # Try to retrieve hdWGCNA degree/connectivity table.
  deg <- tryCatch(
    {
      obj@misc[[cfg$wgcna_name]]$wgcna_degrees |>
        as.data.frame() |>
        rownames_to_column("gene_name")
    },
    error = function(e) {
      tibble(gene_name = character())
    }
  )
  
  # Add module stats to each gene.
  final_tbl <- modules_df |>
    left_join(module_stats, by = "module")
  
  # Join degree table where possible.
  if (nrow(deg) > 0 && "gene_name" %in% colnames(deg)) {
    common_join_cols <- intersect(c("gene_name", "module"), colnames(deg))
    
    if ("module" %in% common_join_cols) {
      final_tbl <- final_tbl |>
        left_join(deg, by = c("gene_name", "module"))
    } else {
      final_tbl <- final_tbl |>
        left_join(deg, by = "gene_name")
    }
  }
  
  final_tbl |>
    arrange(module, gene_name)
}


read_edgeR_gene_set <- function(cfg) {
  if (isTRUE(cfg$use_edgeR_de_csv) && file.exists(cfg$edgeR_de_csv)) {
    de <- read_csv(cfg$edgeR_de_csv, show_col_types = FALSE)
    
    required_cols <- c(
      cfg$edgeR_gene_col,
      cfg$edgeR_fdr_col,
      cfg$edgeR_logfc_col
    )
    
    missing_cols <- setdiff(required_cols, colnames(de))
    
    if (length(missing_cols) > 0) {
      warning(
        "edgeR CSV missing required column(s): ",
        paste(missing_cols, collapse = ", "),
        ". Falling back to hard-coded gene set."
      )
    } else {
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
        pull(.data[[cfg$edgeR_gene_col]]) |>
        unique()
      
      if (length(genes) > 0) {
        return(tibble(
          gene_name = genes,
          source = "edgeR_csv",
          target_cluster = cfg$target_cluster
        ))
      }
      
      warning("No edgeR genes passed configured thresholds. Falling back to hard-coded gene set.")
    }
  } else {
    warning("edgeR CSV not found. Falling back to hard-coded gene set.")
  }
  
  tibble(
    gene_name = unique(cfg$fallback_gene_set),
    source = "fallback_gene_set",
    target_cluster = cfg$target_cluster
  )
}


calculate_module_overlap <- function(gene_module_table, edgeR_genes_df) {
  modules_df <- gene_module_table |>
    filter(!is.na(gene_name), !is.na(module)) |>
    distinct(gene_name, module)
  
  universe <- unique(modules_df$gene_name)
  edgeR_genes <- intersect(unique(edgeR_genes_df$gene_name), universe)
  
  module_sizes <- modules_df |>
    count(module, name = "module_size")
  
  overlap_tbl <- modules_df |>
    mutate(is_edgeR_gene = gene_name %in% edgeR_genes) |>
    group_by(module) |>
    summarise(
      module_size = n_distinct(gene_name),
      n_edgeR_genes_in_module = sum(is_edgeR_gene),
      edgeR_genes_in_module = paste(sort(unique(gene_name[is_edgeR_gene])), collapse = ";"),
      .groups = "drop"
    ) |>
    mutate(
      n_edgeR_genes_total = length(edgeR_genes),
      universe_size = length(universe),
      overlap_fraction_of_module = n_edgeR_genes_in_module / module_size,
      overlap_fraction_of_edgeR = ifelse(
        n_edgeR_genes_total > 0,
        n_edgeR_genes_in_module / n_edgeR_genes_total,
        NA_real_
      )
    )
  
  # Fisher enrichment per module.
  overlap_tbl <- overlap_tbl |>
    rowwise() |>
    mutate(
      fisher_p = {
        a <- n_edgeR_genes_in_module
        b <- module_size - a
        c <- n_edgeR_genes_total - a
        d <- universe_size - a - b - c
        
        if (any(c(a, b, c, d) < 0) || n_edgeR_genes_total == 0) {
          NA_real_
        } else {
          fisher.test(
            matrix(c(a, b, c, d), nrow = 2),
            alternative = "greater"
          )$p.value
        }
      }
    ) |>
    ungroup() |>
    mutate(
      fisher_p_adj = p.adjust(fisher_p, method = "BH")
    ) |>
    arrange(fisher_p_adj, desc(n_edgeR_genes_in_module))
  
  shared_genes <- modules_df |>
    filter(gene_name %in% edgeR_genes) |>
    arrange(module, gene_name)
  
  list(
    overlap = overlap_tbl,
    shared_genes = shared_genes
  )
}


# ----------------------------- #
# 3. Run pipeline
# ----------------------------- #

setup_dirs(config)

message("Loading activated EBV object...")
merged <- load_seurat_object(config$input_rds)

message("Checking required metadata...")
check_required_metadata(merged, config)

message("Defining latent_group...")
merged <- add_latent_group(merged, config)

message("Subsetting stimulated latent target cluster...")
stim_lat_cluster <- subset_stimulated_latent_cluster(merged, config)

message("Preprocessing target cells for hdWGCNA...")
stim_lat_cluster <- preprocess_for_hdWGCNA(stim_lat_cluster, config)

message("Setting up hdWGCNA metacells and expression matrix...")
stim_lat_cluster <- setup_hdWGCNA(stim_lat_cluster, config)

message("Testing soft powers and constructing network...")
stim_lat_cluster <- run_network_construction(stim_lat_cluster, config)

message("Calculating module eigengenes and connectivity...")
stim_lat_cluster <- run_module_eigengenes_and_connectivity(stim_lat_cluster, config)

message("Extracting module eigengenes...")
MEs <- get_module_eigengenes_safe(stim_lat_cluster, config)

ME_tables <- summarise_module_eigengenes(
  obj = stim_lat_cluster,
  MEs = MEs,
  cfg = config
)

write.csv(
  ME_tables$cell,
  config$module_eigengene_cell_csv,
  row.names = FALSE
)

write.csv(
  ME_tables$sample,
  config$module_eigengene_sample_csv,
  row.names = FALSE
)

message("Testing module eigengene differences between MS and Control...")
module_stats <- wilcox_module_eigengenes(ME_tables$sample)

write.csv(
  module_stats,
  config$module_eigengene_stats_csv,
  row.names = FALSE
)

message("Extracting gene module assignments and connectivity...")
gene_module_table <- extract_gene_module_table(
  obj = stim_lat_cluster,
  module_stats = module_stats,
  cfg = config
)

write.csv(
  gene_module_table,
  config$gene_module_table_csv,
  row.names = FALSE
)

message("Reading edgeR DEG set for overlap analysis...")
edgeR_gene_set <- read_edgeR_gene_set(config)

write.csv(
  edgeR_gene_set,
  config$edgeR_gene_set_csv,
  row.names = FALSE
)

message("Calculating DEG overlap across WGCNA modules...")
overlap_results <- calculate_module_overlap(
  gene_module_table = gene_module_table,
  edgeR_genes_df = edgeR_gene_set
)

write.csv(
  overlap_results$overlap,
  config$module_overlap_csv,
  row.names = FALSE
)

write.csv(
  overlap_results$shared_genes,
  config$shared_genes_csv,
  row.names = FALSE
)

message("Saving hdWGCNA object...")
saveRDS(
  stim_lat_cluster,
  config$output_wgcna_rds
)

writeLines(
  capture.output(sessionInfo()),
  config$session_info_file
)

message("\nhdWGCNA stimulated latent cluster module analysis complete.")
message("Target matched_cluster: ", config$target_cluster)
message("Saved hdWGCNA object to: ", config$output_wgcna_rds)
message("Saved soft-power plot to: ", config$soft_power_plot_png)
message("Saved dendrogram to: ", config$dendrogram_png)
message("Saved module eigengene stats to: ", config$module_eigengene_stats_csv)
message("Saved gene module table to: ", config$gene_module_table_csv)
message("Saved edgeR overlap table to: ", config$module_overlap_csv)
message("Saved shared gene list to: ", config$shared_genes_csv)
message("Saved session info to: ", config$session_info_file)