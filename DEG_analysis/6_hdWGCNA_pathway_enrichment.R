#!/usr/bin/env Rscript

# ==============================================================================
# hdWGCNA analysis of stimulated latent EBV-specific T cells
# ==============================================================================
#
# Purpose:
#   This script identifies weighted co-expression modules in stimulated latent
#   EBV-specific T cells from the activated/baseline integrated Seurat object.
#
#   The default analysis focuses on:
#     - Stimulated Latent cells
#     - new_cluster == "2"
#     - MS and Control participants
#
#   It then:
#     1. Runs hdWGCNA on the selected cells.
#     2. Builds metacells by sample.
#     3. Tests soft powers and constructs a signed co-expression network.
#     4. Calculates module eigengenes.
#     5. Compares module eigengene activity between MS and Control at the
#        sample level.
#     6. Exports genes from a selected module, default turquoise.
#     7. Optionally overlaps module genes with an external DE table.
#     8. Runs Enrichr pathway enrichment.
#     9. Saves enrichment tables and publication-style dot plots.
#
# Input:
#   A Seurat RDS object, for example:
#     3_activated_EBV_module_scored_reannotated.rds
#
# Required metadata columns:
#   - batch
#   - lifecycle
#   - cohort
#   - sample
#   - new_cluster
#
# Notes:
#   - Paths are intentionally generic for GitHub.
#   - Update the config block before running.
#   - Proxy credentials are NOT stored in this script. If Enrichr access requires
#     a proxy, use environment variables instead.
#   - Enrichr requires internet access.
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
  library(ggplot2)
  library(viridis)
  library(hdWGCNA)
  library(enrichR)
  library(GeneOverlap)
  library(readr)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  
  # Input RDS
  input_rds = "path/to/Publication_RDS/Activated_EBV_integrated_seurat.rds",
  
  # Output directory
  output_dir = "path/to/Publication_data/EBV_activated/hdWGCNA_stimulated_latent_cluster2",
  
  # Metadata columns
  batch_col = "batch",
  lifecycle_col = "lifecycle",
  cohort_col = "cohort",
  sample_col = "sample",
  cluster_col = "new_cluster",
  
  # Group definitions
  baseline_batch_pattern = "^GEMEBV",
  stimulated_batch_pattern = "^EBVLCL",
  latent_lifecycle_label = "Latent",
  baseline_label = "Baseline Latent",
  stimulated_label = "Stimulated Latent",
  non_latent_label = "Non-Latent",
  
  # Analysis subset
  target_latent_group = "Stimulated Latent",
  target_cluster = "2",
  cohorts_to_keep = c("MS", "Control"),
  
  # Seurat processing
  assay = "RNA",
  reduction = "pca",
  n_variable_features = 5000,
  n_pcs = 30,
  
  # hdWGCNA settings
  wgcna_name = "stimLat_cluster2",
  gene_select = "fraction",
  fraction = 0.05,
  metacell_group_by = "sample",
  metacell_ident_group = "sample",
  metacell_k = 3,
  metacell_max_shared = 10,
  metacell_min_cells = 3,
  network_type = "signed",
  soft_power = 9,
  tom_name = "stimLat_cluster2_TOM",
  set_dat_expr_layer = "data",
  
  # Module to export/enrich
  target_module = "turquoise",
  
  # Module eigengene statistics
  min_samples_per_cohort_for_ME_test = 2,
  
  # Optional overlap with an external DE table.
  # Set de_results_csv to a valid CSV path to enable.
  # The DE table should contain at least a gene column.
  de_results_csv = NA_character_,
  de_gene_col = "gene",
  de_fdr_col = "p_val_adj",
  de_logfc_col = "avg_log2FC",
  de_cluster_col = "cluster",
  de_target_cluster = "2",
  de_fdr_threshold = 0.05,
  de_abs_logfc_threshold = 0,
  
  # Enrichr settings
  run_enrichr = TRUE,
  enrichr_dbs = c(
    "GO_Biological_Process_2025",
    "KEGG_2026",
    "MSigDB_Hallmark_2020",
    "Reactome_Pathways_2024"
  ),
  enrichr_plot_module = "turquoise",
  enrichr_dotplot_database = "GO_Biological_Process_2025",
  enrichr_dotplot_n_terms = 10,
  
  # Custom pathway plot settings
  reactome_database = "Reactome_Pathways_2024",
  combined_plot_databases = c("GO_Biological_Process_2025", "Reactome_Pathways_2024"),
  adjusted_p_threshold_for_pathway_plot = 0.05,
  n_reactome_terms = 10,
  n_combined_terms = 20,
  
  # Optional proxy configuration for Enrichr access.
  # Do NOT hard-code credentials. Set these environment variables instead:
  #   HTTPS_PROXY_HOST, HTTPS_PROXY_PORT, HTTPS_PROXY_USER, HTTPS_PROXY_PASS
  use_proxy_from_env = FALSE,
  
  # Output files
  subset_metadata_csv = "stimulated_latent_cluster_subset_metadata.csv",
  soft_power_plot_pdf = "hdWGCNA_soft_power_plots.pdf",
  module_eigengenes_cells_csv = "module_eigengenes_cell_level.csv",
  module_eigengenes_sample_csv = "module_eigengenes_sample_level.csv",
  module_eigengene_wilcox_csv = "module_eigengene_MS_vs_Control_wilcox_BH.csv",
  modules_all_csv = "hdWGCNA_modules_all_genes.csv",
  target_module_genes_csv = "hdWGCNA_target_module_genes.csv",
  target_module_genes_txt = "hdWGCNA_target_module_genes.txt",
  de_overlap_csv = "hdWGCNA_target_module_overlap_with_DE_genes.csv",
  enrichr_table_csv = "hdWGCNA_enrichr_all_modules.csv",
  enrichr_dotplot_png = "hdWGCNA_target_module_enrichr_dotplot.png",
  reactome_plot_png = "Reactome_pathway_analysis_target_module.png",
  combined_pathway_plot_png = "GO_Reactome_pathway_analysis_target_module.png",
  session_info_file = "sessionInfo_hdWGCNA_stimulated_latent_cluster.txt",
  
  # Plot sizes
  enrichr_dotplot_width = 9,
  enrichr_dotplot_height = 7,
  pathway_plot_width = 12,
  pathway_plot_height = 15,
  plot_dpi = 400
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


clean_chr <- function(x) {
  
  x <- str_trim(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL", "null", "None", "none")] <- NA_character_
  x
}


safe_name <- function(x) {
  
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+$", "") %>%
    str_replace_all("^_+", "")
}


configure_proxy_from_env <- function(use_proxy = FALSE) {
  
  if (!isTRUE(use_proxy)) {
    return(invisible(FALSE))
  }
  
  proxy_host <- Sys.getenv("HTTPS_PROXY_HOST", unset = NA_character_)
  proxy_port <- Sys.getenv("HTTPS_PROXY_PORT", unset = NA_character_)
  proxy_user <- Sys.getenv("HTTPS_PROXY_USER", unset = NA_character_)
  proxy_pass <- Sys.getenv("HTTPS_PROXY_PASS", unset = NA_character_)
  
  if (is.na(proxy_host) || is.na(proxy_port)) {
    warning(
      "Proxy requested, but HTTPS_PROXY_HOST or HTTPS_PROXY_PORT is missing. ",
      "Continuing without proxy configuration."
    )
    return(invisible(FALSE))
  }
  
  rcurl_options <- list(
    proxy = proxy_host,
    proxyport = as.integer(proxy_port),
    proxyauth = "basic"
  )
  
  if (!is.na(proxy_user) && proxy_user != "") {
    rcurl_options$proxyusername <- proxy_user
  }
  
  if (!is.na(proxy_pass) && proxy_pass != "") {
    rcurl_options$proxypassword <- proxy_pass
  }
  
  options(RCurlOptions = rcurl_options)
  
  invisible(TRUE)
}


get_logfc_col <- function(df) {
  
  possible_cols <- c("avg_log2FC", "avg_logFC", "logFC")
  hit <- intersect(possible_cols, colnames(df))
  
  if (length(hit) == 0) {
    return(NA_character_)
  }
  
  hit[1]
}


run_safe_wilcox <- function(values, groups, min_n_per_group = 2) {
  
  test_df <- tibble(
    value = as.numeric(values),
    group = as.character(groups)
  ) %>%
    filter(!is.na(value), !is.na(group))
  
  if (n_distinct(test_df$group) < 2) {
    return(NA_real_)
  }
  
  group_counts <- table(test_df$group)
  
  if (any(group_counts < min_n_per_group)) {
    return(NA_real_)
  }
  
  if (length(unique(test_df$value)) < 2) {
    return(1)
  }
  
  tryCatch(
    wilcox.test(value ~ group, data = test_df, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}


make_pathway_plot <- function(df_plot, x_label = "log10(Enrichment Score)") {
  
  ggplot(
    df_plot,
    aes(
      x = log_combined,
      y = reorder(Term, log_combined)
    )
  ) +
    geom_point(
      aes(size = gene_count, fill = log10_fdr),
      shape = 21,
      color = "black",
      stroke = 0.8,
      alpha = 0.95
    ) +
    scale_size(range = c(6, 16)) +
    scale_fill_gradientn(
      colours = c("#f0d3d3", "#e0a7a7", "#d17a7a", "#c14e4e", "#b22222")
    ) +
    theme_classic() +
    labs(
      x = x_label,
      y = "",
      fill = "-log10(FDR)",
      size = "Gene count"
    ) +
    theme(
      axis.text.y = element_text(size = 13),
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
}


# ----------------------------- #
# 3. Configure proxy if requested
# ----------------------------- #

configure_proxy_from_env(config$use_proxy_from_env)


# ----------------------------- #
# 4. Load object and define latent group
# ----------------------------- #

create_dir(config$output_dir)

message("Loading activated EBV object:")
message(config$input_rds)

merged <- load_seurat_object(config$input_rds)

required_cols <- c(
  config$batch_col,
  config$lifecycle_col,
  config$cohort_col,
  config$sample_col,
  config$cluster_col
)

check_required_columns(
  df = merged@meta.data,
  required_cols = required_cols,
  object_name = "activated EBV metadata"
)

if (!config$assay %in% Assays(merged)) {
  stop(
    "Assay not found: ",
    config$assay,
    ". Available assays: ",
    paste(Assays(merged), collapse = ", ")
  )
}

meta_clean <- merged@meta.data %>%
  as.data.frame() %>%
  mutate(
    batch_for_wgcna = clean_chr(.data[[config$batch_col]]),
    lifecycle_for_wgcna = clean_chr(.data[[config$lifecycle_col]]),
    cohort_for_wgcna = clean_chr(.data[[config$cohort_col]]),
    sample_for_wgcna = clean_chr(.data[[config$sample_col]]),
    cluster_for_wgcna = clean_chr(.data[[config$cluster_col]]),
    
    latent_group = case_when(
      !is.na(batch_for_wgcna) &
        grepl(config$baseline_batch_pattern, batch_for_wgcna) &
        lifecycle_for_wgcna == config$latent_lifecycle_label ~ config$baseline_label,
      
      !is.na(batch_for_wgcna) &
        grepl(config$stimulated_batch_pattern, batch_for_wgcna) &
        lifecycle_for_wgcna == config$latent_lifecycle_label ~ config$stimulated_label,
      
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

rownames(meta_clean) <- colnames(merged)
merged@meta.data <- meta_clean

message("Latent group counts:")
print(table(merged$latent_group, useNA = "ifany"))


# ----------------------------- #
# 5. Subset target cells
# ----------------------------- #

cells_keep <- rownames(merged@meta.data)[
  merged$latent_group == config$target_latent_group &
    merged$cohort_for_wgcna %in% config$cohorts_to_keep &
    merged$cluster_for_wgcna %in% config$target_cluster
]

if (length(cells_keep) == 0) {
  stop("No cells matched target latent group / cohort / cluster filters.")
}

cl_obj <- subset(
  merged,
  cells = cells_keep
)

DefaultAssay(cl_obj) <- config$assay

# Keep standard names expected by hdWGCNA workflow.
cl_obj$sample <- cl_obj$sample_for_wgcna
cl_obj$cohort <- cl_obj$cohort_for_wgcna
cl_obj$cluster <- cl_obj$cluster_for_wgcna

subset_metadata_path <- file.path(
  config$output_dir,
  config$subset_metadata_csv
)

write_csv(
  cl_obj@meta.data %>%
    as.data.frame() %>%
    rownames_to_column("cell_barcode"),
  subset_metadata_path
)

message("Cells in hdWGCNA subset:")
print(ncol(cl_obj))

message("Cells by cohort:")
print(table(cl_obj$cohort, useNA = "ifany"))

message("Cells by sample:")
print(table(cl_obj$sample, useNA = "ifany"))


# ----------------------------- #
# 6. Standard preprocessing for hdWGCNA
# ----------------------------- #

# JoinLayers is needed for Seurat v5 objects with split layers.
cl_obj <- tryCatch(
  {
    JoinLayers(cl_obj, assay = config$assay)
  },
  error = function(e) {
    message("JoinLayers skipped or failed: ", conditionMessage(e))
    cl_obj
  }
)

cl_obj <- NormalizeData(cl_obj, verbose = FALSE)
cl_obj <- FindVariableFeatures(
  cl_obj,
  nfeatures = config$n_variable_features,
  verbose = FALSE
)

cl_obj <- ScaleData(
  cl_obj,
  features = VariableFeatures(cl_obj),
  verbose = FALSE
)

cl_obj <- RunPCA(
  cl_obj,
  npcs = config$n_pcs,
  verbose = FALSE
)


# ----------------------------- #
# 7. hdWGCNA setup and metacells
# ----------------------------- #

cl_obj <- SetupForWGCNA(
  cl_obj,
  gene_select = config$gene_select,
  fraction = config$fraction,
  wgcna_name = config$wgcna_name
)

cl_obj <- MetacellsByGroups(
  seurat_obj = cl_obj,
  group.by = config$metacell_group_by,
  reduction = config$reduction,
  k = config$metacell_k,
  max_shared = config$metacell_max_shared,
  min_cells = config$metacell_min_cells,
  ident.group = config$metacell_ident_group
)

cl_obj <- NormalizeMetacells(cl_obj)

metacell_obj <- GetMetacellObject(cl_obj)
metacell_samples <- sort(unique(metacell_obj$sample))

message("Metacell samples:")
print(metacell_samples)

cl_obj <- SetDatExpr(
  cl_obj,
  group_name = metacell_samples,
  group.by = "sample",
  assay = config$assay,
  layer = config$set_dat_expr_layer
)


# ----------------------------- #
# 8. Soft power testing and network construction
# ----------------------------- #

cl_obj <- TestSoftPowers(
  cl_obj,
  networkType = config$network_type
)

soft_power_plots <- PlotSoftPowers(cl_obj)

soft_power_plot_path <- file.path(
  config$output_dir,
  config$soft_power_plot_pdf
)

pdf(soft_power_plot_path, width = 12, height = 8, useDingbats = FALSE)

if (is.list(soft_power_plots)) {
  for (p in soft_power_plots) {
    print(p)
  }
} else {
  print(soft_power_plots)
}

dev.off()

message("Using soft power: ", config$soft_power)

cl_obj <- ConstructNetwork(
  cl_obj,
  soft_power = config$soft_power,
  tom_name = config$tom_name
)


# ----------------------------- #
# 9. Module eigengenes and module testing
# ----------------------------- #

cl_obj <- ModuleEigengenes(
  cl_obj,
  group.by.vars = "sample",
  wgcna_name = config$wgcna_name
)

MEs <- GetMEs(
  cl_obj,
  wgcna_name = config$wgcna_name,
  harmonized = TRUE
)

MEs <- MEs[colnames(cl_obj), , drop = FALSE]

meta_me <- cl_obj@meta.data[colnames(cl_obj), c("sample", "cohort"), drop = FALSE]

ME_cells <- as.data.frame(MEs) %>%
  rownames_to_column("cell_barcode") %>%
  mutate(
    sample = meta_me$sample,
    cohort = meta_me$cohort
  )

module_eigengenes_cells_path <- file.path(
  config$output_dir,
  config$module_eigengenes_cells_csv
)

write_csv(ME_cells, module_eigengenes_cells_path)

ME_sample <- as.data.frame(MEs) %>%
  mutate(
    sample = meta_me$sample,
    cohort = meta_me$cohort
  ) %>%
  group_by(sample, cohort) %>%
  summarise(
    across(where(is.numeric), mean),
    .groups = "drop"
  )

module_eigengenes_sample_path <- file.path(
  config$output_dir,
  config$module_eigengenes_sample_csv
)

write_csv(ME_sample, module_eigengenes_sample_path)

module_cols <- setdiff(colnames(ME_sample), c("sample", "cohort"))

ME_stats <- bind_rows(lapply(module_cols, function(module_name) {
  
  p_val <- run_safe_wilcox(
    values = ME_sample[[module_name]],
    groups = ME_sample$cohort,
    min_n_per_group = config$min_samples_per_cohort_for_ME_test
  )
  
  tibble(
    module = module_name,
    n_MS = sum(ME_sample$cohort == "MS"),
    n_Control = sum(ME_sample$cohort == "Control"),
    median_MS = median(ME_sample[[module_name]][ME_sample$cohort == "MS"], na.rm = TRUE),
    median_Control = median(ME_sample[[module_name]][ME_sample$cohort == "Control"], na.rm = TRUE),
    diff_median_MS_minus_Control = median_MS - median_Control,
    p_value = p_val
  )
})) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_adj, p_value)

module_eigengene_wilcox_path <- file.path(
  config$output_dir,
  config$module_eigengene_wilcox_csv
)

write_csv(ME_stats, module_eigengene_wilcox_path)


# ----------------------------- #
# 10. Module gene export
# ----------------------------- #

modules_tbl <- GetModules(cl_obj)

modules_all_path <- file.path(
  config$output_dir,
  config$modules_all_csv
)

write_csv(modules_tbl, modules_all_path)

target_module_genes <- modules_tbl %>%
  filter(module == config$target_module) %>%
  pull(gene_name) %>%
  unique() %>%
  sort()

target_module_genes_tbl <- tibble(
  module = config$target_module,
  gene = target_module_genes
)

target_module_genes_path <- file.path(
  config$output_dir,
  config$target_module_genes_csv
)

target_module_genes_txt_path <- file.path(
  config$output_dir,
  config$target_module_genes_txt
)

write_csv(target_module_genes_tbl, target_module_genes_path)
writeLines(target_module_genes, target_module_genes_txt_path)

message("Number of genes in target module ", config$target_module, ": ", length(target_module_genes))


# ----------------------------- #
# 11. Optional overlap with DE genes
# ----------------------------- #

de_overlap_tbl <- tibble()

if (!is.na(config$de_results_csv) && file.exists(config$de_results_csv)) {
  
  de_tbl <- read_csv(config$de_results_csv, show_col_types = FALSE)
  
  if (!config$de_gene_col %in% colnames(de_tbl)) {
    stop("DE gene column not found: ", config$de_gene_col)
  }
  
  logfc_col <- if (config$de_logfc_col %in% colnames(de_tbl)) {
    config$de_logfc_col
  } else {
    get_logfc_col(de_tbl)
  }
  
  de_filtered <- de_tbl
  
  if (config$de_cluster_col %in% colnames(de_filtered)) {
    de_filtered <- de_filtered %>%
      filter(as.character(.data[[config$de_cluster_col]]) == as.character(config$de_target_cluster))
  }
  
  if (config$de_fdr_col %in% colnames(de_filtered)) {
    de_filtered <- de_filtered %>%
      filter(!is.na(.data[[config$de_fdr_col]]), .data[[config$de_fdr_col]] < config$de_fdr_threshold)
  }
  
  if (!is.na(logfc_col) && logfc_col %in% colnames(de_filtered)) {
    de_filtered <- de_filtered %>%
      filter(!is.na(.data[[logfc_col]]), abs(.data[[logfc_col]]) > config$de_abs_logfc_threshold)
  }
  
  de_genes <- de_filtered %>%
    pull(.data[[config$de_gene_col]]) %>%
    as.character() %>%
    unique()
  
  overlap_genes <- intersect(target_module_genes, de_genes)
  
  universe_genes <- unique(c(modules_tbl$gene_name, de_tbl[[config$de_gene_col]]))
  
  overlap_test <- tryCatch(
    {
      GeneOverlap::newGeneOverlap(
        listA = target_module_genes,
        listB = de_genes,
        genome.size = length(universe_genes)
      ) %>%
        GeneOverlap::testGeneOverlap()
    },
    error = function(e) NULL
  )
  
  de_overlap_tbl <- tibble(
    module = config$target_module,
    de_results_csv = config$de_results_csv,
    n_module_genes = length(target_module_genes),
    n_de_genes = length(de_genes),
    n_overlap_genes = length(overlap_genes),
    overlap_genes = paste(overlap_genes, collapse = ";"),
    odds_ratio = if (!is.null(overlap_test)) overlap_test@odds.ratio else NA_real_,
    p_value = if (!is.null(overlap_test)) overlap_test@pval else NA_real_
  )
  
} else {
  
  de_overlap_tbl <- tibble(
    module = config$target_module,
    de_results_csv = config$de_results_csv,
    n_module_genes = length(target_module_genes),
    n_de_genes = NA_integer_,
    n_overlap_genes = NA_integer_,
    overlap_genes = NA_character_,
    odds_ratio = NA_real_,
    p_value = NA_real_,
    note = "No DE results CSV supplied or file does not exist."
  )
}

de_overlap_path <- file.path(
  config$output_dir,
  config$de_overlap_csv
)

write_csv(de_overlap_tbl, de_overlap_path)


# ----------------------------- #
# 12. Enrichr pathway enrichment
# ----------------------------- #

enrichr_table <- tibble()

if (isTRUE(config$run_enrichr)) {
  
  message("Checking Enrichr databases...")
  
  dbs_available <- tryCatch(
    {
      enrichR::listEnrichrDbs()
    },
    error = function(e) {
      warning("Could not retrieve Enrichr databases: ", conditionMessage(e))
      NULL
    }
  )
  
  if (!is.null(dbs_available)) {
    missing_dbs <- setdiff(config$enrichr_dbs, dbs_available$libraryName)
    
    if (length(missing_dbs) > 0) {
      warning(
        "These Enrichr databases were not found in listEnrichrDbs(): ",
        paste(missing_dbs, collapse = ", "),
        ". RunEnrichr may fail if database names have changed."
      )
    }
  }
  
  cl_obj <- ModuleConnectivity(
    cl_obj,
    group.by = "sample"
  )
  
  cl_obj <- RunEnrichr(
    cl_obj,
    dbs = config$enrichr_dbs,
    max_genes = Inf
  )
  
  enrichr_table <- GetEnrichrTable(cl_obj)
  
  enrichr_table_path <- file.path(
    config$output_dir,
    config$enrichr_table_csv
  )
  
  write_csv(enrichr_table, enrichr_table_path)
  
  # hdWGCNA Enrichr dotplot.
  enrichr_dotplot <- tryCatch(
    {
      EnrichrDotPlot(
        cl_obj,
        mods = config$enrichr_plot_module,
        database = config$enrichr_dotplot_database,
        n_terms = config$enrichr_dotplot_n_terms,
        term_size = 8,
        p_adj = TRUE
      ) +
        scale_color_stepsn(colors = rev(viridis::magma(256)))
    },
    error = function(e) {
      warning("Could not create EnrichrDotPlot: ", conditionMessage(e))
      NULL
    }
  )
  
  enrichr_dotplot_path <- file.path(
    config$output_dir,
    config$enrichr_dotplot_png
  )
  
  if (!is.null(enrichr_dotplot)) {
    ggsave(
      filename = enrichr_dotplot_path,
      plot = enrichr_dotplot,
      width = config$enrichr_dotplot_width,
      height = config$enrichr_dotplot_height,
      dpi = config$plot_dpi
    )
  }
  
  # Reactome-only plot for target module.
  reactome_plot_path <- file.path(
    config$output_dir,
    config$reactome_plot_png
  )
  
  reactome_df <- enrichr_table %>%
    filter(
      module == config$target_module,
      db == config$reactome_database,
      Adjusted.P.value < config$adjusted_p_threshold_for_pathway_plot
    ) %>%
    arrange(desc(Combined.Score)) %>%
    slice_head(n = config$n_reactome_terms) %>%
    mutate(
      gene_count = as.numeric(sub("/.*", "", Overlap)),
      log10_fdr = -log10(Adjusted.P.value),
      log_combined = log10(Combined.Score),
      Term = str_wrap(Term, width = 45)
    )
  
  if (nrow(reactome_df) > 0) {
    p_reactome <- make_pathway_plot(
      reactome_df,
      x_label = "log10(Enrichment Score)"
    )
    
    ggsave(
      filename = reactome_plot_path,
      plot = p_reactome,
      width = config$pathway_plot_width,
      height = config$pathway_plot_height,
      dpi = config$plot_dpi
    )
  } else {
    warning("No significant Reactome terms available for the Reactome plot.")
  }
  
  # Combined GO BP + Reactome plot.
  combined_pathway_plot_path <- file.path(
    config$output_dir,
    config$combined_pathway_plot_png
  )
  
  combined_df <- enrichr_table %>%
    filter(
      module == config$target_module,
      db %in% config$combined_plot_databases,
      Adjusted.P.value < config$adjusted_p_threshold_for_pathway_plot
    ) %>%
    mutate(
      source = case_when(
        db == "GO_Biological_Process_2025" ~ "GO",
        db == "Reactome_Pathways_2024" ~ "Reactome",
        TRUE ~ db
      )
    ) %>%
    arrange(desc(Combined.Score)) %>%
    slice_head(n = config$n_combined_terms) %>%
    mutate(
      gene_count = as.numeric(sub("/.*", "", Overlap)),
      log10_fdr = -log10(Adjusted.P.value),
      log_combined = log10(Combined.Score),
      Term = str_wrap(paste0(Term, " (", source, ")"), width = 45)
    )
  
  if (nrow(combined_df) > 0) {
    p_combined <- make_pathway_plot(
      combined_df,
      x_label = "log10(Enrichment Score)"
    )
    
    ggsave(
      filename = combined_pathway_plot_path,
      plot = p_combined,
      width = config$pathway_plot_width,
      height = config$pathway_plot_height,
      dpi = config$plot_dpi
    )
  } else {
    warning("No significant combined GO/Reactome terms available for the combined plot.")
  }
  
} else {
  
  enrichr_table_path <- file.path(
    config$output_dir,
    config$enrichr_table_csv
  )
  
  enrichr_dotplot_path <- file.path(
    config$output_dir,
    config$enrichr_dotplot_png
  )
  
  reactome_plot_path <- file.path(
    config$output_dir,
    config$reactome_plot_png
  )
  
  combined_pathway_plot_path <- file.path(
    config$output_dir,
    config$combined_pathway_plot_png
  )
  
  write_csv(enrichr_table, enrichr_table_path)
}


# ----------------------------- #
# 13. Save processed hdWGCNA object
# ----------------------------- #

processed_object_path <- file.path(
  config$output_dir,
  paste0("hdWGCNA_", safe_name(config$wgcna_name), "_processed.rds")
)

saveRDS(cl_obj, processed_object_path)


# ----------------------------- #
# 14. Session information
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
    capture.output(print(table(merged$latent_group, useNA = "ifany"))),
    "",
    "Subset cell count:",
    as.character(ncol(cl_obj)),
    "",
    "Subset cohort counts:",
    capture.output(print(table(cl_obj$cohort, useNA = "ifany"))),
    "",
    "Subset sample counts:",
    capture.output(print(table(cl_obj$sample, useNA = "ifany"))),
    "",
    "Soft power used:",
    as.character(config$soft_power),
    "",
    "Target module:",
    config$target_module,
    "",
    "Number of target module genes:",
    as.character(length(target_module_genes)),
    "",
    "Module eigengene statistics:",
    capture.output(print(ME_stats)),
    "",
    "DE overlap:",
    capture.output(print(de_overlap_tbl)),
    "",
    "Output files:",
    subset_metadata_path,
    soft_power_plot_path,
    module_eigengenes_cells_path,
    module_eigengenes_sample_path,
    module_eigengene_wilcox_path,
    modules_all_path,
    target_module_genes_path,
    target_module_genes_txt_path,
    de_overlap_path,
    enrichr_table_path,
    enrichr_dotplot_path,
    reactome_plot_path,
    combined_pathway_plot_path,
    processed_object_path,
    "",
    "Session information:",
    capture.output(sessionInfo())
  ),
  session_info_path
)


# ----------------------------- #
# 15. Completion messages
# ----------------------------- #

message("\nhdWGCNA stimulated latent cluster workflow complete.")
message("Subset metadata: ", subset_metadata_path)
message("Soft-power plot: ", soft_power_plot_path)
message("Module eigengenes, cell level: ", module_eigengenes_cells_path)
message("Module eigengenes, sample level: ", module_eigengenes_sample_path)
message("Module eigengene stats: ", module_eigengene_wilcox_path)
message("All module genes: ", modules_all_path)
message("Target module genes: ", target_module_genes_path)
message("DE overlap: ", de_overlap_path)
message("Enrichr table: ", enrichr_table_path)
message("Enrichr dotplot: ", enrichr_dotplot_path)
message("Reactome plot: ", reactome_plot_path)
message("GO/Reactome combined plot: ", combined_pathway_plot_path)
message("Processed hdWGCNA object: ", processed_object_path)
message("Session info: ", session_info_path)
