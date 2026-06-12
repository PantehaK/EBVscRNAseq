#!/usr/bin/env Rscript

# ==============================================================================
# Pseudobulk edgeR differential expression by cluster
# ==============================================================================
#
# Purpose:
#   This script performs pseudobulk differential expression using edgeR to compare
#   MS versus Control T cells within each UMAP/Harmony cluster.
#
#   It can be run separately for:
#     - Latent EBV-specific T cells
#     - Lytic EBV-specific T cells
#     - virus-negative / non-EBV CD8+ T cells
#
#   For each cluster, cells are aggregated into sample-level pseudobulk counts,
#   then edgeR quasi-likelihood negative binomial modelling is used to test
#   MS versus Control.
#
# Outputs:
#   1. all edgeR results across clusters,
#   2. significant results only,
#   3. cluster-level donor QC table,
#   4. skipped cluster log,
#   5. session information.
#
# Expected input:
#   A module-scored baseline EBV Seurat object, for example:
#     15_baseline_EBV_module_scored.rds
#
# Notes:
#   - This script uses RNA counts for pseudobulk aggregation.
#   - Clusters with fewer than the configured minimum donors per cohort are skipped.
#   - logFC is MS versus Control.
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
  library(edgeR)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  dataset_id = "baseline_ebv",
  
  input_rds = "path/to/input/15_baseline_EBV_module_scored.rds",
  
  output_dir = "path/to/output/Publication_data/EBV_baseline/pseudobulk_edgeR_by_cluster",
  
  all_de_csv = "path/to/output/Publication_data/EBV_baseline/pseudobulk_edgeR_by_cluster/edgeR_Lytic_MS_vs_Control_by_cluster_all.csv",
  sig_de_csv = "path/to/output/Publication_data/EBV_baseline/pseudobulk_edgeR_by_cluster/edgeR_Lytic_MS_vs_Control_by_cluster_sig.csv",
  qc_csv = "path/to/output/Publication_data/EBV_baseline/pseudobulk_edgeR_by_cluster/edgeR_Lytic_MS_vs_Control_by_cluster_QC.csv",
  skipped_clusters_csv = "path/to/output/Publication_data/EBV_baseline/pseudobulk_edgeR_by_cluster/edgeR_Lytic_MS_vs_Control_by_cluster_skipped.csv",
  session_info_file = "path/to/output/Publication_data/EBV_baseline/pseudobulk_edgeR_by_cluster/sessionInfo_edgeR_Lytic_by_cluster.txt",
  
  # Metadata columns.
  sample_col = "sample",
  cohort_col = "cohort",
  cluster_col = "harmony_clusters",
  lifecycle_col = "lifecycle",
  
  # Group comparison.
  cohort_levels = c("Control", "MS"),
  contrast_coef = "cohortMS",
  
  # Change this for different analyses:
  #   "Latent"   = latent EBV-specific T cells
  #   "Lytic"    = lytic EBV-specific T cells
  #   "negative" = non-EBV/background CD8+ T cells
  lifecycle_pattern = "Lytic",
  lifecycle_label = "Lytic",
  
  # Assay settings.
  assay = "RNA",
  aggregate_slot = "counts",
  
  # edgeR settings.
  min_donors_per_cohort = 3,
  robust = TRUE,
  
  # Optional sample-level covariates.
  # Keep empty to match your original code.
  # Example: covariates = c("age")
  covariates = character(),
  
  # Significance threshold.
  fdr_threshold = 0.05,
  
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
  create_parent_dir(cfg$all_de_csv)
  create_parent_dir(cfg$sig_de_csv)
  create_parent_dir(cfg$qc_csv)
  create_parent_dir(cfg$skipped_clusters_csv)
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
    cfg$sample_col,
    cfg$cohort_col,
    cfg$cluster_col,
    cfg$lifecycle_col,
    cfg$covariates
  )
  
  missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
  
  if (length(missing_cols) > 0) {
    stop(
      "Object metadata is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
}


subset_cells_for_analysis <- function(obj, cfg) {
  meta <- obj@meta.data
  
  keep_cells <- rownames(meta)[
    meta[[cfg$cohort_col]] %in% cfg$cohort_levels &
      !is.na(meta[[cfg$lifecycle_col]]) &
      grepl(cfg$lifecycle_pattern, meta[[cfg$lifecycle_col]])
  ]
  
  if (length(keep_cells) == 0) {
    stop("No cells matched lifecycle_pattern: ", cfg$lifecycle_pattern)
  }
  
  obj_sub <- subset(obj, cells = keep_cells)
  
  obj_sub[[cfg$cohort_col]] <- factor(
    obj_sub[[cfg$cohort_col, drop = TRUE]],
    levels = cfg$cohort_levels
  )
  
  obj_sub[[cfg$cluster_col]] <- as.character(
    obj_sub[[cfg$cluster_col, drop = TRUE]]
  )
  
  obj_sub
}


aggregate_counts_by_sample <- function(obj_cluster, cfg) {
  pb <- AggregateExpression(
    object = obj_cluster,
    assays = cfg$assay,
    group.by = cfg$sample_col,
    slot = cfg$aggregate_slot,
    return.seurat = FALSE
  )[[cfg$assay]]
  
  pb
}


make_sample_metadata <- function(obj_cluster, pb, cfg) {
  meta <- obj_cluster@meta.data
  
  smeta <- meta %>%
    distinct(
      sample = .data[[cfg$sample_col]],
      cohort = .data[[cfg$cohort_col]],
      across(any_of(cfg$covariates))
    ) %>%
    filter(sample %in% colnames(pb)) %>%
    arrange(match(sample, colnames(pb)))
  
  if (!all(smeta$sample == colnames(pb))) {
    stop("Sample metadata does not match pseudobulk columns.")
  }
  
  smeta$cohort <- factor(smeta$cohort, levels = cfg$cohort_levels)
  
  smeta
}


count_donors_per_group <- function(obj_cluster, cluster_id, cfg) {
  obj_cluster@meta.data %>%
    distinct(
      sample = .data[[cfg$sample_col]],
      cohort = .data[[cfg$cohort_col]]
    ) %>%
    count(cohort, name = "n_donors") %>%
    complete(
      cohort = cfg$cohort_levels,
      fill = list(n_donors = 0)
    ) %>%
    mutate(
      cluster = cluster_id,
      lifecycle = cfg$lifecycle_label,
      .before = 1
    )
}


make_design_formula <- function(cfg) {
  rhs <- c("cohort", cfg$covariates)
  as.formula(paste("~", paste(rhs, collapse = " + ")))
}


run_edger_one_cluster <- function(obj_cluster, cluster_id, cfg) {
  donors_per_group <- count_donors_per_group(
    obj_cluster = obj_cluster,
    cluster_id = cluster_id,
    cfg = cfg
  )
  
  donors_wide <- donors_per_group %>%
    select(cohort, n_donors) %>%
    deframe()
  
  n_control <- donors_wide[[cfg$cohort_levels[1]]]
  n_ms <- donors_wide[[cfg$cohort_levels[2]]]
  
  if (
    is.null(n_control) || is.null(n_ms) ||
    n_control < cfg$min_donors_per_cohort ||
    n_ms < cfg$min_donors_per_cohort
  ) {
    return(list(
      result = tibble(),
      qc = donors_per_group,
      skipped = tibble(
        lifecycle = cfg$lifecycle_label,
        cluster = cluster_id,
        reason = "insufficient_donors",
        n_control = n_control,
        n_ms = n_ms
      )
    ))
  }
  
  pb <- aggregate_counts_by_sample(
    obj_cluster = obj_cluster,
    cfg = cfg
  )
  
  smeta <- make_sample_metadata(
    obj_cluster = obj_cluster,
    pb = pb,
    cfg = cfg
  )
  
  y <- DGEList(
    counts = pb,
    samples = smeta
  )
  
  y <- calcNormFactors(y)
  
  keep <- filterByExpr(
    y,
    group = y$samples$cohort
  )
  
  if (sum(keep) == 0) {
    return(list(
      result = tibble(),
      qc = donors_per_group,
      skipped = tibble(
        lifecycle = cfg$lifecycle_label,
        cluster = cluster_id,
        reason = "no_genes_passed_filterByExpr",
        n_control = n_control,
        n_ms = n_ms
      )
    ))
  }
  
  y <- y[keep, , keep.lib.sizes = FALSE]
  
  design <- model.matrix(
    make_design_formula(cfg),
    data = y$samples
  )
  
  if (!cfg$contrast_coef %in% colnames(design)) {
    return(list(
      result = tibble(),
      qc = donors_per_group,
      skipped = tibble(
        lifecycle = cfg$lifecycle_label,
        cluster = cluster_id,
        reason = paste0("contrast_coef_not_found: ", cfg$contrast_coef),
        n_control = n_control,
        n_ms = n_ms
      )
    ))
  }
  
  fit_result <- tryCatch(
    {
      y <- estimateDisp(y, design)
      fit <- glmQLFit(y, design, robust = cfg$robust)
      qlf <- glmQLFTest(fit, coef = cfg$contrast_coef)
      edgeR::topTags(qlf, n = Inf)$table
    },
    error = function(e) {
      attr(e, "edgeR_error") <- TRUE
      e
    }
  )
  
  if (inherits(fit_result, "error")) {
    return(list(
      result = tibble(),
      qc = donors_per_group,
      skipped = tibble(
        lifecycle = cfg$lifecycle_label,
        cluster = cluster_id,
        reason = paste0("edgeR_failed: ", fit_result$message),
        n_control = n_control,
        n_ms = n_ms
      )
    ))
  }
  
  tt <- as.data.frame(fit_result) %>%
    rownames_to_column("gene")
  
  if ("PValue" %in% colnames(tt)) {
    tt <- rename(tt, p_val = PValue)
  }
  
  if ("FDR" %in% colnames(tt)) {
    tt <- rename(tt, p_val_adj = FDR)
  }
  
  if ("logFC" %in% colnames(tt)) {
    tt <- mutate(tt, avg_log2FC = logFC)
  } else {
    tt <- mutate(tt, avg_log2FC = NA_real_)
  }
  
  result <- tt %>%
    mutate(
      lifecycle = cfg$lifecycle_label,
      cluster = cluster_id,
      n_donors_MS = n_ms,
      n_donors_Control = n_control
    ) %>%
    relocate(
      lifecycle,
      cluster,
      gene,
      n_donors_MS,
      n_donors_Control,
      avg_log2FC
    )
  
  list(
    result = result,
    qc = donors_per_group,
    skipped = tibble()
  )
}


run_edger_by_cluster <- function(obj, cfg) {
  clusters <- sort(unique(obj[[cfg$cluster_col, drop = TRUE]]))
  
  all_results <- list()
  all_qc <- list()
  skipped <- list()
  
  for (cluster_id in clusters) {
    message("Running cluster: ", cluster_id)
    
    cells_cluster <- rownames(obj@meta.data)[
      obj[[cfg$cluster_col, drop = TRUE]] == cluster_id
    ]
    
    obj_cluster <- subset(obj, cells = cells_cluster)
    
    res <- run_edger_one_cluster(
      obj_cluster = obj_cluster,
      cluster_id = cluster_id,
      cfg = cfg
    )
    
    all_results[[cluster_id]] <- res$result
    all_qc[[cluster_id]] <- res$qc
    
    if (nrow(res$skipped) > 0) {
      skipped[[cluster_id]] <- res$skipped
    }
  }
  
  list(
    de_all = bind_rows(all_results),
    qc = bind_rows(all_qc),
    skipped = bind_rows(skipped)
  )
}


# ----------------------------- #
# 3. Run pipeline
# ----------------------------- #

setup_dirs(config)

message("Loading module-scored baseline EBV object...")
merged <- load_seurat_object(config$input_rds)

message("Checking metadata...")
check_required_metadata(merged, config)

message("Subsetting cells for lifecycle: ", config$lifecycle_label)
obj_sub <- subset_cells_for_analysis(merged, config)

message("Running pseudobulk edgeR by cluster...")
edgeR_results <- run_edger_by_cluster(
  obj = obj_sub,
  cfg = config
)

de_all <- edgeR_results$de_all
qc_df <- edgeR_results$qc
skipped_df <- edgeR_results$skipped

de_sig <- de_all %>%
  filter(!is.na(p_val_adj), p_val_adj < config$fdr_threshold)

write.csv(
  de_all,
  config$all_de_csv,
  row.names = FALSE
)

write.csv(
  de_sig,
  config$sig_de_csv,
  row.names = FALSE
)

write.csv(
  qc_df,
  config$qc_csv,
  row.names = FALSE
)

if (nrow(skipped_df) == 0) {
  skipped_df <- tibble(
    lifecycle = character(),
    cluster = character(),
    reason = character(),
    n_control = integer(),
    n_ms = integer()
  )
}

write.csv(
  skipped_df,
  config$skipped_clusters_csv,
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  config$session_info_file
)

message("\nPseudobulk edgeR analysis complete.")
message("Lifecycle analysed: ", config$lifecycle_label)
message("Saved all DE results to: ", config$all_de_csv)
message("Saved significant DE results to: ", config$sig_de_csv)
message("Saved QC table to: ", config$qc_csv)
message("Saved skipped cluster log to: ", config$skipped_clusters_csv)
message("Saved session info to: ", config$session_info_file)