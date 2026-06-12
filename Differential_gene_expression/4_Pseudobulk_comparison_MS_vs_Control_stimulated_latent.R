#!/usr/bin/env Rscript

# ==============================================================================
# Pseudobulk edgeR DE in stimulated latent EBV-specific T cells by cluster
# ==============================================================================
#
# Purpose:
#   This script performs pseudobulk differential expression using edgeR to compare
#   MS versus Control within stimulated latent EBV-specific T cells, separately
#   for each matched UMAP/Harmony cluster.
#
#   It performs:
#     1. latent_group annotation,
#     2. subsetting to stimulated latent cells,
#     3. sample-level pseudobulk aggregation within each matched cluster,
#     4. edgeR quasi-likelihood negative-binomial modelling,
#     5. age-adjusted MS versus Control testing,
#     6. donor/cell QC export,
#     7. skipped-cluster logging.
#
# Expected input:
#   Activated paired EBV CD8+ T cell object, for example:
#     15_activated_EBV_MS_stimLat_signature_scored.rds
#
# Notes:
#   - logFC is MS versus Control.
#   - Clusters with fewer than min_donors_per_cohort per group are skipped.
#   - Default model is:
#       ~ age + cohort
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
  dataset_id = "activated_ebv_stimulated_latent",
  
  input_rds = "path/to/input/15_activated_EBV_MS_stimLat_signature_scored.rds",
  
  output_dir = "path/to/output/Publication_data/EBV_activated/pseudobulk_edgeR_stimulated_latent_by_cluster",
  
  all_de_csv = "path/to/output/Publication_data/EBV_activated/pseudobulk_edgeR_stimulated_latent_by_cluster/Age_adjusted_EdgeR_StimulatedLatent_MS_vs_Control_by_matched_cluster_ALL.csv",
  sig_de_csv = "path/to/output/Publication_data/EBV_activated/pseudobulk_edgeR_stimulated_latent_by_cluster/Age_adjusted_EdgeR_StimulatedLatent_MS_vs_Control_by_matched_cluster_FDR0.05.csv",
  qc_donors_csv = "path/to/output/Publication_data/EBV_activated/pseudobulk_edgeR_stimulated_latent_by_cluster/Age_adjusted_EdgeR_QC_StimulatedLatent_donors_per_group_by_matched_cluster.csv",
  qc_cells_csv = "path/to/output/Publication_data/EBV_activated/pseudobulk_edgeR_stimulated_latent_by_cluster/Age_adjusted_EdgeR_QC_StimulatedLatent_cells_per_sample_by_matched_cluster.csv",
  skipped_clusters_csv = "path/to/output/Publication_data/EBV_activated/pseudobulk_edgeR_stimulated_latent_by_cluster/Age_adjusted_EdgeR_StimulatedLatent_skipped_clusters.csv",
  session_info_file = "path/to/output/Publication_data/EBV_activated/pseudobulk_edgeR_stimulated_latent_by_cluster/sessionInfo_edgeR_stimLat_by_cluster.txt",
  
  # Metadata columns.
  sample_col = "sample",
  cohort_col = "cohort",
  cluster_col = "matched_cluster",
  batch_col = "batch",
  lifecycle_col = "lifecycle",
  
  # Group labels.
  latent_label = "Latent",
  non_latent_label = "Non-Latent",
  baseline_latent_label = "Baseline Latent",
  stimulated_latent_label = "Stimulated Latent",
  
  baseline_batch_pattern = "^GEMEBV",
  stimulated_batch_pattern = "^EBVLCL",
  
  # Cohort comparison.
  cohort_levels = c("Control", "MS"),
  contrast_coef = "cohortMS",
  
  # Assay settings.
  assay = "RNA",
  aggregate_slot = "counts",
  
  # edgeR settings.
  min_donors_per_cohort = 3,
  robust = TRUE,
  
  # Covariates.
  # Original pasted code effectively used age only:
  #   model.matrix(~ age + cohort)
  #
  # To include sex as well, change this to:
  #   covariates = c("age", "sex")
  covariates = c("age"),
  
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
  create_parent_dir(cfg$qc_donors_csv)
  create_parent_dir(cfg$qc_cells_csv)
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
    cfg$batch_col,
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


subset_stimulated_latent_cells <- function(obj, cfg) {
  md <- obj@meta.data
  
  keep_cells <- rownames(md)[
    md$latent_group == cfg$stimulated_latent_label &
      md[[cfg$cohort_col]] %in% cfg$cohort_levels
  ]
  
  if (length(keep_cells) == 0) {
    stop("No stimulated latent cells found for MS/Control comparison.")
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


count_cells_per_sample <- function(obj_cluster, cluster_id, cfg) {
  obj_cluster@meta.data %>%
    count(
      matched_cluster = .data[[cfg$cluster_col]],
      cohort = .data[[cfg$cohort_col]],
      sample = .data[[cfg$sample_col]],
      name = "n_cells"
    ) %>%
    mutate(
      matched_cluster = cluster_id
    ) %>%
    arrange(matched_cluster, cohort, desc(n_cells))
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
      latent_group = cfg$stimulated_latent_label,
      matched_cluster = cluster_id,
      .before = 1
    ) %>%
    arrange(cohort)
}


aggregate_counts_by_sample <- function(obj_cluster, cfg) {
  DefaultAssay(obj_cluster) <- cfg$assay
  
  AggregateExpression(
    object = obj_cluster,
    assays = cfg$assay,
    group.by = cfg$sample_col,
    slot = cfg$aggregate_slot,
    return.seurat = FALSE
  )[[cfg$assay]]
}


make_sample_metadata <- function(obj_cluster, pb, cfg) {
  smeta <- obj_cluster@meta.data %>%
    select(
      sample = .data[[cfg$sample_col]],
      cohort = .data[[cfg$cohort_col]],
      any_of(cfg$covariates)
    ) %>%
    distinct(sample, .keep_all = TRUE) %>%
    filter(sample %in% colnames(pb)) %>%
    arrange(match(sample, colnames(pb)))
  
  if (!all(smeta$sample == colnames(pb))) {
    stop("Sample metadata does not match pseudobulk matrix columns.")
  }
  
  smeta <- smeta %>%
    mutate(
      cohort = factor(cohort, levels = cfg$cohort_levels)
    )
  
  for (covar in cfg$covariates) {
    if (covar == "sex") {
      smeta[[covar]] <- factor(smeta[[covar]])
    } else {
      smeta[[covar]] <- as.numeric(smeta[[covar]])
    }
  }
  
  smeta <- smeta %>%
    filter(
      !if_any(
        all_of(c("cohort", cfg$covariates)),
        is.na
      )
    )
  
  smeta
}


count_donors_after_covariate_filter <- function(smeta, cluster_id, cfg) {
  smeta %>%
    distinct(sample, cohort) %>%
    count(cohort, name = "n_donors_after_covariate_filter") %>%
    complete(
      cohort = cfg$cohort_levels,
      fill = list(n_donors_after_covariate_filter = 0)
    ) %>%
    mutate(
      matched_cluster = cluster_id,
      .before = 1
    )
}


make_design_formula <- function(cfg) {
  rhs <- c(cfg$covariates, "cohort")
  
  as.formula(
    paste("~", paste(rhs, collapse = " + "))
  )
}


run_edger_one_cluster <- function(obj_cluster, cluster_id, cfg) {
  qc_cells <- count_cells_per_sample(
    obj_cluster = obj_cluster,
    cluster_id = cluster_id,
    cfg = cfg
  )
  
  qc_donors_pre <- count_donors_per_group(
    obj_cluster = obj_cluster,
    cluster_id = cluster_id,
    cfg = cfg
  )
  
  donors_wide_pre <- qc_donors_pre %>%
    select(cohort, n_donors) %>%
    deframe()
  
  n_control_pre <- donors_wide_pre[[cfg$cohort_levels[1]]]
  n_ms_pre <- donors_wide_pre[[cfg$cohort_levels[2]]]
  
  if (
    is.null(n_control_pre) || is.null(n_ms_pre) ||
    n_control_pre < cfg$min_donors_per_cohort ||
    n_ms_pre < cfg$min_donors_per_cohort
  ) {
    return(list(
      result = tibble(),
      qc_donors = qc_donors_pre,
      qc_cells = qc_cells,
      skipped = tibble(
        latent_group = cfg$stimulated_latent_label,
        matched_cluster = cluster_id,
        reason = "insufficient_donors_before_covariate_filter",
        n_control = n_control_pre,
        n_ms = n_ms_pre
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
  
  pb <- pb[, smeta$sample, drop = FALSE]
  
  qc_donors_post <- count_donors_after_covariate_filter(
    smeta = smeta,
    cluster_id = cluster_id,
    cfg = cfg
  )
  
  donors_wide_post <- qc_donors_post %>%
    select(cohort, n_donors_after_covariate_filter) %>%
    deframe()
  
  n_control_post <- donors_wide_post[[cfg$cohort_levels[1]]]
  n_ms_post <- donors_wide_post[[cfg$cohort_levels[2]]]
  
  qc_donors <- qc_donors_pre %>%
    left_join(qc_donors_post, by = c("matched_cluster", "cohort"))
  
  if (
    is.null(n_control_post) || is.null(n_ms_post) ||
    n_control_post < cfg$min_donors_per_cohort ||
    n_ms_post < cfg$min_donors_per_cohort
  ) {
    return(list(
      result = tibble(),
      qc_donors = qc_donors,
      qc_cells = qc_cells,
      skipped = tibble(
        latent_group = cfg$stimulated_latent_label,
        matched_cluster = cluster_id,
        reason = "insufficient_donors_after_covariate_filter",
        n_control = n_control_post,
        n_ms = n_ms_post
      )
    ))
  }
  
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
      qc_donors = qc_donors,
      qc_cells = qc_cells,
      skipped = tibble(
        latent_group = cfg$stimulated_latent_label,
        matched_cluster = cluster_id,
        reason = "no_genes_passed_filterByExpr",
        n_control = n_control_post,
        n_ms = n_ms_post
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
      qc_donors = qc_donors,
      qc_cells = qc_cells,
      skipped = tibble(
        latent_group = cfg$stimulated_latent_label,
        matched_cluster = cluster_id,
        reason = paste0("contrast_coef_not_found: ", cfg$contrast_coef),
        n_control = n_control_post,
        n_ms = n_ms_post
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
      qc_donors = qc_donors,
      qc_cells = qc_cells,
      skipped = tibble(
        latent_group = cfg$stimulated_latent_label,
        matched_cluster = cluster_id,
        reason = paste0("edgeR_failed: ", fit_result$message),
        n_control = n_control_post,
        n_ms = n_ms_post
      )
    ))
  }
  
  tab <- as.data.frame(fit_result) %>%
    rownames_to_column("gene") %>%
    rename(
      avg_log2FC = logFC,
      p_val = PValue,
      p_val_adj = FDR
    ) %>%
    mutate(
      latent_group = cfg$stimulated_latent_label,
      matched_cluster = cluster_id,
      n_donors_MS = n_ms_post,
      n_donors_Control = n_control_post
    ) %>%
    relocate(
      latent_group,
      matched_cluster,
      gene,
      n_donors_MS,
      n_donors_Control,
      avg_log2FC
    )
  
  list(
    result = tab,
    qc_donors = qc_donors,
    qc_cells = qc_cells,
    skipped = tibble()
  )
}


run_edger_by_cluster <- function(obj, cfg) {
  clusters <- sort(unique(obj[[cfg$cluster_col, drop = TRUE]]))
  
  all_results <- list()
  all_qc_donors <- list()
  all_qc_cells <- list()
  skipped <- list()
  
  for (cluster_id in clusters) {
    message("Running matched_cluster: ", cluster_id)
    
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
    all_qc_donors[[cluster_id]] <- res$qc_donors
    all_qc_cells[[cluster_id]] <- res$qc_cells
    
    if (nrow(res$skipped) > 0) {
      skipped[[cluster_id]] <- res$skipped
    }
  }
  
  list(
    de_all = bind_rows(all_results),
    qc_donors = bind_rows(all_qc_donors),
    qc_cells = bind_rows(all_qc_cells),
    skipped = bind_rows(skipped)
  )
}


# ----------------------------- #
# 3. Run pipeline
# ----------------------------- #

setup_dirs(config)

message("Loading activated EBV object...")
merged <- load_seurat_object(config$input_rds)

message("Checking metadata...")
check_required_metadata(merged, config)

message("Defining latent_group...")
merged <- add_latent_group(merged, config)

message("Subsetting to stimulated latent cells...")
stim_lat <- subset_stimulated_latent_cells(merged, config)

message("Running pseudobulk edgeR by matched cluster...")
edgeR_results <- run_edger_by_cluster(
  obj = stim_lat,
  cfg = config
)

de_all <- edgeR_results$de_all
de_sig <- de_all %>%
  filter(!is.na(p_val_adj), p_val_adj < config$fdr_threshold)

qc_donors <- edgeR_results$qc_donors
qc_cells <- edgeR_results$qc_cells
skipped_clusters <- edgeR_results$skipped

if (nrow(skipped_clusters) == 0) {
  skipped_clusters <- tibble(
    latent_group = character(),
    matched_cluster = character(),
    reason = character(),
    n_control = integer(),
    n_ms = integer()
  )
}

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
  qc_donors,
  config$qc_donors_csv,
  row.names = FALSE
)

write.csv(
  qc_cells,
  config$qc_cells_csv,
  row.names = FALSE
)

write.csv(
  skipped_clusters,
  config$skipped_clusters_csv,
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  config$session_info_file
)

message("\nStimulated latent pseudobulk edgeR analysis complete.")
message("Saved all DE results to: ", config$all_de_csv)
message("Saved significant DE results to: ", config$sig_de_csv)
message("Saved donor QC to: ", config$qc_donors_csv)
message("Saved cell QC to: ", config$qc_cells_csv)
message("Saved skipped-cluster log to: ", config$skipped_clusters_csv)
message("Saved session info to: ", config$session_info_file)