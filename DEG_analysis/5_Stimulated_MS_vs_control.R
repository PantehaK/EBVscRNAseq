#!/usr/bin/env Rscript

# ==============================================================================
# edgeR pseudobulk: MS vs Control in stimulated latent EBV-specific T cells
# ==============================================================================

# Purpose:
#   This script performs pseudobulk differential expression comparing MS versus
#   Control participants within stimulated latent EBV-specific T cells, separately
#   for each activated-data cluster.
#
#   Starting object:
#     3_activated_EBV_module_scored_reannotated.rds
#
#   Group definition:
#     Baseline Latent:
#       batch starts with GEMEBV AND lifecycle == Latent
#
#     Stimulated Latent:
#       batch starts with EBVLCL AND lifecycle == Latent
#
#     Non-Latent:
#       all other cells
#
#   Analysis subset:
#     Stimulated Latent cells only, MS and Control participants.
#
#   Pseudobulk unit:
#     donor/sample within each cluster.
#
#   Model:
#     ~ covariates + cohort
#
#   Default model:
#     ~ age + cohort
#
#   Because Control is the reference:
#     logFC > 0 = higher in MS
#     logFC < 0 = higher in Control
#
# Notes:
#   - Paths are intentionally generic for GitHub.
#   - Update the config block before running.
#   - The cluster column defaults to `new_cluster`.
#   - The donor/sample column defaults to `sample`; change to `id` if needed.
#   - The default covariates are age only, matching the supplied analysis.
#     To adjust for sex as well, set config$covariates = c("age", "sex").
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
  library(edgeR)
  library(Matrix)
  library(readr)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  
  # Input RDS
  input_rds = "path/to/Publication_RDS/Activated_EBV_integrated_seurat.rds",
  
  # Output directory
  output_dir = "path/to/Publication_data/EBV_activated/DEG/edgeR_StimulatedLatent_MS_vs_Control_by_cluster",
  
  # Assay and count layer/slot
  assay = "RNA",
  counts_layer = "counts",
  
  # Metadata columns
  donor_col = "sample",
  cohort_col = "cohort",
  cluster_col = "new_cluster",
  batch_col = "batch",
  lifecycle_col = "lifecycle",
  age_col = "age",
  sex_col = "sex",
  
  # Group assignment
  baseline_batch_pattern = "^GEMEBV",
  stimulated_batch_pattern = "^EBVLCL",
  latent_lifecycle_label = "Latent",
  baseline_label = "Baseline Latent",
  stimulated_label = "Stimulated Latent",
  non_latent_label = "Non-Latent",
  
  # Cohort labels and reference
  control_label = "Control",
  ms_label = "MS",
  
  # Covariates in the edgeR design.
  # Default matches the original age-adjusted analysis.
  # Change to c("age", "sex") if sex adjustment is desired and supported.
  covariates = c("age"),
  
  # Minimum donors per cohort required per cluster
  min_donors_per_cohort = 3,
  
  # edgeR settings
  robust_ql_fit = TRUE,
  
  # Significance threshold
  fdr_threshold = 0.05,
  
  # Optional cluster order. Leave as NULL to use sorted values from the data.
  cluster_order = NULL,
  
  # Output files
  filtered_metadata_csv = "Age_adjusted_edgeR_StimulatedLatent_filtered_metadata.csv",
  qc_donors_csv = "Age_adjusted_edgeR_QC_StimulatedLatent_donors_per_group_by_cluster.csv",
  qc_cells_csv = "Age_adjusted_edgeR_QC_StimulatedLatent_cells_per_sample_by_cluster.csv",
  skipped_clusters_csv = "Age_adjusted_edgeR_StimulatedLatent_skipped_clusters.csv",
  all_results_csv = "Age_adjusted_edgeR_StimulatedLatent_MS_vs_Control_by_cluster_ALL.csv",
  significant_results_csv = "Age_adjusted_edgeR_StimulatedLatent_MS_vs_Control_by_cluster_FDR0.05.csv",
  significant_summary_csv = "Age_adjusted_edgeR_StimulatedLatent_significant_gene_summary_by_cluster.csv",
  session_info_file = "sessionInfo_edgeR_StimulatedLatent_MS_vs_Control_by_cluster.txt"
)


# ----------------------------- #
# 2. Helper functions
# ----------------------------- #

create_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}


safe_name <- function(x) {
  
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+$", "") %>%
    str_replace_all("^_+", "")
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


get_assay_counts_compat <- function(object, assay = "RNA", counts_layer = "counts") {
  
  # Compatible with Seurat v5 and Seurat v4.
  tryCatch(
    {
      GetAssayData(object, assay = assay, layer = counts_layer)
    },
    error = function(e) {
      GetAssayData(object, assay = assay, slot = counts_layer)
    }
  )
}


clean_chr <- function(x) {
  
  x <- str_trim(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL", "null", "None", "none")] <- NA_character_
  x
}


clean_cohort <- function(x, control_label = "Control", ms_label = "MS") {
  
  x_raw <- as.character(x)
  x_clean <- str_to_lower(str_trim(x_raw))
  
  case_when(
    x_clean %in% c(
      "control", "controls", "healthy control", "healthy controls",
      "ctrl", "hc", "nms", "non-ms", "nonms", "non_ms"
    ) ~ control_label,
    
    x_clean %in% c(
      "ms", "rrms", "pwms", "multiple sclerosis"
    ) ~ ms_label,
    
    TRUE ~ x_raw
  )
}


assign_direction <- function(logfc, fdr, fdr_threshold) {
  
  case_when(
    !is.na(fdr) & fdr < fdr_threshold & !is.na(logfc) & logfc > 0 ~ "MS_up",
    !is.na(fdr) & fdr < fdr_threshold & !is.na(logfc) & logfc < 0 ~ "Control_up",
    TRUE ~ "NS"
  )
}


coerce_covariates <- function(df, covariates) {
  
  for (cv in covariates) {
    if (!cv %in% colnames(df)) {
      stop("Covariate column is missing from sample metadata: ", cv)
    }
    
    if (cv %in% c("age", "Age")) {
      df[[cv]] <- suppressWarnings(as.numeric(df[[cv]]))
    } else {
      # Treat non-numeric covariates as factors unless already numeric.
      if (!is.numeric(df[[cv]])) {
        df[[cv]] <- factor(df[[cv]])
      }
    }
  }
  
  df
}


make_design_formula <- function(covariates) {
  
  rhs <- c(covariates, "cohort")
  
  if (length(rhs) == 0) {
    rhs <- "1"
  }
  
  as.formula(paste("~", paste(rhs, collapse = " + ")))
}


make_empty_result_table <- function() {
  
  tibble(
    latent_group = character(),
    cluster = character(),
    comparison = character(),
    gene = character(),
    n_donors_MS = numeric(),
    n_donors_Control = numeric(),
    n_pseudobulk_samples = numeric(),
    n_cells = numeric(),
    avg_log2FC = numeric(),
    logCPM = numeric(),
    F = numeric(),
    p_val = numeric(),
    p_val_adj = numeric(),
    direction = character()
  )
}


# ----------------------------- #
# 3. Load object and prepare metadata
# ----------------------------- #

create_dir(config$output_dir)

message("Loading activated EBV object:")
message(config$input_rds)

merged <- load_seurat_object(config$input_rds)

required_cols <- c(
  config$donor_col,
  config$cohort_col,
  config$cluster_col,
  config$batch_col,
  config$lifecycle_col,
  config$age_col
)

# Only require sex if requested as a covariate or explicitly present in config.
if ("sex" %in% config$covariates) {
  required_cols <- unique(c(required_cols, config$sex_col))
}

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

DefaultAssay(merged) <- config$assay

counts_all <- get_assay_counts_compat(
  object = merged,
  assay = config$assay,
  counts_layer = config$counts_layer
)

meta_all <- merged@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell_barcode") %>%
  mutate(
    donor_id = clean_chr(.data[[config$donor_col]]),
    cohort = clean_cohort(
      .data[[config$cohort_col]],
      control_label = config$control_label,
      ms_label = config$ms_label
    ),
    cluster = clean_chr(.data[[config$cluster_col]]),
    batch_for_de = clean_chr(.data[[config$batch_col]]),
    lifecycle_for_de = clean_chr(.data[[config$lifecycle_col]]),
    age = suppressWarnings(as.numeric(.data[[config$age_col]])),
    
    sex = if (config$sex_col %in% colnames(.)) {
      clean_chr(.data[[config$sex_col]])
    } else {
      NA_character_
    },
    
    latent_group = case_when(
      !is.na(batch_for_de) &
        grepl(config$baseline_batch_pattern, batch_for_de) &
        lifecycle_for_de == config$latent_lifecycle_label ~ config$baseline_label,
      
      !is.na(batch_for_de) &
        grepl(config$stimulated_batch_pattern, batch_for_de) &
        lifecycle_for_de == config$latent_lifecycle_label ~ config$stimulated_label,
      
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

if (!identical(colnames(counts_all), meta_all$cell_barcode)) {
  stop("Counts columns and metadata cell barcodes are not aligned.")
}

message("Latent group counts:")
print(table(meta_all$latent_group, useNA = "ifany"))


# ----------------------------- #
# 4. Filter to stimulated latent MS/Control cells
# ----------------------------- #

meta_use <- meta_all %>%
  filter(
    latent_group == config$stimulated_label,
    cohort %in% c(config$control_label, config$ms_label),
    !is.na(donor_id), donor_id != "",
    !is.na(cluster), cluster != ""
  )

if (nrow(meta_use) == 0) {
  stop("No cells remained after stimulated latent/cohort/donor/cluster filtering.")
}

counts_use <- counts_all[, meta_use$cell_barcode, drop = FALSE]

filtered_metadata_path <- file.path(
  config$output_dir,
  config$filtered_metadata_csv
)

write_csv(meta_use, filtered_metadata_path)

message("Cells retained for stimulated latent MS vs Control:")
print(table(meta_use$cohort, useNA = "ifany"))

message("Cells by cluster and cohort:")
print(table(meta_use$cluster, meta_use$cohort, useNA = "ifany"))


# ----------------------------- #
# 5. Determine clusters to test
# ----------------------------- #

if (is.null(config$cluster_order)) {
  clusters_to_test <- sort(unique(meta_use$cluster))
} else {
  clusters_to_test <- config$cluster_order[
    config$cluster_order %in% unique(meta_use$cluster)
  ]
  
  unexpected_clusters <- setdiff(
    sort(unique(meta_use$cluster)),
    config$cluster_order
  )
  
  clusters_to_test <- c(clusters_to_test, unexpected_clusters)
}


# ----------------------------- #
# 6. edgeR function per cluster
# ----------------------------- #

run_edger_one_cluster <- function(cluster_label) {
  
  message("Running cluster: ", cluster_label)
  
  meta_ct_cells <- meta_use %>%
    filter(cluster == cluster_label)
  
  qc_cells_ct <- meta_ct_cells %>%
    count(cluster, cohort, donor_id, name = "n_cells") %>%
    arrange(cluster, cohort, desc(n_cells))
  
  donors_per_group <- meta_ct_cells %>%
    distinct(donor_id, cohort) %>%
    count(cohort, name = "n_donors") %>%
    complete(
      cohort = c(config$control_label, config$ms_label),
      fill = list(n_donors = 0)
    ) %>%
    arrange(cohort)
  
  donors_wide <- donors_per_group %>%
    deframe()
  
  n_control <- unname(donors_wide[config$control_label])
  n_ms <- unname(donors_wide[config$ms_label])
  
  qc_donors_ct <- donors_per_group %>%
    mutate(cluster = cluster_label) %>%
    relocate(cluster)
  
  if (
    is.na(n_control) || is.na(n_ms) ||
    n_control < config$min_donors_per_cohort ||
    n_ms < config$min_donors_per_cohort
  ) {
    message(
      "  Skipping ",
      cluster_label,
      " (need >= ",
      config$min_donors_per_cohort,
      " donors/group; have Control=",
      n_control,
      ", MS=",
      n_ms,
      ")"
    )
    
    return(list(
      result = NULL,
      qc_donors = qc_donors_ct,
      qc_cells = qc_cells_ct,
      skipped = tibble(
        cluster = cluster_label,
        reason = "insufficient_donors_per_cohort",
        n_Control = n_control,
        n_MS = n_ms,
        n_cells = nrow(meta_ct_cells)
      )
    ))
  }
  
  # Pseudobulk counts per donor within this cluster.
  donor_factor <- factor(meta_ct_cells$donor_id)
  
  aggregation_matrix <- Matrix::sparseMatrix(
    i = seq_along(donor_factor),
    j = as.integer(donor_factor),
    x = 1,
    dims = c(length(donor_factor), nlevels(donor_factor)),
    dimnames = list(meta_ct_cells$cell_barcode, levels(donor_factor))
  )
  
  pb_counts <- counts_use[, meta_ct_cells$cell_barcode, drop = FALSE] %*% aggregation_matrix
  pb_counts <- as.matrix(pb_counts)
  
  smeta <- meta_ct_cells %>%
    select(donor_id, cohort, age, sex) %>%
    distinct(donor_id, .keep_all = TRUE) %>%
    filter(donor_id %in% colnames(pb_counts)) %>%
    arrange(match(donor_id, colnames(pb_counts))) %>%
    mutate(
      cohort = factor(
        cohort,
        levels = c(config$control_label, config$ms_label)
      )
    )
  
  smeta <- coerce_covariates(smeta, config$covariates)
  
  # Drop samples missing covariates.
  smeta <- smeta %>%
    filter(if_all(all_of(config$covariates), ~ !is.na(.x)))
  
  pb_counts <- pb_counts[, smeta$donor_id, drop = FALSE]
  
  if (!identical(colnames(pb_counts), smeta$donor_id)) {
    stop("Pseudobulk columns and sample metadata are not aligned for cluster: ", cluster_label)
  }
  
  # Re-check donor counts after covariate filtering.
  donors_after_covariate_filter <- smeta %>%
    count(cohort, name = "n_donors") %>%
    complete(
      cohort = factor(
        c(config$control_label, config$ms_label),
        levels = c(config$control_label, config$ms_label)
      ),
      fill = list(n_donors = 0)
    )
  
  donors_after_wide <- donors_after_covariate_filter %>%
    mutate(cohort = as.character(cohort)) %>%
    deframe()
  
  n_control_after <- unname(donors_after_wide[config$control_label])
  n_ms_after <- unname(donors_after_wide[config$ms_label])
  
  if (
    is.na(n_control_after) || is.na(n_ms_after) ||
    n_control_after < config$min_donors_per_cohort ||
    n_ms_after < config$min_donors_per_cohort
  ) {
    message(
      "  Skipping ",
      cluster_label,
      " after covariate filtering (Control=",
      n_control_after,
      ", MS=",
      n_ms_after,
      ")"
    )
    
    return(list(
      result = NULL,
      qc_donors = qc_donors_ct,
      qc_cells = qc_cells_ct,
      skipped = tibble(
        cluster = cluster_label,
        reason = "insufficient_donors_after_covariate_filtering",
        n_Control = n_control_after,
        n_MS = n_ms_after,
        n_cells = nrow(meta_ct_cells)
      )
    ))
  }
  
  y <- DGEList(
    counts = pb_counts,
    samples = smeta
  )
  
  y <- calcNormFactors(y)
  
  keep <- filterByExpr(
    y,
    group = y$samples$cohort
  )
  
  y <- y[keep, , keep.lib.sizes = FALSE]
  
  if (nrow(y) == 0) {
    message("  Skipping ", cluster_label, ": no genes retained after filterByExpr.")
    
    return(list(
      result = NULL,
      qc_donors = qc_donors_ct,
      qc_cells = qc_cells_ct,
      skipped = tibble(
        cluster = cluster_label,
        reason = "no_genes_retained_after_filterByExpr",
        n_Control = n_control_after,
        n_MS = n_ms_after,
        n_cells = nrow(meta_ct_cells)
      )
    ))
  }
  
  design_formula <- make_design_formula(config$covariates)
  design <- model.matrix(design_formula, data = y$samples)
  
  if (qr(design)$rank < ncol(design)) {
    message("  Skipping ", cluster_label, ": design matrix is not full rank.")
    
    return(list(
      result = NULL,
      qc_donors = qc_donors_ct,
      qc_cells = qc_cells_ct,
      skipped = tibble(
        cluster = cluster_label,
        reason = "design_matrix_not_full_rank",
        n_Control = n_control_after,
        n_MS = n_ms_after,
        n_cells = nrow(meta_ct_cells)
      )
    ))
  }
  
  y <- estimateDisp(y, design)
  
  fit <- glmQLFit(
    y,
    design,
    robust = config$robust_ql_fit
  )
  
  coef_name <- paste0("cohort", config$ms_label)
  
  if (!(coef_name %in% colnames(design))) {
    message("  Skipping ", cluster_label, ": coefficient ", coef_name, " not found.")
    
    return(list(
      result = NULL,
      qc_donors = qc_donors_ct,
      qc_cells = qc_cells_ct,
      skipped = tibble(
        cluster = cluster_label,
        reason = paste0("coefficient_not_found_", coef_name),
        n_Control = n_control_after,
        n_MS = n_ms_after,
        n_cells = nrow(meta_ct_cells)
      )
    ))
  }
  
  qlf <- glmQLFTest(
    fit,
    coef = which(colnames(design) == coef_name)
  )
  
  tab <- topTags(qlf, n = Inf)$table %>%
    as.data.frame() %>%
    rownames_to_column("gene")
  
  if ("logFC" %in% colnames(tab)) {
    tab <- rename(tab, avg_log2FC = logFC)
  } else {
    tab <- mutate(tab, avg_log2FC = NA_real_)
  }
  
  if ("PValue" %in% colnames(tab)) {
    tab <- rename(tab, p_val = PValue)
  }
  
  if ("FDR" %in% colnames(tab)) {
    tab <- rename(tab, p_val_adj = FDR)
  }
  
  results_ct <- tab %>%
    mutate(
      latent_group = config$stimulated_label,
      cluster = cluster_label,
      comparison = paste0(config$ms_label, "_vs_", config$control_label),
      covariates = paste(config$covariates, collapse = " + "),
      n_donors_MS = n_ms_after,
      n_donors_Control = n_control_after,
      n_pseudobulk_samples = ncol(pb_counts),
      n_cells = nrow(meta_ct_cells),
      direction = assign_direction(
        logfc = avg_log2FC,
        fdr = p_val_adj,
        fdr_threshold = config$fdr_threshold
      )
    ) %>%
    relocate(
      latent_group,
      cluster,
      comparison,
      covariates,
      gene,
      n_donors_MS,
      n_donors_Control,
      n_pseudobulk_samples,
      n_cells,
      avg_log2FC
    ) %>%
    arrange(p_val_adj, p_val)
  
  per_cluster_dir <- file.path(config$output_dir, "per_cluster")
  create_dir(per_cluster_dir)
  
  output_file <- file.path(
    per_cluster_dir,
    paste0(
      "edgeR_",
      safe_name(config$stimulated_label),
      "_",
      config$ms_label,
      "_vs_",
      config$control_label,
      "_cluster_",
      safe_name(cluster_label),
      ".csv"
    )
  )
  
  write_csv(results_ct, output_file)
  
  list(
    result = results_ct,
    qc_donors = qc_donors_ct,
    qc_cells = qc_cells_ct,
    skipped = tibble()
  )
}


# ----------------------------- #
# 7. Run edgeR across clusters
# ----------------------------- #

run_list <- lapply(
  clusters_to_test,
  run_edger_one_cluster
)

de_all <- bind_rows(
  lapply(run_list, `[[`, "result")
)

qc_donors <- bind_rows(
  lapply(run_list, `[[`, "qc_donors")
)

qc_cells <- bind_rows(
  lapply(run_list, `[[`, "qc_cells")
)

skipped_clusters <- bind_rows(
  lapply(run_list, `[[`, "skipped")
)

if (nrow(de_all) == 0) {
  warning("No edgeR results were generated for any cluster.")
  de_all <- make_empty_result_table()
}

if (nrow(skipped_clusters) == 0) {
  skipped_clusters <- tibble(
    cluster = character(),
    reason = character(),
    n_Control = numeric(),
    n_MS = numeric(),
    n_cells = numeric()
  )
}


# ----------------------------- #
# 8. Significant genes
# ----------------------------- #

de_sig <- de_all %>%
  filter(
    !is.na(p_val_adj),
    p_val_adj < config$fdr_threshold
  )

sig_summary <- de_sig %>%
  count(
    latent_group,
    cluster,
    direction,
    name = "n_genes"
  ) %>%
  arrange(latent_group, cluster, direction)


# ----------------------------- #
# 9. Save outputs
# ----------------------------- #

all_results_path <- file.path(
  config$output_dir,
  config$all_results_csv
)

sig_results_path <- file.path(
  config$output_dir,
  config$significant_results_csv
)

qc_donors_path <- file.path(
  config$output_dir,
  config$qc_donors_csv
)

qc_cells_path <- file.path(
  config$output_dir,
  config$qc_cells_csv
)

skipped_clusters_path <- file.path(
  config$output_dir,
  config$skipped_clusters_csv
)

sig_summary_path <- file.path(
  config$output_dir,
  config$significant_summary_csv
)

write_csv(de_all, all_results_path)
write_csv(de_sig, sig_results_path)
write_csv(qc_donors, qc_donors_path)
write_csv(qc_cells, qc_cells_path)
write_csv(skipped_clusters, skipped_clusters_path)
write_csv(sig_summary, sig_summary_path)


# ----------------------------- #
# 10. Session information
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
    capture.output(print(table(meta_all$latent_group, useNA = "ifany"))),
    "",
    "Cells retained for stimulated latent MS vs Control:",
    capture.output(print(table(meta_use$cohort, useNA = "ifany"))),
    "",
    "Cluster x cohort table:",
    capture.output(print(table(meta_use$cluster, meta_use$cohort, useNA = "ifany"))),
    "",
    "QC donors:",
    capture.output(print(qc_donors)),
    "",
    "Skipped clusters:",
    capture.output(print(skipped_clusters)),
    "",
    "Significant gene summary:",
    capture.output(print(sig_summary)),
    "",
    "Output files:",
    filtered_metadata_path,
    all_results_path,
    sig_results_path,
    qc_donors_path,
    qc_cells_path,
    skipped_clusters_path,
    sig_summary_path,
    "",
    "Session information:",
    capture.output(sessionInfo())
  ),
  session_info_path
)


# ----------------------------- #
# 11. Completion messages
# ----------------------------- #

message("\nStimulated latent edgeR MS vs Control by cluster workflow complete.")
message("Filtered metadata: ", filtered_metadata_path)
message("All DE results: ", all_results_path)
message("Significant DE results: ", sig_results_path)
message("QC donors: ", qc_donors_path)
message("QC cells: ", qc_cells_path)
message("Skipped clusters: ", skipped_clusters_path)
message("Significant summary: ", sig_summary_path)
message("Session info: ", session_info_path)

message("\nSignificant gene summary:")
print(sig_summary)
