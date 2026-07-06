#!/usr/bin/env Rscript

# ==============================================================================
# edgeR pseudobulk: MS vs Control across Latent, Lytic and negative cells
# ==============================================================================
#
# Purpose:
#   This script performs pseudobulk differential expression comparing MS versus
#   Control participants within each selected lifecycle/antigen group and within
#   each T cell cluster/cell type.
#
#   Starting object:
#     5_baseline_EBV_cluster_reannotated.rds
#
#   Groups analysed by default:
#     - Latent
#     - Lytic
#     - negative
#
#   Pseudobulk unit:
#     donor/sample within each lifecycle group x celltype_new group
#
#   Model:
#     ~ cohort
#
#   Because Control is the reference:
#     logFC > 0 = higher in MS
#     logFC < 0 = higher in Control
#
# Notes:
#   - Paths are intentionally generic for GitHub.
#   - Update the config block before running.
#   - This version loops over all groups in config$target_groups.
#   - The donor/sample column defaults to `sample`. Change config$donor_col to
#     "id" if your analysis should use the `id` column instead.
#   - This version loops over `celltype_new` by default. The original code used
#     harmony cluster IDs but then filtered by celltype_new, which could mismatch.
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
  library(purrr)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  
  # Input RDS
  input_rds = "path/to/Publication_RDS/Baseline_EBV_integrated_seurat.rds",
  
  # Output directory
  output_dir = "path/to/Publication_data/EBV_baseline/DEG/edgeR_MS_vs_Control_by_lifecycle_celltype",
  
  # Assay and count layer/slot
  assay = "RNA",
  counts_layer = "counts",
  
  # Metadata columns
  donor_col = "sample",
  cohort_col = "cohort",
  celltype_col = "celltype_new",
  lifecycle_col = "lifecycle",
  
  # Lifecycle/antigen groups to run.
  # pattern is matched against lifecycle_col.
  # Set use_regex = FALSE for exact matching.
  target_groups = tibble::tibble(
    target_lifecycle_label = c("Latent", "Lytic", "negative"),
    target_lifecycle_pattern = c("Latent", "Lytic", "negative"),
    use_regex = c(FALSE, FALSE, FALSE),
    ignore_case = c(FALSE, FALSE, FALSE)
  ),
  
  # Cohort labels and reference
  control_label = "Control",
  ms_label = "MS",
  
  # Minimum donors per group required to run edgeR for a cell type
  min_donors_per_cohort = 3,
  
  # edgeR filtering/model settings
  robust_ql_fit = TRUE,
  
  # Significance settings
  fdr_threshold = 0.05,
  abs_logfc_threshold_for_direction = 0,
  
  # Optional fixed cell type order.
  # Leave as NULL to use cell types found in the object.
  celltype_order = c(
    "Naive/early TCM",
    "TEM",
    "Late TEM",
    "CTL",
    "CD69+ early activated T",
    "innate-like T",
    "CD69+ TEM"
  ),
  
  # Output files
  filtered_metadata_csv = "edgeR_MS_vs_Control_filtered_metadata_ALL_lifecycles.csv",
  qc_by_celltype_csv = "edgeR_MS_vs_Control_QC_by_lifecycle_celltype.csv",
  all_results_csv = "edgeR_MS_vs_Control_ALL_by_lifecycle_celltype.csv",
  significant_results_csv = "edgeR_MS_vs_Control_SIGNIFICANT_by_lifecycle_celltype.csv",
  significant_summary_csv = "edgeR_MS_vs_Control_significant_gene_summary_by_lifecycle_celltype.csv",
  skipped_celltypes_csv = "edgeR_MS_vs_Control_skipped_lifecycle_celltypes.csv",
  session_info_file = "sessionInfo_edgeR_MS_vs_Control_by_lifecycle_celltype.txt"
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
  str_trim(as.character(x))
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


match_lifecycle <- function(x, pattern, use_regex = TRUE, ignore_case = FALSE) {
  
  x <- as.character(x)
  
  if (use_regex) {
    grepl(pattern, x, ignore.case = ignore_case)
  } else {
    if (ignore_case) {
      tolower(x) == tolower(pattern)
    } else {
      x == pattern
    }
  }
}


assign_direction <- function(logfc, fdr, fdr_threshold, abs_logfc_threshold) {
  
  case_when(
    !is.na(fdr) &
      fdr < fdr_threshold &
      !is.na(logfc) &
      logfc > abs_logfc_threshold ~ "MS_up",
    
    !is.na(fdr) &
      fdr < fdr_threshold &
      !is.na(logfc) &
      logfc < -abs_logfc_threshold ~ "Control_up",
    
    TRUE ~ "NS"
  )
}


make_empty_result_table <- function() {
  
  tibble(
    target_lifecycle = character(),
    celltype_new = character(),
    comparison = character(),
    gene = character(),
    n_donors_MS = numeric(),
    n_donors_Control = numeric(),
    n_pseudobulk_samples = numeric(),
    n_cells = numeric(),
    avg_log2FC = numeric(),
    logFC = numeric(),
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

message("Loading baseline EBV object:")
message(config$input_rds)

merged <- load_seurat_object(config$input_rds)

required_cols <- c(
  config$donor_col,
  config$cohort_col,
  config$celltype_col,
  config$lifecycle_col
)

check_required_columns(
  df = merged@meta.data,
  required_cols = required_cols,
  object_name = "baseline EBV metadata"
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
    celltype_for_de = clean_chr(.data[[config$celltype_col]]),
    lifecycle_for_de = clean_chr(.data[[config$lifecycle_col]])
  )

if (!identical(colnames(counts_all), meta_all$cell_barcode)) {
  stop("Counts columns and metadata cell barcodes are not aligned.")
}

message("Lifecycle values in object:")
print(table(meta_all$lifecycle_for_de, useNA = "ifany"))


# ----------------------------- #
# 4. edgeR function per lifecycle x celltype
# ----------------------------- #

run_edger_one_celltype <- function(meta_group, counts_group, target_label, celltype_label) {
  
  message("Running lifecycle: ", target_label, " | celltype: ", celltype_label)
  
  meta_ct_cells <- meta_group %>%
    filter(celltype_for_de == celltype_label)
  
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
  
  qc_row <- donors_per_group %>%
    mutate(
      target_lifecycle = target_label,
      celltype_new = celltype_label,
      n_cells = nrow(meta_ct_cells)
    ) %>%
    relocate(target_lifecycle, celltype_new, n_cells)
  
  if (
    is.na(n_control) || is.na(n_ms) ||
    n_control < config$min_donors_per_cohort ||
    n_ms < config$min_donors_per_cohort
  ) {
    message(
      "  Skipping ",
      target_label,
      " | ",
      celltype_label,
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
      qc = qc_row,
      skipped = tibble(
        target_lifecycle = target_label,
        celltype_new = celltype_label,
        reason = "insufficient_donors_per_cohort",
        n_Control = n_control,
        n_MS = n_ms,
        n_cells = nrow(meta_ct_cells)
      )
    ))
  }
  
  # Pseudobulk counts per donor within this lifecycle x cell type.
  donor_factor <- factor(meta_ct_cells$donor_id)
  
  aggregation_matrix <- Matrix::sparseMatrix(
    i = seq_along(donor_factor),
    j = as.integer(donor_factor),
    x = 1,
    dims = c(length(donor_factor), nlevels(donor_factor)),
    dimnames = list(meta_ct_cells$cell_barcode, levels(donor_factor))
  )
  
  pb_counts <- counts_group[, meta_ct_cells$cell_barcode, drop = FALSE] %*% aggregation_matrix
  pb_counts <- as.matrix(pb_counts)
  
  smeta <- meta_ct_cells %>%
    distinct(donor_id, cohort) %>%
    filter(donor_id %in% colnames(pb_counts)) %>%
    arrange(match(donor_id, colnames(pb_counts))) %>%
    mutate(
      cohort = factor(
        cohort,
        levels = c(config$control_label, config$ms_label)
      )
    )
  
  if (!identical(smeta$donor_id, colnames(pb_counts))) {
    stop(
      "Pseudobulk columns and sample metadata are not aligned for: ",
      target_label,
      " | ",
      celltype_label
    )
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
    message("  Skipping ", target_label, " | ", celltype_label, ": no genes retained after filterByExpr.")
    
    return(list(
      result = NULL,
      qc = qc_row,
      skipped = tibble(
        target_lifecycle = target_label,
        celltype_new = celltype_label,
        reason = "no_genes_retained_after_filterByExpr",
        n_Control = n_control,
        n_MS = n_ms,
        n_cells = nrow(meta_ct_cells)
      )
    ))
  }
  
  design <- model.matrix(
    ~ cohort,
    data = y$samples
  )
  
  if (qr(design)$rank < ncol(design)) {
    message("  Skipping ", target_label, " | ", celltype_label, ": design matrix is not full rank.")
    
    return(list(
      result = NULL,
      qc = qc_row,
      skipped = tibble(
        target_lifecycle = target_label,
        celltype_new = celltype_label,
        reason = "design_matrix_not_full_rank",
        n_Control = n_control,
        n_MS = n_ms,
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
    message("  Skipping ", target_label, " | ", celltype_label, ": coefficient ", coef_name, " not found.")
    
    return(list(
      result = NULL,
      qc = qc_row,
      skipped = tibble(
        target_lifecycle = target_label,
        celltype_new = celltype_label,
        reason = paste0("coefficient_not_found_", coef_name),
        n_Control = n_control,
        n_MS = n_ms,
        n_cells = nrow(meta_ct_cells)
      )
    ))
  }
  
  qlf <- glmQLFTest(
    fit,
    coef = which(colnames(design) == coef_name)
  )
  
  tt <- topTags(qlf, n = Inf)$table %>%
    as.data.frame() %>%
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
  
  results_ct <- tt %>%
    mutate(
      target_lifecycle = target_label,
      celltype_new = celltype_label,
      comparison = paste0(config$ms_label, "_vs_", config$control_label),
      n_donors_MS = n_ms,
      n_donors_Control = n_control,
      n_pseudobulk_samples = ncol(pb_counts),
      n_cells = nrow(meta_ct_cells),
      direction = assign_direction(
        logfc = avg_log2FC,
        fdr = p_val_adj,
        fdr_threshold = config$fdr_threshold,
        abs_logfc_threshold = config$abs_logfc_threshold_for_direction
      )
    ) %>%
    relocate(
      target_lifecycle,
      celltype_new,
      comparison,
      gene,
      n_donors_MS,
      n_donors_Control,
      n_pseudobulk_samples,
      n_cells,
      avg_log2FC
    ) %>%
    arrange(p_val_adj, p_val)
  
  lifecycle_dir <- file.path(
    config$output_dir,
    safe_name(target_label)
  )
  
  create_dir(lifecycle_dir)
  
  output_file <- file.path(
    lifecycle_dir,
    paste0(
      "edgeR_",
      safe_name(target_label),
      "_",
      config$ms_label,
      "_vs_",
      config$control_label,
      "_",
      safe_name(celltype_label),
      ".csv"
    )
  )
  
  write_csv(results_ct, output_file)
  
  list(
    result = results_ct,
    qc = qc_row,
    skipped = tibble()
  )
}


# ----------------------------- #
# 5. Run one lifecycle group
# ----------------------------- #

run_one_lifecycle_group <- function(target_label, target_pattern, use_regex, ignore_case) {
  
  message("\n========================================")
  message("Running lifecycle group: ", target_label)
  message("========================================")
  
  meta_group <- meta_all %>%
    mutate(
      target_group = match_lifecycle(
        lifecycle_for_de,
        pattern = target_pattern,
        use_regex = use_regex,
        ignore_case = ignore_case
      )
    ) %>%
    filter(
      target_group,
      cohort %in% c(config$control_label, config$ms_label),
      !is.na(donor_id), donor_id != "",
      !is.na(celltype_for_de), celltype_for_de != ""
    )
  
  if (nrow(meta_group) == 0) {
    warning("No cells remained for lifecycle group: ", target_label)
    
    return(list(
      result = NULL,
      qc = tibble(),
      skipped = tibble(
        target_lifecycle = target_label,
        celltype_new = NA_character_,
        reason = "no_cells_after_filtering",
        n_Control = NA_real_,
        n_MS = NA_real_,
        n_cells = 0
      ),
      filtered_metadata = tibble()
    ))
  }
  
  counts_group <- counts_all[, meta_group$cell_barcode, drop = FALSE]
  
  message("Cells retained for ", target_label, ":")
  print(table(meta_group$cohort, useNA = "ifany"))
  
  message("Cells by celltype for ", target_label, ":")
  print(table(meta_group$celltype_for_de, meta_group$cohort, useNA = "ifany"))
  
  if (is.null(config$celltype_order)) {
    celltypes_to_test <- meta_group %>%
      distinct(celltype_for_de) %>%
      arrange(celltype_for_de) %>%
      pull(celltype_for_de)
  } else {
    celltypes_to_test <- config$celltype_order[
      config$celltype_order %in% unique(meta_group$celltype_for_de)
    ]
    
    unexpected_celltypes <- setdiff(
      sort(unique(meta_group$celltype_for_de)),
      config$celltype_order
    )
    
    celltypes_to_test <- c(celltypes_to_test, unexpected_celltypes)
  }
  
  run_list <- lapply(
    celltypes_to_test,
    function(ct) {
      run_edger_one_celltype(
        meta_group = meta_group,
        counts_group = counts_group,
        target_label = target_label,
        celltype_label = ct
      )
    }
  )
  
  result_tbl <- bind_rows(
    lapply(run_list, `[[`, "result")
  )
  
  qc_tbl <- bind_rows(
    lapply(run_list, `[[`, "qc")
  )
  
  skipped_tbl <- bind_rows(
    lapply(run_list, `[[`, "skipped")
  )
  
  if (nrow(result_tbl) == 0) {
    result_tbl <- make_empty_result_table()
  }
  
  if (nrow(skipped_tbl) == 0) {
    skipped_tbl <- tibble(
      target_lifecycle = character(),
      celltype_new = character(),
      reason = character(),
      n_Control = numeric(),
      n_MS = numeric(),
      n_cells = numeric()
    )
  }
  
  lifecycle_dir <- file.path(
    config$output_dir,
    safe_name(target_label)
  )
  create_dir(lifecycle_dir)
  
  filtered_metadata_path <- file.path(
    lifecycle_dir,
    paste0("edgeR_", safe_name(target_label), "_filtered_metadata.csv")
  )
  
  all_results_path <- file.path(
    lifecycle_dir,
    paste0("edgeR_", safe_name(target_label), "_MS_vs_Control_ALL_by_celltype.csv")
  )
  
  significant_results_path <- file.path(
    lifecycle_dir,
    paste0("edgeR_", safe_name(target_label), "_MS_vs_Control_SIGNIFICANT_by_celltype.csv")
  )
  
  qc_path <- file.path(
    lifecycle_dir,
    paste0("edgeR_", safe_name(target_label), "_QC_by_celltype.csv")
  )
  
  skipped_path <- file.path(
    lifecycle_dir,
    paste0("edgeR_", safe_name(target_label), "_skipped_celltypes.csv")
  )
  
  significant_tbl <- result_tbl %>%
    filter(
      !is.na(p_val_adj),
      p_val_adj < config$fdr_threshold
    )
  
  write_csv(meta_group, filtered_metadata_path)
  write_csv(result_tbl, all_results_path)
  write_csv(significant_tbl, significant_results_path)
  write_csv(qc_tbl, qc_path)
  write_csv(skipped_tbl, skipped_path)
  
  list(
    result = result_tbl,
    qc = qc_tbl,
    skipped = skipped_tbl,
    filtered_metadata = meta_group
  )
}


# ----------------------------- #
# 6. Run all lifecycle groups
# ----------------------------- #

group_results <- pmap(
  list(
    target_label = config$target_groups$target_lifecycle_label,
    target_pattern = config$target_groups$target_lifecycle_pattern,
    use_regex = config$target_groups$use_regex,
    ignore_case = config$target_groups$ignore_case
  ),
  run_one_lifecycle_group
)

de_all <- bind_rows(
  lapply(group_results, `[[`, "result")
)

qc_by_celltype <- bind_rows(
  lapply(group_results, `[[`, "qc")
)

skipped_celltypes <- bind_rows(
  lapply(group_results, `[[`, "skipped")
)

filtered_metadata_all <- bind_rows(
  lapply(group_results, `[[`, "filtered_metadata")
)

if (nrow(de_all) == 0) {
  warning("No edgeR results were generated for any lifecycle/celltype.")
  de_all <- make_empty_result_table()
}

if (nrow(skipped_celltypes) == 0) {
  skipped_celltypes <- tibble(
    target_lifecycle = character(),
    celltype_new = character(),
    reason = character(),
    n_Control = numeric(),
    n_MS = numeric(),
    n_cells = numeric()
  )
}


# ----------------------------- #
# 7. Significant genes
# ----------------------------- #

de_sig <- de_all %>%
  filter(
    !is.na(p_val_adj),
    p_val_adj < config$fdr_threshold
  )

sig_summary <- de_sig %>%
  count(
    target_lifecycle,
    celltype_new,
    direction,
    name = "n_genes"
  ) %>%
  arrange(target_lifecycle, celltype_new, direction)


# ----------------------------- #
# 8. Save combined outputs
# ----------------------------- #

filtered_metadata_path <- file.path(
  config$output_dir,
  config$filtered_metadata_csv
)

qc_path <- file.path(
  config$output_dir,
  config$qc_by_celltype_csv
)

all_results_path <- file.path(
  config$output_dir,
  config$all_results_csv
)

sig_results_path <- file.path(
  config$output_dir,
  config$significant_results_csv
)

sig_summary_path <- file.path(
  config$output_dir,
  config$significant_summary_csv
)

skipped_path <- file.path(
  config$output_dir,
  config$skipped_celltypes_csv
)

write_csv(filtered_metadata_all, filtered_metadata_path)
write_csv(qc_by_celltype, qc_path)
write_csv(de_all, all_results_path)
write_csv(de_sig, sig_results_path)
write_csv(sig_summary, sig_summary_path)
write_csv(skipped_celltypes, skipped_path)


# ----------------------------- #
# 9. Session information
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
    "Lifecycle values in object:",
    capture.output(print(table(meta_all$lifecycle_for_de, useNA = "ifany"))),
    "",
    "Combined filtered metadata rows:",
    as.character(nrow(filtered_metadata_all)),
    "",
    "QC by lifecycle/celltype:",
    capture.output(print(qc_by_celltype)),
    "",
    "Skipped lifecycle/celltypes:",
    capture.output(print(skipped_celltypes)),
    "",
    "Significant gene summary:",
    capture.output(print(sig_summary)),
    "",
    "Output files:",
    filtered_metadata_path,
    qc_path,
    all_results_path,
    sig_results_path,
    sig_summary_path,
    skipped_path,
    "",
    "Session information:",
    capture.output(sessionInfo())
  ),
  session_info_path
)


# ----------------------------- #
# 10. Completion messages
# ----------------------------- #

message("\nedgeR MS vs Control by lifecycle and celltype workflow complete.")
message("Combined filtered metadata: ", filtered_metadata_path)
message("Combined QC by lifecycle/celltype: ", qc_path)
message("Combined all DE results: ", all_results_path)
message("Combined significant DE results: ", sig_results_path)
message("Combined significant summary: ", sig_summary_path)
message("Combined skipped lifecycle/celltypes: ", skipped_path)
message("Session info: ", session_info_path)

message("\nSignificant gene summary:")
print(sig_summary)
