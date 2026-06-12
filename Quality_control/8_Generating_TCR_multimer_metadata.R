#!/usr/bin/env Rscript

# ==============================================================================
# Multimer cutoff filtering and VDJ export from final merged Seurat object
# ==============================================================================
#
# Purpose:
#   This script loads the final merged Seurat object, rescans Cell Ranger VDJ
#   output folders to create merged contig-level TCR CSV files, visualises
#   tetramer CLR distributions, applies manually defined tetramer-specific CLR
#   cutoffs, and saves:
#     1. final merged object with cutoff-updated tetramer calls,
#     2. EBV CD8-only subset,
#     3. merged VDJ contig CSV for all matched cells,
#     4. merged VDJ contig CSV for EBV CD8-only cells,
#     5. multimer CLR violin plots,
#     6. cutoff and subset summaries,
#     7. session information.
#
# Expected input:
#   A final merged Seurat object from the previous integration step:
#     8_merged_all_objects_final.rds
#
# Expected metadata:
#   The Seurat object should contain:
#     - id
#     - tetramer
#     - tetramer_CLR
#     - sample and batch, either already present or joinable from sample metadata
#
# Notes:
#   - This script does not re-call tetramers from ADT.
#   - It applies manual CLR cutoffs to existing tetramer calls.
#   - Low-CLR tetramer calls are changed to "negative".
#   - Original calls are preserved in tetramer_raw.
#   - Keep project-specific paths in the CONFIG section only.
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
  library(readr)
  library(readxl)
  library(stringr)
  library(purrr)
  library(ggplot2)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  input_rds = "path/to/input/8_merged_all_objects_final.rds",
  
  output_dir = "path/to/output",
  plot_dir = "path/to/output/multimer_CLR_cutoff_plots",
  vdj_output_dir = "path/to/output/VDJ_analysis",
  
  output_cutoff_rds = "path/to/output/9_merged_all_objects_multimer_cutoff.rds",
  output_ebv_cd8_rds = "path/to/output/10_merged_EBV_CD8_samples_only.rds",
  
  merged_all_vdj_csv = "path/to/output/VDJ_analysis/merged_all_vdjs.csv",
  merged_ebv_cd8_vdj_csv = "path/to/output/VDJ_analysis/merged_ebv_cd8_vdjs.csv",
  
  cutoff_summary_csv = "path/to/output/multimer_cutoff_summary.csv",
  vdj_scan_summary_csv = "path/to/output/VDJ_scan_summary.csv",
  subset_summary_csv = "path/to/output/EBV_CD8_subset_summary.csv",
  session_info_file = "path/to/output/sessionInfo_multimer_cutoff_VDJ_export.txt",
  
  # Optional sample-level metadata.
  # If sample/batch already exist in the object, this can remain TRUE and it will
  # only fill/join metadata using id.
  use_sample_metadata = TRUE,
  sample_metadata_file = "path/to/input/Sample_information.xlsx",
  sample_id_column = "id",
  sample_column = "sample",
  batch_column = "batch",
  
  # Metadata columns used for multimer filtering.
  tetramer_column = "tetramer",
  tetramer_clr_column = "tetramer_CLR",
  tetramer_raw_column = "tetramer_raw",
  negative_label = "negative",
  
  # Fix old mislabelled tetramer names before plotting/cutoff.
  apply_tetramer_name_fix = TRUE,
  tetramer_fix_pattern = "^RPPK",
  tetramer_fix_replacement = "RPQK",
  
  # If a metadata column called "multimer" exists, use it to define multimer+
  # plotting cells. Otherwise, use tetramer != "negative".
  optional_multimer_status_column = "multimer",
  optional_multimer_positive_value = "positive",
  
  # Tetramer plotting.
  tetramer_order = c(
    "GLCT", "CLGG", "FLYA", "YLQQ", "LLDF", "GILG", "SSCS*", "AVFD*",
    "ATIG", "PYLF", "TYGP", "TYSA", "TYPV", "RYGF", "RPQG", "RPQK*",
    "RPPI", "LPRR", "YNLR*", "FLRG", "QAKW", "RAKF", "ELRS", "HPVG",
    "MGSL", "YPLH*", "EPLP*", "LPFE", "LEKA", "FENI", "IEDP", "VEDL", "QEIR"
  ),
  tetramers_per_plot_page = 6,
  plot_width = 25,
  plot_height = 20,
  plot_dpi = 300,
  
  # Manual tetramer CLR cutoffs.
  cutoffs = tibble::tribble(
    ~tetramer, ~cutoff,
    "GLCT",    2.5,
    "CLGG",    2.5,
    "FLYA",    3,
    "YLQQ",    2,
    "LLDF",    2.5,
    "GILG",    3,
    "SSCS*",   2.5,
    "AVFD*",   3,
    "ATIG",    2,
    "PYLF",    3,
    "TYGP",    2.5,
    "TYSA",    2,
    "TYPV",    4,
    "RYGF",    1,
    "RPQG",    3,
    "RPQK*",   2.5,
    "RPPI",    2.5,
    "LPRR",    3,
    "YNLR*",   3,
    "FLRG",    3,
    "QAKW",    3,
    "RAKF",    3,
    "ELRS",    4,
    "HPVG",    2,
    "MGSL",    2,
    "YPLH*",   2,
    "EPLP*",   3,
    "LPFE",    2,
    "LEKA",    2,
    "FENI",    2,
    "IEDP",    5,
    "VEDL",    4.5,
    "QEIR",    5
  ),
  
  # VDJ scanning configuration.
  # Each pattern should point to vdj_t folders.
  # Use * for sample folder, not {sample}, because this script scans all VDJ
  # folders and then keeps only samples present in the Seurat object.
  vdj_sources = list(
    list(
      source_id = "global_pbmc_batch1",
      base_dir = "path/to/vdj/base_main",
      patterns = c("GEM*/*/outs/per_sample_outs/*/vdj_t")
    ),
    list(
      source_id = "global_pbmc_batch2",
      base_dir = "path/to/vdj/New_global_batches",
      patterns = c(
        "*/outs/per_sample_outs/*/vdj_t",
        "*/*/outs/per_sample_outs/*/vdj_t"
      )
    ),
    list(
      source_id = "ebv_enrichment_batch1",
      base_dir = "path/to/vdj/EBV_data",
      patterns = c(
        "*/outs/per_sample_outs/*/vdj_t",
        "*/*/outs/per_sample_outs/*/vdj_t"
      )
    ),
    list(
      source_id = "ebv_enrichment_batch2",
      base_dir = "path/to/vdj/new_EBV_data",
      patterns = c(
        "*/outs/per_sample_outs/*/vdj_t",
        "*/*/outs/per_sample_outs/*/vdj_t"
      )
    ),
    list(
      source_id = "ebv_activated",
      base_dir = "path/to/vdj/LCL_data",
      patterns = c(
        "*/outs/per_sample_outs/*/vdj_t",
        "*/*/outs/per_sample_outs/*/vdj_t"
      )
    )
  ),
  
  vdj_contig_file = "filtered_contig_annotations.csv",
  
  # EBV CD8 subset settings.
  ebv_batch_pattern = "GEMEBV",
  exclude_batch_pattern = "CD4",
  
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


load_seurat_object <- function(input_rds) {
  if (!file.exists(input_rds)) {
    stop("Input RDS file does not exist: ", input_rds)
  }
  
  obj <- readRDS(input_rds)
  
  if (!inherits(obj, "Seurat")) {
    stop("Input RDS must contain a Seurat object.")
  }
  
  obj
}


normalise_tetramer_names <- function(x, cfg) {
  x <- as.character(x)
  
  if (isTRUE(cfg$apply_tetramer_name_fix)) {
    x <- ifelse(
      grepl(cfg$tetramer_fix_pattern, x),
      sub(cfg$tetramer_fix_pattern, cfg$tetramer_fix_replacement, x),
      x
    )
  }
  
  x
}


ensure_required_metadata <- function(obj, cfg) {
  metadata_cols <- colnames(obj@meta.data)
  
  required_cols <- c(
    cfg$sample_id_column,
    cfg$tetramer_column,
    cfg$tetramer_clr_column
  )
  
  missing_cols <- setdiff(required_cols, metadata_cols)
  
  if (length(missing_cols) > 0) {
    stop(
      "Merged object is missing required metadata column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  obj[[cfg$tetramer_column]] <- normalise_tetramer_names(
    obj[[cfg$tetramer_column, drop = TRUE]],
    cfg
  )
  
  obj[[cfg$tetramer_clr_column]] <- as.numeric(
    obj[[cfg$tetramer_clr_column, drop = TRUE]]
  )
  
  obj
}


join_sample_metadata_if_needed <- function(obj, cfg) {
  if (!isTRUE(cfg$use_sample_metadata)) {
    return(obj)
  }
  
  if (!file.exists(cfg$sample_metadata_file)) {
    stop("Sample metadata file does not exist: ", cfg$sample_metadata_file)
  }
  
  sample_info <- read_xlsx(cfg$sample_metadata_file) |>
    as.data.frame()
  
  if (!cfg$sample_id_column %in% colnames(sample_info)) {
    stop(
      "Sample metadata file must contain column: ",
      cfg$sample_id_column
    )
  }
  
  sample_info_unique <- sample_info |>
    group_by(.data[[cfg$sample_id_column]]) |>
    slice(1) |>
    ungroup()
  
  meta <- obj@meta.data |>
    rownames_to_column("cell")
  
  # Drop overlapping sample-info columns except the join key, so the sample sheet
  # becomes the clean source of truth.
  overlapping_cols <- intersect(colnames(meta), colnames(sample_info_unique))
  overlapping_cols <- setdiff(overlapping_cols, c("cell", cfg$sample_id_column))
  
  if (length(overlapping_cols) > 0) {
    meta <- meta |>
      select(-all_of(overlapping_cols))
  }
  
  meta_joined <- meta |>
    left_join(sample_info_unique, by = cfg$sample_id_column) |>
    column_to_rownames("cell")
  
  obj@meta.data <- meta_joined
  
  obj
}


extract_batch_from_path <- function(vdj_dir, base_dir) {
  vdj_dir <- normalizePath(vdj_dir, winslash = "/", mustWork = FALSE)
  base_dir <- normalizePath(base_dir, winslash = "/", mustWork = FALSE)
  
  relative_path <- sub(paste0("^", base_dir, "/?"), "", vdj_dir)
  batch <- strsplit(relative_path, "/", fixed = TRUE)[[1]][1]
  
  if (is.na(batch) || batch == "") {
    batch <- NA_character_
  }
  
  batch
}


find_vdj_dirs_for_source <- function(source_cfg) {
  vdj_dirs <- unlist(
    lapply(
      source_cfg$patterns,
      function(pattern) {
        Sys.glob(file.path(source_cfg$base_dir, pattern))
      }
    ),
    use.names = FALSE
  )
  
  unique(vdj_dirs[dir.exists(vdj_dirs)])
}


scan_one_vdj_source <- function(source_cfg, sample_names, cfg) {
  message("Scanning VDJ source: ", source_cfg$source_id)
  
  vdj_dirs <- find_vdj_dirs_for_source(source_cfg)
  
  if (length(vdj_dirs) == 0) {
    return(list(
      vdj_data = tibble(),
      summary = tibble(
        source_id = source_cfg$source_id,
        base_dir = source_cfg$base_dir,
        status = "no_vdj_dirs_found",
        sample = NA_character_,
        vdj_dir = NA_character_,
        n_rows = 0
      )
    ))
  }
  
  vdj_entries <- list()
  summaries <- list()
  
  for (vdj_dir in vdj_dirs) {
    sample_name <- basename(dirname(vdj_dir))
    
    if (!sample_name %in% sample_names) {
      next
    }
    
    vdj_file <- file.path(vdj_dir, cfg$vdj_contig_file)
    
    if (!file.exists(vdj_file)) {
      summaries[[vdj_dir]] <- tibble(
        source_id = source_cfg$source_id,
        base_dir = source_cfg$base_dir,
        status = "missing_contig_file",
        sample = sample_name,
        vdj_dir = vdj_dir,
        n_rows = 0
      )
      next
    }
    
    dat <- tryCatch(
      read_csv(vdj_file, show_col_types = FALSE),
      error = function(e) NULL
    )
    
    if (is.null(dat) || nrow(dat) == 0) {
      summaries[[vdj_dir]] <- tibble(
        source_id = source_cfg$source_id,
        base_dir = source_cfg$base_dir,
        status = "empty_or_failed_read",
        sample = sample_name,
        vdj_dir = vdj_dir,
        n_rows = 0
      )
      next
    }
    
    if (!"barcode" %in% colnames(dat)) {
      summaries[[vdj_dir]] <- tibble(
        source_id = source_cfg$source_id,
        base_dir = source_cfg$base_dir,
        status = "missing_barcode_column",
        sample = sample_name,
        vdj_dir = vdj_dir,
        n_rows = nrow(dat)
      )
      next
    }
    
    batch <- extract_batch_from_path(vdj_dir, source_cfg$base_dir)
    
    dat <- dat |>
      mutate(
        barcode = paste0(sample_name, "_", .data$barcode),
        id = sample_name,
        batch = batch,
        source_id = source_cfg$source_id,
        source_base = source_cfg$base_dir,
        vdj_dir = vdj_dir
      )
    
    key <- paste(source_cfg$source_id, sample_name, vdj_dir, sep = "||")
    vdj_entries[[key]] <- dat
    
    summaries[[key]] <- tibble(
      source_id = source_cfg$source_id,
      base_dir = source_cfg$base_dir,
      status = "processed",
      sample = sample_name,
      vdj_dir = vdj_dir,
      n_rows = nrow(dat)
    )
    
    message(
      "  Added ", nrow(dat),
      " rows [sample: ", sample_name,
      ", batch: ", batch, "]"
    )
  }
  
  list(
    vdj_data = bind_rows(vdj_entries),
    summary = bind_rows(summaries)
  )
}


scan_all_vdj_sources <- function(obj, cfg) {
  sample_names <- unique(obj@meta.data[[cfg$sample_id_column]])
  sample_names <- sample_names[!is.na(sample_names)]
  
  results <- lapply(
    cfg$vdj_sources,
    scan_one_vdj_source,
    sample_names = sample_names,
    cfg = cfg
  )
  
  vdj_data <- bind_rows(lapply(results, `[[`, "vdj_data"))
  vdj_summary <- bind_rows(lapply(results, `[[`, "summary"))
  
  list(
    vdj_data = vdj_data,
    vdj_summary = vdj_summary
  )
}


make_plot_dataframe <- function(obj, cfg) {
  meta <- obj@meta.data |>
    rownames_to_column("cell")
  
  if (!cfg$sample_column %in% colnames(meta)) {
    meta[[cfg$sample_column]] <- meta[[cfg$sample_id_column]]
  }
  
  if (all(is.na(meta[[cfg$sample_column]]))) {
    meta[[cfg$sample_column]] <- meta[[cfg$sample_id_column]]
  }
  
  meta <- meta |>
    mutate(
      tetramer = normalise_tetramer_names(.data[[cfg$tetramer_column]], cfg),
      tetramer_CLR = as.numeric(.data[[cfg$tetramer_clr_column]]),
      sample = .data[[cfg$sample_column]]
    )
  
  if (cfg$optional_multimer_status_column %in% colnames(meta)) {
    meta <- meta |>
      filter(.data[[cfg$optional_multimer_status_column]] == cfg$optional_multimer_positive_value)
  } else {
    meta <- meta |>
      filter(tetramer != cfg$negative_label)
  }
  
  meta |>
    filter(
      !is.na(tetramer),
      !is.na(tetramer_CLR),
      !is.na(sample),
      tetramer %in% cfg$tetramer_order
    ) |>
    mutate(
      tetramer = factor(tetramer, levels = cfg$tetramer_order)
    )
}


plot_tetramer_clr_pages <- function(plot_df, cfg) {
  create_dir(cfg$plot_dir)
  
  if (nrow(plot_df) == 0) {
    warning("No multimer-positive cells available for plotting.")
    return(invisible(NULL))
  }
  
  tetramer_chunks <- split(
    cfg$tetramer_order,
    ceiling(seq_along(cfg$tetramer_order) / cfg$tetramers_per_plot_page)
  )
  
  for (i in seq_along(tetramer_chunks)) {
    tetramer_group <- tetramer_chunks[[i]]
    
    subset_df <- plot_df |>
      filter(tetramer %in% tetramer_group) |>
      mutate(tetramer = factor(tetramer, levels = tetramer_group))
    
    if (nrow(subset_df) == 0) {
      next
    }
    
    p <- ggplot(
      subset_df,
      aes(x = .data$sample, y = .data$tetramer_CLR, fill = .data$sample)
    ) +
      geom_violin(scale = "width", trim = FALSE) +
      geom_jitter(width = 0.2, size = 0.6, alpha = 0.6) +
      facet_wrap(~ tetramer, scales = "free_y") +
      theme_minimal(base_size = 16) +
      labs(
        title = paste("Multimer-positive tetramer CLR expression: page", i),
        x = "Sample",
        y = "Tetramer CLR"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.position = "none"
      )
    
    ggsave(
      filename = file.path(
        cfg$plot_dir,
        paste0("multimer_positive_tetramer_CLR_page_", i, ".jpeg")
      ),
      plot = p,
      width = cfg$plot_width,
      height = cfg$plot_height,
      dpi = cfg$plot_dpi
    )
  }
  
  invisible(NULL)
}


apply_tetramer_cutoffs <- function(obj, cfg) {
  meta <- obj@meta.data |>
    rownames_to_column("cell") |>
    mutate(
      tetramer = normalise_tetramer_names(.data[[cfg$tetramer_column]], cfg),
      tetramer_CLR = as.numeric(.data[[cfg$tetramer_clr_column]])
    )
  
  if (!cfg$tetramer_raw_column %in% colnames(meta)) {
    meta[[cfg$tetramer_raw_column]] <- meta$tetramer
  }
  
  before_counts <- meta |>
    count(tetramer, name = "n_before")
  
  meta_cut <- meta |>
    left_join(cfg$cutoffs, by = "tetramer") |>
    mutate(
      tetramer_pass_cutoff = case_when(
        tetramer == cfg$negative_label ~ FALSE,
        is.na(cutoff) ~ TRUE,
        is.na(tetramer_CLR) ~ FALSE,
        tetramer_CLR >= cutoff ~ TRUE,
        TRUE ~ FALSE
      ),
      tetramer = ifelse(
        tetramer != cfg$negative_label &
          !is.na(cutoff) &
          !is.na(tetramer_CLR) &
          tetramer_CLR < cutoff,
        cfg$negative_label,
        tetramer
      )
    )
  
  after_counts <- meta_cut |>
    count(tetramer, name = "n_after")
  
  cutoff_summary <- before_counts |>
    full_join(after_counts, by = "tetramer") |>
    mutate(
      n_before = replace_na(n_before, 0L),
      n_after = replace_na(n_after, 0L),
      n_removed_by_cutoff = n_before - n_after
    ) |>
    arrange(match(tetramer, c(cfg$tetramer_order, cfg$negative_label)))
  
  meta_cut <- meta_cut |>
    select(-cutoff) |>
    column_to_rownames("cell")
  
  obj@meta.data <- meta_cut
  
  list(
    obj = obj,
    cutoff_summary = cutoff_summary
  )
}


subset_ebv_cd8_cells <- function(obj, cfg) {
  meta <- obj@meta.data |>
    rownames_to_column("cell")
  
  required_cols <- c(cfg$sample_column, cfg$batch_column)
  missing_cols <- setdiff(required_cols, colnames(meta))
  
  if (length(missing_cols) > 0) {
    stop(
      "Cannot create EBV CD8 subset. Missing metadata column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  keep_cells <- meta |>
    filter(
      !is.na(.data[[cfg$batch_column]]),
      grepl(cfg$ebv_batch_pattern, .data[[cfg$batch_column]]),
      !grepl(cfg$exclude_batch_pattern, .data[[cfg$batch_column]], ignore.case = TRUE)
    ) |>
    pull(cell)
  
  subset_obj <- subset(obj, cells = keep_cells)
  
  subset_summary <- tibble(
    n_cells_before_subset = ncol(obj),
    n_cells_after_subset = ncol(subset_obj),
    ebv_batch_pattern = cfg$ebv_batch_pattern,
    exclude_batch_pattern = cfg$exclude_batch_pattern,
    n_samples_after_subset = n_distinct(subset_obj@meta.data[[cfg$sample_column]]),
    n_batches_after_subset = n_distinct(subset_obj@meta.data[[cfg$batch_column]])
  )
  
  list(
    subset_obj = subset_obj,
    subset_summary = subset_summary
  )
}


# ----------------------------- #
# 3. Run pipeline
# ----------------------------- #

create_dir(config$output_dir)
create_dir(config$plot_dir)
create_dir(config$vdj_output_dir)
create_parent_dir(config$output_cutoff_rds)
create_parent_dir(config$output_ebv_cd8_rds)
create_parent_dir(config$merged_all_vdj_csv)
create_parent_dir(config$merged_ebv_cd8_vdj_csv)
create_parent_dir(config$cutoff_summary_csv)
create_parent_dir(config$vdj_scan_summary_csv)
create_parent_dir(config$subset_summary_csv)
create_parent_dir(config$session_info_file)

message("Loading final merged object...")
merged_obj <- load_seurat_object(config$input_rds)

message("Checking required metadata...")
merged_obj <- ensure_required_metadata(merged_obj, config)

message("Joining sample metadata if configured...")
merged_obj <- join_sample_metadata_if_needed(merged_obj, config)

message("Scanning VDJ folders and creating merged VDJ table...")
vdj_scan <- scan_all_vdj_sources(merged_obj, config)

final_merged_vdj <- vdj_scan$vdj_data
vdj_scan_summary <- vdj_scan$vdj_summary

write.csv(
  final_merged_vdj,
  config$merged_all_vdj_csv,
  row.names = FALSE
)

write.csv(
  vdj_scan_summary,
  config$vdj_scan_summary_csv,
  row.names = FALSE
)

message("Creating tetramer CLR plots...")
plot_df <- make_plot_dataframe(merged_obj, config)
plot_tetramer_clr_pages(plot_df, config)

message("Applying tetramer CLR cutoffs...")
cutoff_result <- apply_tetramer_cutoffs(merged_obj, config)

merged_obj <- cutoff_result$obj
cutoff_summary <- cutoff_result$cutoff_summary

write.csv(
  cutoff_summary,
  config$cutoff_summary_csv,
  row.names = FALSE
)

message("Creating EBV CD8-only subset...")
subset_result <- subset_ebv_cd8_cells(merged_obj, config)

merged_ebv_cd8 <- subset_result$subset_obj
subset_summary <- subset_result$subset_summary

write.csv(
  subset_summary,
  config$subset_summary_csv,
  row.names = FALSE
)

message("Exporting EBV CD8-only VDJ table...")
cells_keep <- colnames(merged_ebv_cd8)

vdj_ebv_cd8_subset <- final_merged_vdj |>
  filter(barcode %in% cells_keep)

write.csv(
  vdj_ebv_cd8_subset,
  config$merged_ebv_cd8_vdj_csv,
  row.names = FALSE
)

message("Saving RDS outputs...")
saveRDS(
  merged_obj,
  config$output_cutoff_rds
)

saveRDS(
  merged_ebv_cd8,
  config$output_ebv_cd8_rds
)

writeLines(
  capture.output(sessionInfo()),
  config$session_info_file
)

message("\nMultimer cutoff filtering and VDJ export complete.")
message("Saved cutoff-updated merged object to: ", config$output_cutoff_rds)
message("Saved EBV CD8-only object to: ", config$output_ebv_cd8_rds)
message("Saved all VDJ CSV to: ", config$merged_all_vdj_csv)
message("Saved EBV CD8 VDJ CSV to: ", config$merged_ebv_cd8_vdj_csv)
message("Saved cutoff summary to: ", config$cutoff_summary_csv)
message("Saved VDJ scan summary to: ", config$vdj_scan_summary_csv)
message("Saved subset summary to: ", config$subset_summary_csv)
message("Saved session info to: ", config$session_info_file)