#!/usr/bin/env Rscript

# ==============================================================================
# TCR quality control, multimer refinement, and publication object export
# ==============================================================================
#
# Purpose:
#   This script performs TCR-level quality control and multimer refinement for
#   EBV-specific CD8+ T cell analyses. It:
#
#     1. collapses Cell Ranger VDJ contig rows to one row per cell,
#     2. resolves cells with multiple TRA/TRB chains using UMI dominance,
#     3. joins curated Seurat metadata,
#     4. performs within-sample clonotype-level tetramer QC,
#     5. optionally applies manual VDJdb-resolved tetramer calls,
#     6. performs global clonotype-level tetramer QC across samples,
#     7. optionally applies manual global VDJdb-resolved tetramer calls,
#     8. removes HLA-mismatched tetramer calls,
#     9. projects baseline GEM tetramer calls onto activated non-GEM cells,
#    10. adds curated TCR and tetramer metadata back to the Seurat object,
#    11. annotates virus/lifecycle/latency/antigen information,
#    12. gates global CD8+ cells using ADT CD3/CD4/CD8 expression,
#    13. saves final combined, baseline-only, and activated-paired objects.
#
# Expected inputs:
#   - Final EBV CD8-only Seurat object from the previous step:
#       10_merged_EBV_CD8_samples_only.rds
#   - Merged VDJ contig CSV:
#       merged_ebv_cd8_vdjs.csv
#   - Optional manual VDJdb-resolution sheets:
#       Conflicting_clones.xlsx
#       Global_conflicting_clones.xlsx
#   - HLA and tetramer reference sheets.
#
# Notes:
#   - Manual VDJdb calls are optional but recommended for final manuscript outputs.
#   - Original tetramer calls are preserved where possible.
#   - Keep project-specific paths in the CONFIG section only.
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
  library(readxl)
  library(purrr)
  library(ggplot2)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  input_seurat_rds = "path/to/input/10_merged_EBV_CD8_samples_only.rds",
  input_vdj_csv = "path/to/input/VDJ_analysis/merged_ebv_cd8_vdjs.csv",
  
  output_dir = "path/to/output/TCR_QC",
  output_rds_dir = "path/to/output/RDS",
  output_publication_dir = "path/to/output/publication_RDS",
  
  # Main RDS outputs.
  output_vdj_seurat_rds = "path/to/output/RDS/11_merged_VDJ_EBV_CD8_samples_only.rds",
  output_annotated_seurat_rds = "path/to/output/RDS/12_merged_VDJ_EBV_CD8_samples_annotated.rds",
  output_combined_publication_rds = "path/to/output/publication_RDS/13_EBV_CD8_combined_TCR_publication.rds",
  output_baseline_publication_rds = "path/to/output/publication_RDS/13_EBV_baseline_TCR_publication.rds",
  output_activated_publication_rds = "path/to/output/publication_RDS/13_EBV_LCL_activated_TCR_publication.rds",
  
  # CSV outputs.
  collapsed_multi_trb_csv = "path/to/output/TCR_QC/collapsed_CD8_EBV_cells_with_multiple_TRBs_labeled.csv",
  tcr_gem_only_csv = "path/to/output/TCR_QC/EBV_TCR_filtered_GEM_only.csv",
  within_sample_conflicts_csv = "path/to/output/TCR_QC/conflicting_clones_WITHIN_SAMPLE.csv",
  within_sample_auto_calls_csv = "path/to/output/TCR_QC/conflicting_clones_WITHIN_SAMPLE_AUTO.csv",
  tcr_qc1_csv = "path/to/output/TCR_QC/EBV_TCR_filtered_GEM_only_TETRAMER_QC1.csv",
  global_conflicts_csv = "path/to/output/TCR_QC/global_conflicting_TCRs.csv",
  global_auto_calls_csv = "path/to/output/TCR_QC/QC2_global_clone_calls_highest_or_vdjdbflag.csv",
  tcr_qc2_csv = "path/to/output/TCR_QC/EBV_TCR_filtered_GEM_only_TETRAMER_FINAL_QC2.csv",
  hla_removed_cells_csv = "path/to/output/TCR_QC/HLA_mismatched_cells_removed.csv",
  projected_all_tcr_csv = "path/to/output/TCR_QC/EBV_TCR_baseline_and_activated_projected.csv",
  qc_summary_csv = "path/to/output/TCR_QC/TCR_QC_summary.csv",
  session_info_file = "path/to/output/TCR_QC/sessionInfo_TCR_QC.txt",
  
  # Optional manual VDJdb-resolved files.
  use_manual_within_sample_calls = TRUE,
  manual_within_sample_xlsx = "path/to/input/Conflicting_clones.xlsx",
  within_sample_manual_final_col = "tetramer_final",
  
  use_manual_global_calls = TRUE,
  manual_global_xlsx = "path/to/input/Global_conflicting_clones.xlsx",
  global_manual_final_col = "new_tetramer",
  
  # If TRUE, missing manual files stop the script.
  # For public GitHub example scripts, FALSE is easier.
  # For final manuscript reproduction, set this to TRUE.
  manual_files_required = FALSE,
  
  # Reference files.
  tetramer_hla_xlsx = "path/to/input/Tetramer_HLA_match.xlsx",
  hla_information_xlsx = "path/to/input/HLA_information.xlsx",
  multimer_information_xlsx = "path/to/input/EBV_multimer_information.xlsx",
  
  # Metadata columns.
  barcode_col = "barcode",
  sample_col = "sample",
  id_col = "id",
  batch_col = "batch",
  cohort_col = "cohort",
  tetramer_col = "tetramer",
  tetramer_clr_col = "tetramer_CLR",
  tetramer_raw_col = "tetramer_raw",
  
  # Chain handling.
  chains_keep = c("TRA", "TRB"),
  tcr_columns = c("cdr3", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt"),
  tcr_numeric_columns = c("umis", "reads"),
  
  # Batch logic.
  baseline_batch_pattern = "^GEM",
  activated_batch_pattern = "^EBVLCL",
  ebv_enriched_batch_pattern = "GEMEBV|DEX",
  lcl_responsive_batch_pattern = "CD8",
  cd4_batch_exclude = c("EBVLCL1CD4", "EBVLCL2CD4"),
  
  # HLA filtering.
  apply_hla_filter = TRUE,
  negative_label = "negative",
  
  # Activated-cell projection.
  activated_enriched_label = "enriched",
  activated_lcl_responsive_label = "LCL-responsive",
  
  # ADT CD8 gating for global PBMC/global GEM cells.
  perform_adt_cd8_gating = TRUE,
  adt_assay = "ADT",
  adt_cd3_feature = "CD3",
  adt_cd4_feature = "CD4.1",
  adt_cd8_feature = "CD8",
  adt_cd3_min = 2.5,
  adt_cd4_max = 3,
  adt_cd8_min = 2,
  
  # Metadata cleanup before final save.
  metadata_columns_to_remove = c(
    "sample.x", "sample_id", "seurat_clusters", "clonotype_id",
    "is_cell", "contig_id", "high_confidence", "length", "multimer",
    "chain", "v_gene", "d_gene", "j_gene", "c_gene", "full_length",
    "productive", "fwr1", "fwr1_nt", "cdr1", "cdr1_nt",
    "fwr2", "fwr2_nt", "cdr2", "cdr2_nt",
    "fwr3", "fwr3_nt", "cdr3", "cdr3_nt",
    "fwr4", "fwr4_nt", "reads", "umis",
    "raw_consensus_id", "exact_subclonotype_id", "cdr3s_aa",
    "sample.y"
  ),
  
  verbose = TRUE
)


# ----------------------------- #
# 2. General helper functions
# ----------------------------- #

create_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

create_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

setup_output_dirs <- function(cfg) {
  create_dir(cfg$output_dir)
  create_dir(cfg$output_rds_dir)
  create_dir(cfg$output_publication_dir)
  
  output_paths <- cfg[names(cfg)[str_detect(names(cfg), "csv$|rds$|file$")]]
  walk(output_paths, create_parent_dir)
}

load_seurat_object <- function(path) {
  if (!file.exists(path)) {
    stop("Input Seurat RDS does not exist: ", path)
  }
  
  obj <- readRDS(path)
  
  if (!inherits(obj, "Seurat")) {
    stop("Input RDS must contain a Seurat object.")
  }
  
  obj
}

read_csv_required <- function(path) {
  if (!file.exists(path)) {
    stop("Input CSV does not exist: ", path)
  }
  
  read_csv(path, show_col_types = FALSE)
}

read_xlsx_optional <- function(path, required = FALSE) {
  if (is.null(path) || is.na(path) || path == "") {
    if (isTRUE(required)) stop("Required Excel path is empty.")
    return(NULL)
  }
  
  if (!file.exists(path)) {
    if (isTRUE(required)) stop("Required Excel file does not exist: ", path)
    warning("Optional Excel file not found: ", path)
    return(NULL)
  }
  
  read_xlsx(path)
}

get_assay_data_compat <- function(seurat_obj, assay, layer_or_slot = "data") {
  tryCatch(
    GetAssayData(seurat_obj, assay = assay, layer = layer_or_slot),
    error = function(e) {
      GetAssayData(seurat_obj, assay = assay, slot = layer_or_slot)
    }
  )
}

make_clone_key <- function(df, group_col = "sample") {
  required <- c(group_col, "TRB_v_gene", "TRB_j_gene", "TRB_cdr3")
  missing <- setdiff(required, colnames(df))
  
  if (length(missing) > 0) {
    stop("Missing columns for clone_key: ", paste(missing, collapse = ", "))
  }
  
  df |>
    mutate(
      clone_key = paste(
        .data[[group_col]],
        TRB_v_gene,
        TRB_j_gene,
        TRB_cdr3,
        sep = "|"
      )
    )
}

make_clone_key_global <- function(df) {
  required <- c("TRB_v_gene", "TRB_j_gene", "TRB_cdr3")
  missing <- setdiff(required, colnames(df))
  
  if (length(missing) > 0) {
    stop("Missing columns for clone_key_global: ", paste(missing, collapse = ", "))
  }
  
  df |>
    mutate(
      clone_key_global = ifelse(
        !is.na(TRB_cdr3) & TRB_cdr3 != "" &
          !is.na(TRB_v_gene) & TRB_v_gene != "" &
          !is.na(TRB_j_gene) & TRB_j_gene != "",
        paste(TRB_v_gene, TRB_j_gene, TRB_cdr3, sep = "|"),
        NA_character_
      )
    )
}


# ----------------------------- #
# 3. Collapse and clean VDJ chains
# ----------------------------- #

ensure_chain_slot_columns <- function(df, chains = c("TRA", "TRB")) {
  fields <- c("cdr3", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "umis", "reads")
  
  for (chain in chains) {
    for (slot in 1:2) {
      for (field in fields) {
        col_name <- paste0(chain, slot, "_", field)
        
        if (!col_name %in% colnames(df)) {
          if (field %in% c("umis", "reads")) {
            df[[col_name]] <- NA_real_
          } else {
            df[[col_name]] <- NA_character_
          }
        }
      }
    }
  }
  
  df
}

collapse_vdj_to_cell_rows <- function(vdj_df, cfg) {
  tcr_columns <- cfg$tcr_columns
  numeric_columns <- cfg$tcr_numeric_columns
  
  required_cols <- c("barcode", "chain")
  missing_required <- setdiff(required_cols, colnames(vdj_df))
  
  if (length(missing_required) > 0) {
    stop("VDJ CSV is missing required column(s): ", paste(missing_required, collapse = ", "))
  }
  
  existing_tcr_cols <- intersect(tcr_columns, colnames(vdj_df))
  existing_numeric_cols <- intersect(numeric_columns, colnames(vdj_df))
  
  meta_columns <- setdiff(
    colnames(vdj_df),
    c(tcr_columns, numeric_columns, "chain", "barcode", "contig_id", "length", "productive")
  )
  
  df_labeled <- vdj_df |>
    filter(chain %in% cfg$chains_keep) |>
    mutate(
      across(all_of(existing_tcr_cols), as.character),
      across(all_of(existing_numeric_cols), as.numeric)
    ) |>
    group_by(barcode, chain) |>
    arrange(desc(umis), desc(reads), .by_group = TRUE) |>
    mutate(
      chain_index = row_number(),
      chain_label = paste0(chain, chain_index)
    ) |>
    ungroup()
  
  df_wide <- df_labeled |>
    select(barcode, chain_label, all_of(existing_tcr_cols), all_of(existing_numeric_cols)) |>
    pivot_wider(
      names_from = chain_label,
      values_from = c(all_of(existing_tcr_cols), all_of(existing_numeric_cols)),
      names_glue = "{chain_label}_{.value}"
    )
  
  meta_df <- vdj_df |>
    select(barcode, all_of(meta_columns)) |>
    distinct(barcode, .keep_all = TRUE)
  
  df_collapsed <- df_wide |>
    left_join(meta_df, by = "barcode") |>
    ensure_chain_slot_columns(chains = cfg$chains_keep)
  
  df_collapsed
}

resolve_double_chain_by_umi <- function(df, chain) {
  slot1_prefix <- paste0(chain, "1_")
  slot2_prefix <- paste0(chain, "2_")
  
  fields <- c("cdr3", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "umis", "reads")
  
  slot1_cols <- paste0(slot1_prefix, fields)
  slot2_cols <- paste0(slot2_prefix, fields)
  
  cdr3_1 <- paste0(chain, "1_cdr3")
  cdr3_2 <- paste0(chain, "2_cdr3")
  umis_1 <- paste0(chain, "1_umis")
  umis_2 <- paste0(chain, "2_umis")
  
  multi_targets <- df |>
    filter(!is.na(.data[[cdr3_1]]) & !is.na(.data[[cdr3_2]])) |>
    select(barcode, all_of(c(umis_1, umis_2)))
  
  drop_barcodes <- multi_targets |>
    filter(
      is.na(.data[[umis_1]]) |
        is.na(.data[[umis_2]]) |
        .data[[umis_1]] == .data[[umis_2]]
    ) |>
    pull(barcode)
  
  best_chain <- multi_targets |>
    filter(!(barcode %in% drop_barcodes)) |>
    mutate(
      keep_slot = if_else(.data[[umis_1]] > .data[[umis_2]], 1L, 2L)
    ) |>
    select(barcode, keep_slot)
  
  df_resolved <- df |>
    left_join(best_chain, by = "barcode")
  
  for (i in seq_len(nrow(best_chain))) {
    bc <- best_chain$barcode[[i]]
    keep_slot <- best_chain$keep_slot[[i]]
    
    idx <- which(df_resolved$barcode == bc)
    if (length(idx) == 0) next
    
    if (keep_slot == 1L) {
      df_resolved[idx, slot2_cols] <- NA
    } else {
      df_resolved[idx, slot1_cols] <- df_resolved[idx, slot2_cols]
      df_resolved[idx, slot2_cols] <- NA
    }
  }
  
  df_resolved <- df_resolved |>
    select(-any_of("keep_slot")) |>
    filter(!(barcode %in% drop_barcodes))
  
  summary <- tibble(
    chain = chain,
    n_cells_with_two_chains = nrow(multi_targets),
    n_cells_dropped_equal_or_missing_umi = length(drop_barcodes),
    n_cells_resolved_by_umi = nrow(best_chain)
  )
  
  list(
    df = df_resolved,
    summary = summary,
    dropped_barcodes = drop_barcodes
  )
}

collapse_resolved_chain_columns <- function(df) {
  fields <- c("cdr3", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "umis", "reads")
  
  for (chain in c("TRB", "TRA")) {
    for (field in fields) {
      out_col <- paste0(chain, "_", field)
      col1 <- paste0(chain, "1_", field)
      col2 <- paste0(chain, "2_", field)
      
      df[[out_col]] <- coalesce(df[[col1]], df[[col2]])
    }
  }
  
  df |>
    select(
      -starts_with("TRB1_"),
      -starts_with("TRB2_"),
      -starts_with("TRA1_"),
      -starts_with("TRA2_")
    ) |>
    filter(!is.na(TRB_cdr3), TRB_cdr3 != "")
}

prepare_clean_tcr_table <- function(vdj_df, cfg) {
  df_collapsed <- collapse_vdj_to_cell_rows(vdj_df, cfg)
  
  multi_trb_labeled <- df_collapsed |>
    filter(!is.na(TRB1_cdr3) & !is.na(TRB2_cdr3))
  
  write.csv(
    multi_trb_labeled,
    cfg$collapsed_multi_trb_csv,
    row.names = FALSE
  )
  
  trb_resolved <- resolve_double_chain_by_umi(df_collapsed, chain = "TRB")
  
  df_trb <- trb_resolved$df |>
    filter(!is.na(TRB1_cdr3) | !is.na(TRB2_cdr3))
  
  tra_resolved <- resolve_double_chain_by_umi(df_trb, chain = "TRA")
  
  df_final <- collapse_resolved_chain_columns(tra_resolved$df)
  
  # Remove old sample column if it came from VDJ metadata and will be replaced.
  df_final <- df_final |>
    select(-any_of("sample"))
  
  summary <- bind_rows(trb_resolved$summary, tra_resolved$summary)
  
  list(
    df = df_final,
    chain_summary = summary
  )
}


# ----------------------------- #
# 4. Join Seurat metadata to TCR table
# ----------------------------- #

prepare_seurat_metadata <- function(obj, cfg) {
  obj@meta.data |>
    rownames_to_column("barcode")
}

join_metadata_to_tcr <- function(tcr_df, obj, cfg) {
  meta <- prepare_seurat_metadata(obj, cfg)
  
  meta_to_add <- meta |>
    select(
      barcode,
      any_of(c(
        cfg$sample_col,
        "age",
        "sex",
        "diagnosis",
        "years_since_diagnosis",
        cfg$cohort_col,
        cfg$tetramer_col,
        cfg$tetramer_clr_col,
        cfg$tetramer_raw_col,
        cfg$batch_col
      ))
    ) |>
    distinct(barcode, .keep_all = TRUE)
  
  tcr_df |>
    left_join(meta_to_add, by = "barcode") |>
    mutate(
      tetramer = if_else(is.na(.data[[cfg$tetramer_col]]), cfg$negative_label, .data[[cfg$tetramer_col]]),
      tetramer_CLR = if_else(is.na(.data[[cfg$tetramer_clr_col]]), 0, as.numeric(.data[[cfg$tetramer_clr_col]]))
    )
}


# ----------------------------- #
# 5. Within-sample tetramer QC
# ----------------------------- #

summarise_within_sample_conflicts <- function(df, cfg) {
  df2 <- df |>
    filter(
      !is.na(TRB_cdr3), TRB_cdr3 != "",
      !is.na(TRB_v_gene), TRB_v_gene != "",
      !is.na(TRB_j_gene), TRB_j_gene != "",
      !is.na(tetramer), tetramer != ""
    ) |>
    make_clone_key(group_col = cfg$sample_col)
  
  conflict_keys <- df2 |>
    group_by(clone_key) |>
    summarise(n_tetramers = n_distinct(tetramer), .groups = "drop") |>
    filter(n_tetramers > 1) |>
    pull(clone_key)
  
  df_conflict <- df2 |>
    filter(clone_key %in% conflict_keys)
  
  tet_counts <- df_conflict |>
    count(clone_key, tetramer, name = "n_cells") |>
    group_by(clone_key) |>
    summarise(
      tetramer_counts = paste0(tetramer, "=", n_cells, collapse = "; "),
      .groups = "drop"
    )
  
  id_counts <- df_conflict |>
    count(clone_key, id, name = "n_cells") |>
    group_by(clone_key) |>
    summarise(
      id_counts = paste0(id, "=", n_cells, collapse = "; "),
      .groups = "drop"
    )
  
  tet_by_id <- df_conflict |>
    count(clone_key, id, tetramer, name = "n_cells") |>
    arrange(clone_key, id, desc(n_cells)) |>
    group_by(clone_key, id) |>
    summarise(
      tet_counts_in_id = paste0(tetramer, "=", n_cells, collapse = ", "),
      .groups = "drop"
    ) |>
    group_by(clone_key) |>
    summarise(
      tetramer_by_id = paste0(id, ": ", tet_counts_in_id, collapse = " | "),
      .groups = "drop"
    )
  
  conflict_summary <- df_conflict |>
    group_by(clone_key) |>
    summarise(
      sample = first(sample),
      TRB_v_gene = first(TRB_v_gene),
      TRB_j_gene = first(TRB_j_gene),
      TRB_cdr3 = first(TRB_cdr3),
      n_cells_total = n(),
      n_tetramers = n_distinct(tetramer),
      .groups = "drop"
    ) |>
    left_join(tet_counts, by = "clone_key") |>
    left_join(id_counts, by = "clone_key") |>
    left_join(tet_by_id, by = "clone_key") |>
    arrange(desc(n_cells_total), desc(n_tetramers))
  
  list(
    df_with_keys = df2,
    df_conflict = df_conflict,
    conflict_summary = conflict_summary
  )
}

make_within_sample_auto_calls <- function(df_conflict, tet_col = "tetramer") {
  if (nrow(df_conflict) == 0) {
    return(tibble(
      clone_key = character(),
      tetramer_counts = character(),
      tetramer_auto = character(),
      needs_vdjdb = logical(),
      rule_used = character()
    ))
  }
  
  df_conflict |>
    count(clone_key, tet = .data[[tet_col]], name = "n") |>
    group_by(clone_key) |>
    summarise(
      tetramer_counts = paste0(tet, "=", n, collapse = "; "),
      n_neg = sum(n[tet == "negative"], na.rm = TRUE),
      top_nonneg_n = ifelse(
        any(tet != "negative"),
        max(n[tet != "negative"], na.rm = TRUE),
        NA_integer_
      ),
      top_nonneg_tets = ifelse(
        any(tet != "negative"),
        paste(
          tet[tet != "negative" & n == max(n[tet != "negative"], na.rm = TRUE)],
          collapse = "|"
        ),
        NA_character_
      ),
      n_top_nonneg_tets = ifelse(
        any(tet != "negative"),
        sum(tet != "negative" & n == max(n[tet != "negative"], na.rm = TRUE)),
        NA_integer_
      ),
      second_nonneg_n = {
        nn <- n[tet != "negative"]
        if (length(nn) >= 2) sort(nn, decreasing = TRUE)[2] else 0
      },
      .groups = "drop"
    ) |>
    mutate(
      neg_ties_top_nonneg = !is.na(top_nonneg_n) &
        n_neg == top_nonneg_n &
        top_nonneg_n > 0,
      tetramer_auto = case_when(
        !is.na(top_nonneg_n) & n_neg > top_nonneg_n ~ "negative",
        is.na(top_nonneg_n) & n_neg > 0 ~ "negative",
        neg_ties_top_nonneg ~ NA_character_,
        !is.na(top_nonneg_n) &
          n_top_nonneg_tets == 1 &
          top_nonneg_n > n_neg &
          top_nonneg_n > second_nonneg_n ~ str_split_i(top_nonneg_tets, "\\|", 1),
        !is.na(top_nonneg_n) & n_top_nonneg_tets > 1 ~ NA_character_,
        TRUE ~ "negative"
      ),
      needs_vdjdb = neg_ties_top_nonneg |
        (!is.na(top_nonneg_n) & n_top_nonneg_tets > 1),
      rule_used = case_when(
        !is.na(top_nonneg_n) & n_neg > top_nonneg_n ~ "negatives_dominate",
        is.na(top_nonneg_n) & n_neg > 0 ~ "negatives_only",
        neg_ties_top_nonneg ~ "neg_vs_nonneg_tie_needs_vdjdb",
        !is.na(top_nonneg_n) &
          n_top_nonneg_tets == 1 &
          top_nonneg_n > n_neg &
          top_nonneg_n > second_nonneg_n ~ "majority_nonneg",
        !is.na(top_nonneg_n) & n_top_nonneg_tets > 1 ~ "tie_needs_vdjdb",
        TRUE ~ "fallback_negative"
      )
    ) |>
    select(clone_key, tetramer_counts, tetramer_auto, needs_vdjdb, rule_used)
}

read_within_manual_map <- function(cfg) {
  if (!isTRUE(cfg$use_manual_within_sample_calls)) {
    return(tibble(clone_key = character(), tetramer_manual = character()))
  }
  
  conf <- read_xlsx_optional(
    cfg$manual_within_sample_xlsx,
    required = cfg$manual_files_required
  )
  
  if (is.null(conf)) {
    return(tibble(clone_key = character(), tetramer_manual = character()))
  }
  
  if (!"clone_key" %in% colnames(conf)) {
    needed <- c("sample", "TRB_v_gene", "TRB_j_gene", "TRB_cdr3")
    if (!all(needed %in% colnames(conf))) {
      stop(
        "Within-sample manual file needs either clone_key or columns: ",
        paste(needed, collapse = ", ")
      )
    }
    
    conf <- make_clone_key(conf, group_col = "sample")
  }
  
  if (!cfg$within_sample_manual_final_col %in% colnames(conf)) {
    stop(
      "Within-sample manual file missing final call column: ",
      cfg$within_sample_manual_final_col
    )
  }
  
  conf |>
    filter(
      !is.na(.data[[cfg$within_sample_manual_final_col]]),
      .data[[cfg$within_sample_manual_final_col]] != ""
    ) |>
    transmute(
      clone_key = as.character(clone_key),
      tetramer_manual = as.character(.data[[cfg$within_sample_manual_final_col]])
    ) |>
    distinct(clone_key, .keep_all = TRUE)
}

apply_within_sample_qc <- function(df, cfg) {
  conflict_result <- summarise_within_sample_conflicts(df, cfg)
  
  write.csv(
    conflict_result$conflict_summary,
    cfg$within_sample_conflicts_csv,
    row.names = FALSE
  )
  
  auto_calls <- make_within_sample_auto_calls(conflict_result$df_conflict)
  
  conflict_auto_summary <- conflict_result$conflict_summary |>
    left_join(auto_calls, by = "clone_key")
  
  write.csv(
    conflict_auto_summary,
    cfg$within_sample_auto_calls_csv,
    row.names = FALSE
  )
  
  manual_map <- read_within_manual_map(cfg)
  
  df_qc1 <- df |>
    make_clone_key(group_col = cfg$sample_col) |>
    left_join(auto_calls |> select(clone_key, tetramer_auto), by = "clone_key") |>
    left_join(manual_map, by = "clone_key") |>
    mutate(
      tetramer_qc1_final = case_when(
        !is.na(tetramer_manual) ~ tetramer_manual,
        clone_key %in% auto_calls$clone_key & !is.na(tetramer_auto) ~ tetramer_auto,
        TRUE ~ tetramer
      ),
      tetramer_qc1_final = if_else(
        is.na(tetramer_qc1_final),
        cfg$negative_label,
        tetramer_qc1_final
      ),
      qc1_was_overridden = tetramer_qc1_final != tetramer
    )
  
  write.csv(df_qc1, cfg$tcr_qc1_csv, row.names = FALSE)
  
  list(
    df = df_qc1,
    auto_calls = auto_calls,
    manual_map = manual_map
  )
}


# ----------------------------- #
# 6. Global tetramer QC
# ----------------------------- #

summarise_global_conflicts <- function(df, tet_col, cfg) {
  df_qc2_base <- df |>
    filter(
      !is.na(TRB_cdr3), TRB_cdr3 != "",
      !is.na(TRB_v_gene), TRB_v_gene != "",
      !is.na(TRB_j_gene), TRB_j_gene != "",
      !is.na(.data[[tet_col]]), .data[[tet_col]] != ""
    ) |>
    mutate(tetramer = .data[[tet_col]]) |>
    make_clone_key_global()
  
  global_conflict_keys <- df_qc2_base |>
    group_by(clone_key_global) |>
    summarise(n_tetramers = n_distinct(tetramer), .groups = "drop") |>
    filter(n_tetramers > 1) |>
    pull(clone_key_global)
  
  df_global_conflict <- df_qc2_base |>
    filter(clone_key_global %in% global_conflict_keys)
  
  tet_counts_global <- df_global_conflict |>
    count(clone_key_global, tetramer, name = "n_cells") |>
    group_by(clone_key_global) |>
    summarise(
      tetramer_counts = paste0(tetramer, "=", n_cells, collapse = "; "),
      .groups = "drop"
    )
  
  sample_counts_global <- df_global_conflict |>
    count(clone_key_global, sample, name = "n_cells") |>
    group_by(clone_key_global) |>
    summarise(
      sample_counts = paste0(sample, "=", n_cells, collapse = "; "),
      .groups = "drop"
    )
  
  tet_by_sample <- df_global_conflict |>
    count(clone_key_global, sample, tetramer, name = "n_cells") |>
    arrange(clone_key_global, sample, desc(n_cells)) |>
    group_by(clone_key_global, sample) |>
    summarise(
      tet_counts_in_sample = paste0(tetramer, "=", n_cells, collapse = ", "),
      .groups = "drop"
    ) |>
    group_by(clone_key_global) |>
    summarise(
      tetramer_by_sample = paste0(sample, ": ", tet_counts_in_sample, collapse = " | "),
      .groups = "drop"
    )
  
  global_summary <- df_global_conflict |>
    group_by(clone_key_global) |>
    summarise(
      TRB_v_gene = first(TRB_v_gene),
      TRB_j_gene = first(TRB_j_gene),
      TRB_cdr3 = first(TRB_cdr3),
      n_cells_total = n(),
      n_tetramers = n_distinct(tetramer),
      n_samples = n_distinct(sample),
      .groups = "drop"
    ) |>
    left_join(tet_counts_global, by = "clone_key_global") |>
    left_join(sample_counts_global, by = "clone_key_global") |>
    left_join(tet_by_sample, by = "clone_key_global") |>
    arrange(desc(n_cells_total), desc(n_samples), desc(n_tetramers))
  
  list(
    df_base = df_qc2_base,
    df_conflict = df_global_conflict,
    global_summary = global_summary
  )
}

make_global_auto_calls <- function(df_global_conflict) {
  if (nrow(df_global_conflict) == 0) {
    return(tibble(
      clone_key_global = character(),
      tetramer_counts = character(),
      tetramer_qc2_auto = character(),
      qc2_needs_vdjdb = logical(),
      qc2_rule = character()
    ))
  }
  
  df_global_conflict |>
    count(clone_key_global, tetramer, name = "n_cells") |>
    group_by(clone_key_global) |>
    summarise(
      tetramer_counts = paste0(tetramer, "=", n_cells, collapse = "; "),
      top_n = max(n_cells),
      n_top = sum(n_cells == max(n_cells)),
      top_tets = paste(tetramer[n_cells == max(n_cells)], collapse = "|"),
      tetramer_qc2_auto = ifelse(
        n_top == 1,
        str_split_i(top_tets, "\\|", 1),
        NA_character_
      ),
      qc2_needs_vdjdb = n_top > 1,
      qc2_rule = ifelse(n_top == 1, "qc2_highest_unique", "qc2_tie_needs_vdjdb"),
      .groups = "drop"
    ) |>
    select(
      clone_key_global,
      tetramer_counts,
      tetramer_qc2_auto,
      qc2_needs_vdjdb,
      qc2_rule
    )
}

read_global_manual_map <- function(cfg) {
  if (!isTRUE(cfg$use_manual_global_calls)) {
    return(tibble(clone_key_global = character(), tetramer_global_manual = character()))
  }
  
  global_map_raw <- read_xlsx_optional(
    cfg$manual_global_xlsx,
    required = cfg$manual_files_required
  )
  
  if (is.null(global_map_raw)) {
    return(tibble(clone_key_global = character(), tetramer_global_manual = character()))
  }
  
  if (!"clone_key_global" %in% colnames(global_map_raw)) {
    needed <- c("TRB_v_gene", "TRB_j_gene", "TRB_cdr3")
    if (!all(needed %in% colnames(global_map_raw))) {
      stop(
        "Global manual file needs either clone_key_global or columns: ",
        paste(needed, collapse = ", ")
      )
    }
    
    global_map_raw <- global_map_raw |>
      mutate(clone_key_global = paste(TRB_v_gene, TRB_j_gene, TRB_cdr3, sep = "|"))
  }
  
  if (!cfg$global_manual_final_col %in% colnames(global_map_raw)) {
    stop(
      "Global manual file missing final call column: ",
      cfg$global_manual_final_col
    )
  }
  
  global_map_raw |>
    filter(
      !is.na(.data[[cfg$global_manual_final_col]]),
      .data[[cfg$global_manual_final_col]] != ""
    ) |>
    transmute(
      clone_key_global = as.character(clone_key_global),
      tetramer_global_manual = as.character(.data[[cfg$global_manual_final_col]])
    ) |>
    distinct(clone_key_global, .keep_all = TRUE)
}

apply_global_qc <- function(df, cfg) {
  df <- df |>
    mutate(
      tetramer = tetramer_qc1_final,
      tetramer = if_else(is.na(tetramer), cfg$negative_label, tetramer)
    )
  
  global_result <- summarise_global_conflicts(
    df = df,
    tet_col = "tetramer",
    cfg = cfg
  )
  
  write.csv(
    global_result$global_summary,
    cfg$global_conflicts_csv,
    row.names = FALSE
  )
  
  global_auto <- make_global_auto_calls(global_result$df_conflict)
  
  write.csv(
    global_auto,
    cfg$global_auto_calls_csv,
    row.names = FALSE
  )
  
  global_manual <- read_global_manual_map(cfg)
  
  df_qc2 <- df |>
    make_clone_key_global() |>
    left_join(global_auto |> select(clone_key_global, tetramer_qc2_auto), by = "clone_key_global") |>
    left_join(global_manual, by = "clone_key_global") |>
    mutate(
      tetramer_qc2_final = case_when(
        !is.na(tetramer_global_manual) ~ tetramer_global_manual,
        clone_key_global %in% global_auto$clone_key_global &
          !is.na(tetramer_qc2_auto) ~ tetramer_qc2_auto,
        TRUE ~ tetramer
      ),
      tetramer_qc2_final = if_else(
        is.na(tetramer_qc2_final),
        cfg$negative_label,
        tetramer_qc2_final
      ),
      qc2_was_overridden = tetramer_qc2_final != tetramer
    )
  
  write.csv(
    df_qc2,
    cfg$tcr_qc2_csv,
    row.names = FALSE
  )
  
  list(
    df = df_qc2,
    global_auto = global_auto,
    global_manual = global_manual
  )
}


# ----------------------------- #
# 7. HLA mismatch filtering
# ----------------------------- #

apply_hla_mismatch_filter <- function(df, cfg, tet_col = "tetramer_qc2_final") {
  if (!isTRUE(cfg$apply_hla_filter)) {
    return(list(df = df, removed = tibble()))
  }
  
  tetramer_hla <- read_xlsx_optional(cfg$tetramer_hla_xlsx, required = TRUE)
  hla <- read_xlsx_optional(cfg$hla_information_xlsx, required = TRUE)
  
  tetramer_long <- tetramer_hla |>
    pivot_longer(
      cols = -tetramer,
      names_to = "allele",
      values_to = "tetramer_hla_status"
    ) |>
    filter(tetramer_hla_status == "positive") |>
    transmute(
      tetramer = as.character(tetramer),
      allele = as.character(allele)
    ) |>
    distinct()
  
  hla_long <- hla |>
    pivot_longer(
      cols = -sample,
      names_to = "allele",
      values_to = "sample_hla_status"
    ) |>
    filter(sample_hla_status == "positive") |>
    transmute(
      sample = as.character(sample),
      allele = as.character(allele)
    ) |>
    distinct()
  
  valid_pairs <- tetramer_long |>
    inner_join(hla_long, by = "allele", relationship = "many-to-many") |>
    distinct(tetramer, sample)
  
  known_tetramers <- unique(valid_pairs$tetramer)
  known_samples <- unique(valid_pairs$sample)
  
  df2 <- df |>
    mutate(
      tetramer_raw_final = .data[[tet_col]],
      tetramer_key = str_replace(str_trim(.data[[tet_col]]), "\\*$", "")
    )
  
  bad_calls <- df2 |>
    filter(!is.na(tetramer_key), tetramer_key != cfg$negative_label) |>
    filter(tetramer_key %in% known_tetramers, sample %in% known_samples) |>
    anti_join(
      valid_pairs,
      by = c("tetramer_key" = "tetramer", "sample" = "sample")
    )
  
  df_filtered <- df2 |>
    anti_join(
      bad_calls |> select(barcode) |> distinct(),
      by = "barcode"
    )
  
  write.csv(
    bad_calls,
    cfg$hla_removed_cells_csv,
    row.names = FALSE
  )
  
  list(
    df = df_filtered,
    removed = bad_calls
  )
}


# ----------------------------- #
# 8. Project GEM calls onto activated non-GEM cells
# ----------------------------- #

project_gem_tetramers_to_non_gem <- function(df_gem_filtered, df_non_gem_raw, obj, cfg) {
  meta <- prepare_seurat_metadata(obj, cfg)
  
  meta_to_add <- meta |>
    select(
      barcode,
      any_of(c(
        cfg$sample_col,
        "age",
        "sex",
        "diagnosis",
        "years_since_diagnosis",
        cfg$cohort_col,
        cfg$batch_col
      ))
    ) |>
    distinct(barcode, .keep_all = TRUE)
  
  df_non <- df_non_gem_raw |>
    left_join(meta_to_add, by = "barcode") |>
    filter(
      !is.na(sample), sample != "",
      !is.na(TRB_v_gene), TRB_v_gene != "",
      !is.na(TRB_j_gene), TRB_j_gene != "",
      !is.na(TRB_cdr3), TRB_cdr3 != ""
    ) |>
    mutate(
      clone_key = paste(sample, TRB_v_gene, TRB_j_gene, TRB_cdr3, sep = "|")
    )
  
  gem_tet_col <- "tetramer_qc2_final"
  
  df_gem <- df_gem_filtered |>
    filter(
      !is.na(sample), sample != "",
      !is.na(TRB_v_gene), TRB_v_gene != "",
      !is.na(TRB_j_gene), TRB_j_gene != "",
      !is.na(TRB_cdr3), TRB_cdr3 != "",
      !is.na(.data[[gem_tet_col]]), .data[[gem_tet_col]] != ""
    ) |>
    mutate(
      clone_key = paste(sample, TRB_v_gene, TRB_j_gene, TRB_cdr3, sep = "|")
    )
  
  gem_clone_map_check <- df_gem |>
    group_by(sample, clone_key) |>
    summarise(n_tet = n_distinct(.data[[gem_tet_col]]), .groups = "drop") |>
    filter(n_tet > 1)
  
  if (nrow(gem_clone_map_check) > 0) {
    print(head(gem_clone_map_check, 50))
    stop("Baseline GEM still has >1 tetramer_qc2_final for some sample-level clones.")
  }
  
  gem_clone_map <- df_gem |>
    distinct(sample, clone_key, .data[[gem_tet_col]]) |>
    rename(tetramer_from_GEM = !!sym(gem_tet_col))
  
  df_non_projected <- df_non |>
    left_join(gem_clone_map, by = c("sample", "clone_key")) |>
    mutate(
      tetramer = case_when(
        !is.na(tetramer_from_GEM) ~ tetramer_from_GEM,
        str_detect(batch, "EBVLCL1DEX|EBVLCL2DEX") ~ cfg$activated_enriched_label,
        str_detect(batch, "EBVLCL1CD8|EBVLCL2CD8") ~ cfg$activated_lcl_responsive_label,
        TRUE ~ cfg$negative_label
      ),
      tetramer_source = case_when(
        !is.na(tetramer_from_GEM) ~ "projected_from_GEM_same_sample",
        str_detect(batch, "EBVLCL1DEX|EBVLCL2DEX") ~ "batch_inferred",
        str_detect(batch, "EBVLCL1CD8|EBVLCL2CD8") ~ "batch_inferred",
        TRUE ~ "negative"
      ),
      data_source = "nonGEM_activated"
    ) |>
    select(-tetramer_from_GEM)
  
  df_gem_final <- df_gem_filtered |>
    mutate(
      tetramer = tetramer_qc2_final,
      tetramer = if_else(is.na(tetramer), cfg$negative_label, tetramer),
      data_source = "GEM_baseline",
      tetramer_source = "GEM_measured_QC1_QC2_HLA"
    )
  
  df_all <- bind_rows(df_gem_final, df_non_projected)
  
  write.csv(
    df_all,
    cfg$projected_all_tcr_csv,
    row.names = FALSE
  )
  
  df_all
}


# ----------------------------- #
# 9. Add VDJ/TCR metadata back to Seurat
# ----------------------------- #

add_tcr_metadata_to_seurat <- function(obj, df_all, cfg) {
  vdj_to_add <- df_all |>
    transmute(
      barcode = as.character(barcode),
      tetramer = as.character(tetramer),
      data_source = as.character(data_source),
      tetramer_source = as.character(tetramer_source),
      TRA_v_gene,
      TRA_d_gene,
      TRA_j_gene,
      TRA_c_gene,
      TRA_cdr3,
      TRA_cdr3_nt,
      TRB_v_gene,
      TRB_d_gene,
      TRB_j_gene,
      TRB_c_gene,
      TRB_cdr3,
      TRB_cdr3_nt
    ) |>
    distinct(barcode, .keep_all = TRUE)
  
  meta <- obj@meta.data |>
    rownames_to_column("seurat_cell") |>
    mutate(barcode = seurat_cell)
  
  meta_matched <- meta |>
    inner_join(vdj_to_add |> select(barcode), by = "barcode")
  
  cells_keep <- meta_matched$seurat_cell
  
  obj_vdj <- subset(obj, cells = cells_keep)
  
  vdj_meta_for_seurat <- meta_matched |>
    select(seurat_cell, barcode) |>
    left_join(vdj_to_add, by = "barcode") |>
    select(-barcode) |>
    distinct(seurat_cell, .keep_all = TRUE)
  
  rownames(vdj_meta_for_seurat) <- vdj_meta_for_seurat$seurat_cell
  vdj_meta_for_seurat$seurat_cell <- NULL
  
  cols_to_drop <- intersect(colnames(obj_vdj@meta.data), colnames(vdj_meta_for_seurat))
  
  if (length(cols_to_drop) > 0) {
    obj_vdj@meta.data <- obj_vdj@meta.data |>
      select(-all_of(cols_to_drop))
  }
  
  AddMetaData(obj_vdj, metadata = vdj_meta_for_seurat)
}


# ----------------------------- #
# 10. Add virus/tetramer annotation
# ----------------------------- #

add_multimer_information <- function(obj, cfg) {
  multimer_info <- read_xlsx_optional(cfg$multimer_information_xlsx, required = TRUE)
  
  multimer_info <- multimer_info |>
    add_row(
      tetramer = cfg$activated_enriched_label,
      virus = cfg$activated_enriched_label,
      lifecycle = cfg$activated_enriched_label,
      latency = cfg$activated_enriched_label,
      antigen = cfg$activated_enriched_label
    ) |>
    add_row(
      tetramer = cfg$activated_lcl_responsive_label,
      virus = cfg$activated_lcl_responsive_label,
      lifecycle = cfg$activated_lcl_responsive_label,
      latency = cfg$activated_lcl_responsive_label,
      antigen = cfg$activated_lcl_responsive_label
    ) |>
    add_row(
      tetramer = cfg$negative_label,
      virus = cfg$negative_label,
      lifecycle = cfg$negative_label,
      latency = cfg$negative_label,
      antigen = cfg$negative_label
    )
  
  lookup <- multimer_info |>
    distinct(tetramer, .keep_all = TRUE) |>
    select(tetramer, virus, lifecycle, latency, antigen)
  
  meta <- obj@meta.data |>
    rownames_to_column("cell_barcode")
  
  # Drop old annotation columns before clean join.
  meta <- meta |>
    select(-any_of(c("virus", "lifecycle", "latency", "antigen")))
  
  meta2 <- meta |>
    left_join(lookup, by = "tetramer")
  
  # CD19.1 is treated as negative before batch inference.
  is_cd19 <- meta2$tetramer == "CD19.1"
  
  meta2 <- meta2 |>
    mutate(
      tetramer = if_else(is_cd19, cfg$negative_label, tetramer),
      virus = if_else(is_cd19, cfg$negative_label, virus),
      latency = if_else(is_cd19, cfg$negative_label, latency),
      antigen = if_else(is_cd19, cfg$negative_label, antigen),
      lifecycle = if_else(is_cd19, cfg$negative_label, lifecycle)
    ) |>
    mutate(
      tetramer = case_when(
        tetramer == cfg$negative_label &
          str_detect(batch, cfg$ebv_enriched_batch_pattern) ~ cfg$activated_enriched_label,
        tetramer == cfg$negative_label &
          str_detect(batch, cfg$lcl_responsive_batch_pattern) ~ cfg$activated_lcl_responsive_label,
        TRUE ~ tetramer
      ),
      virus = case_when(
        (is.na(virus) | virus == cfg$negative_label) &
          str_detect(batch, cfg$ebv_enriched_batch_pattern) ~ cfg$activated_enriched_label,
        (is.na(virus) | virus == cfg$negative_label) &
          str_detect(batch, cfg$lcl_responsive_batch_pattern) ~ cfg$activated_lcl_responsive_label,
        TRUE ~ virus
      ),
      lifecycle = case_when(
        (is.na(lifecycle) | lifecycle == cfg$negative_label) &
          str_detect(batch, cfg$ebv_enriched_batch_pattern) ~ cfg$activated_enriched_label,
        (is.na(lifecycle) | lifecycle == cfg$negative_label) &
          str_detect(batch, cfg$lcl_responsive_batch_pattern) ~ cfg$activated_lcl_responsive_label,
        TRUE ~ lifecycle
      ),
      latency = case_when(
        (is.na(latency) | latency == cfg$negative_label) &
          str_detect(batch, cfg$ebv_enriched_batch_pattern) ~ cfg$activated_enriched_label,
        (is.na(latency) | latency == cfg$negative_label) &
          str_detect(batch, cfg$lcl_responsive_batch_pattern) ~ cfg$activated_lcl_responsive_label,
        TRUE ~ latency
      ),
      antigen = case_when(
        (is.na(antigen) | antigen == cfg$negative_label) &
          str_detect(batch, cfg$ebv_enriched_batch_pattern) ~ cfg$activated_enriched_label,
        (is.na(antigen) | antigen == cfg$negative_label) &
          str_detect(batch, cfg$lcl_responsive_batch_pattern) ~ cfg$activated_lcl_responsive_label,
        TRUE ~ antigen
      )
    )
  
  meta2 <- meta2 |>
    column_to_rownames("cell_barcode")
  
  obj@meta.data <- meta2
  
  obj
}


# ----------------------------- #
# 11. ADT CD8 gating and publication subsets
# ----------------------------- #

add_adt_marker_metadata <- function(obj, cfg) {
  required_features <- c(cfg$adt_cd4_feature, cfg$adt_cd8_feature, cfg$adt_cd3_feature)
  
  if (!cfg$adt_assay %in% Assays(obj)) {
    warning("ADT assay not found. Skipping ADT marker metadata.")
    return(obj)
  }
  
  adt_features <- rownames(obj[[cfg$adt_assay]])
  missing_features <- setdiff(required_features, adt_features)
  
  if (length(missing_features) > 0) {
    warning(
      "Missing ADT feature(s): ",
      paste(missing_features, collapse = ", "),
      ". Skipping ADT marker metadata."
    )
    return(obj)
  }
  
  adt_mat <- get_assay_data_compat(
    seurat_obj = obj,
    assay = cfg$adt_assay,
    layer_or_slot = "data"
  )[required_features, , drop = FALSE]
  
  adt_df <- t(as.matrix(adt_mat))
  colnames(adt_df) <- c("ADT_CD4.1_norm", "ADT_CD8_norm", "ADT_CD3_norm")
  
  AddMetaData(obj, adt_df)
}

apply_final_cd8_gating <- function(obj, cfg) {
  if (!isTRUE(cfg$perform_adt_cd8_gating)) {
    return(obj)
  }
  
  obj <- add_adt_marker_metadata(obj, cfg)
  
  required_cols <- c("ADT_CD3_norm", "ADT_CD4.1_norm", "ADT_CD8_norm")
  
  if (!all(required_cols %in% colnames(obj@meta.data))) {
    warning("ADT gating columns missing. Returning ungated object.")
    return(obj)
  }
  
  meta <- obj@meta.data
  
  is_gemebv_ebvlcl <- grepl("^GEMEBV", meta[[cfg$batch_col]]) |
    grepl("^EBVLCL", meta[[cfg$batch_col]])
  
  is_other_gated <- !is_gemebv_ebvlcl &
    meta$ADT_CD3_norm > cfg$adt_cd3_min &
    meta$ADT_CD4.1_norm < cfg$adt_cd4_max &
    meta$ADT_CD8_norm > cfg$adt_cd8_min
  
  keep_cells <- rownames(meta)[is_gemebv_ebvlcl | is_other_gated]
  
  subset(obj, cells = keep_cells)
}

remove_cd4_activated_batches <- function(obj, cfg) {
  subset(
    obj,
    subset = !(batch %in% cfg$cd4_batch_exclude)
  )
}

make_publication_subsets <- function(obj, cfg) {
  merged_cd8 <- apply_final_cd8_gating(obj, cfg)
  merged_cd8 <- remove_cd4_activated_batches(merged_cd8, cfg)
  
  merged_gem_only <- subset(
    merged_cd8,
    subset = grepl("^GEM", batch)
  )
  
  meta <- merged_cd8@meta.data
  
  ebvlcl_samples <- unique(meta$sample[grepl("^EBVLCL", meta$batch)])
  gem_like_samples <- unique(meta$sample[grepl("^GEM", meta$batch)])
  shared_samples <- intersect(ebvlcl_samples, gem_like_samples)
  
  merged_counterpart <- subset(
    merged_cd8,
    subset = sample %in% shared_samples
  )
  
  list(
    combined = merged_cd8,
    baseline = merged_gem_only,
    activated_counterpart = merged_counterpart,
    shared_samples = shared_samples
  )
}

clean_final_metadata <- function(obj, cfg) {
  obj@meta.data <- obj@meta.data |>
    select(-any_of(cfg$metadata_columns_to_remove))
  
  obj
}


# ----------------------------- #
# 12. Run pipeline
# ----------------------------- #

setup_output_dirs(config)

message("Loading inputs...")
merged_EBV_noCD4 <- load_seurat_object(config$input_seurat_rds)
vdj_df <- read_csv_required(config$input_vdj_csv)

message("Preparing collapsed TCR table...")
tcr_prepared <- prepare_clean_tcr_table(vdj_df, config)
df_final <- tcr_prepared$df
chain_summary <- tcr_prepared$chain_summary

message("Splitting baseline GEM and activated non-GEM VDJ rows...")
df_gem_raw <- df_final |>
  filter(grepl(config$baseline_batch_pattern, batch))

df_non_gem_raw <- df_final |>
  filter(!grepl(config$baseline_batch_pattern, batch))

message("Joining Seurat metadata to baseline GEM TCR rows...")
df_gem_with_meta <- join_metadata_to_tcr(
  tcr_df = df_gem_raw,
  obj = merged_EBV_noCD4,
  cfg = config
) |>
  filter(!is.na(sample))

df_gem_with_meta <- df_gem_with_meta |>
  select(
    barcode,
    batch,
    id,
    tetramer,
    tetramer_CLR,
    TRA_v_gene,
    TRA_d_gene,
    TRA_j_gene,
    TRA_c_gene,
    TRA_cdr3,
    TRA_cdr3_nt,
    TRB_v_gene,
    TRB_d_gene,
    TRB_j_gene,
    TRB_c_gene,
    TRB_cdr3,
    TRB_cdr3_nt,
    sample,
    cohort,
    diagnosis,
    age,
    sex,
    everything()
  )

write.csv(
  df_gem_with_meta,
  config$tcr_gem_only_csv,
  row.names = FALSE
)

message("Running within-sample tetramer QC...")
qc1_result <- apply_within_sample_qc(df_gem_with_meta, config)
df_qc1 <- qc1_result$df

message("Running global tetramer QC...")
qc2_result <- apply_global_qc(df_qc1, config)
df_qc2 <- qc2_result$df

message("Applying HLA mismatch filter...")
hla_result <- apply_hla_mismatch_filter(
  df = df_qc2,
  cfg = config,
  tet_col = "tetramer_qc2_final"
)

df_gem_filtered <- hla_result$df

message("Projecting baseline GEM tetramer calls onto activated non-GEM cells...")
df_all <- project_gem_tetramers_to_non_gem(
  df_gem_filtered = df_gem_filtered,
  df_non_gem_raw = df_non_gem_raw,
  obj = merged_EBV_noCD4,
  cfg = config
)

message("Adding curated TCR metadata back to Seurat object...")
merged_EBV_noCD4_vdj <- add_tcr_metadata_to_seurat(
  obj = merged_EBV_noCD4,
  df_all = df_all,
  cfg = config
)

saveRDS(
  merged_EBV_noCD4_vdj,
  config$output_vdj_seurat_rds
)

message("Adding virus/lifecycle/latency/antigen annotations...")
annotated_obj <- add_multimer_information(
  obj = merged_EBV_noCD4_vdj,
  cfg = config
)

annotated_obj <- clean_final_metadata(annotated_obj, config)

saveRDS(
  annotated_obj,
  config$output_annotated_seurat_rds
)

message("Creating publication subsets...")
publication_objects <- make_publication_subsets(
  obj = annotated_obj,
  cfg = config
)

saveRDS(
  publication_objects$combined,
  config$output_combined_publication_rds
)

saveRDS(
  publication_objects$baseline,
  config$output_baseline_publication_rds
)

saveRDS(
  publication_objects$activated_counterpart,
  config$output_activated_publication_rds
)

message("Saving QC summary...")
qc_summary <- tibble(
  metric = c(
    "input_vdj_rows",
    "collapsed_tcr_cells",
    "baseline_gem_tcr_cells",
    "within_sample_manual_clone_calls",
    "global_manual_clone_calls",
    "hla_mismatched_cells_removed",
    "combined_publication_cells",
    "baseline_publication_cells",
    "activated_counterpart_publication_cells",
    "activated_counterpart_shared_samples"
  ),
  value = c(
    nrow(vdj_df),
    nrow(df_final),
    nrow(df_gem_with_meta),
    nrow(qc1_result$manual_map),
    nrow(qc2_result$global_manual),
    nrow(hla_result$removed),
    ncol(publication_objects$combined),
    ncol(publication_objects$baseline),
    ncol(publication_objects$activated_counterpart),
    length(publication_objects$shared_samples)
  )
)

qc_summary_full <- bind_rows(
  qc_summary,
  chain_summary |>
    mutate(
      metric = paste0(chain, "_", c(
        "n_cells_with_two_chains",
        "n_cells_dropped_equal_or_missing_umi",
        "n_cells_resolved_by_umi"
      )[row_number()]
      ) |>
        as.character(),
      value = c(
        n_cells_with_two_chains,
        n_cells_dropped_equal_or_missing_umi,
        n_cells_resolved_by_umi
      )
    ) |>
    select(metric, value)
)

write.csv(
  qc_summary_full,
  config$qc_summary_csv,
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  config$session_info_file
)

message("\nTCR QC and publication object export complete.")
message("Saved VDJ-annotated object to: ", config$output_vdj_seurat_rds)
message("Saved annotated object to: ", config$output_annotated_seurat_rds)
message("Saved combined publication object to: ", config$output_combined_publication_rds)
message("Saved baseline publication object to: ", config$output_baseline_publication_rds)
message("Saved activated-paired publication object to: ", config$output_activated_publication_rds)
message("Saved QC summary to: ", config$qc_summary_csv)
message("Saved session info to: ", config$session_info_file)