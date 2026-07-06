#!/usr/bin/env Rscript

# ==============================================================================
# Public TRB clonotype sharing from baseline EBV Seurat RDS
# ==============================================================================
#
# Purpose:
#   This script analyses public TRB clonotype sharing across EBV tetramer
#   specificities directly from the final baseline EBV Seurat RDS object.
#
#   It replaces older Excel-based workflows by reading:
#     5_baseline_EBV_cluster_reannotated.rds
#
#   The script:
#     1. Loads the baseline EBV Seurat object.
#     2. Extracts metadata from merged@meta.data.
#     3. Filters to GEMEBV batches and target tetramers.
#     4. Builds TRB-only clonotypes:
#          TRB_cdr3 + TRB_v_gene + TRB_j_gene
#     5. Creates pairwise participant-sharing files for each tetramer.
#     6. Defines public clonotypes within each tetramer.
#     7. Calculates participant-level public sharing burden.
#     8. Summarises public and cross-cohort public sharing per tetramer.
#     9. Runs one-vs-rest Fisher and Wilcoxon tests.
#    10. Calculates pairwise participant overlap metrics.
#    11. Runs pairwise tetramer comparisons.
#    12. Generates summary plots and lower-triangle heatmap.
#    13. Exports CSV outputs, plots, an Excel workbook, and session info.
#
# Public clonotype definition:
#   A TRB clonotype is public within a tetramer if:
#     - it is detected in at least two independent participants/samples, and
#     - it has at least config$min_public_clone_cells total cells across those
#       participants/samples.
#
# Notes:
#   - Paths are intentionally generic for GitHub.
#   - Update the config block before running.
#   - The default participant/sample identifier is `id`, because this project has
#     used `id` for RDS-based participant-level analyses. If your object should
#     use `sample` instead, change config$sample_col to "sample".
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
  library(purrr)
  library(ggplot2)
  library(readr)
  library(openxlsx)
  library(scales)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  
  # Input RDS
  input_rds = "path/to/Publication_RDS/Baseline_EBV_integrated_seurat.rds",
  
  # Output directory
  outdir = "path/to/Publication_data/EBV_baseline/Public_TRB_sharing",
  
  # Metadata columns
  sample_col = "id",          # change to "sample" if required
  cohort_col = "cohort",
  batch_col = "batch",
  tetramer_col = "tetramer",
  
  trb_cdr3_col = "TRB_cdr3",
  trb_v_col = "TRB_v_gene",
  trb_j_col = "TRB_j_gene",
  
  # GEMEBV filtering
  gemebv_pattern = "GEMEBV",
  
  # Tetramers to analyse
  target_tetramers = c(
    "GLCT",
    "RAKF",
    "RPPI",
    "FLRG",
    "QAKW",
    "LLDF",
    "YLQQ",
    "CLGG",
    "FLYA",
    "RPQK*",
    "YNLR*"
  ),
  
  # Public clonotype definition
  min_public_clone_cells = 3,
  
  # Pairwise participant-sharing files
  min_pair_cells_shared = 3,
  log_transform_pairwise_percent = TRUE,
  
  # Cohort labels
  control_label = "Control",
  ms_label = "MS",
  
  # Plot settings
  plot_width = 7,
  plot_height = 5,
  plot_dpi = 300,
  heatmap_width = 7.5,
  heatmap_height = 6.5,
  heatmap_dpi = 300,
  
  # Output files
  filtered_metadata_csv = "GEMEBV_TRB_public_sharing_filtered_metadata.csv",
  tetramer_summary_csv = "EBV_tetramer_public_TRB_sharing_summary.csv",
  sample_level_sharing_csv = "EBV_sample_level_public_TRB_sharing.csv",
  public_clone_lookup_csv = "EBV_public_TRB_clone_lookup.csv",
  clone_counts_csv = "EBV_TRB_clone_counts_by_tetramer_sample.csv",
  pairwise_summary_csv = "EBV_pairwise_sample_overlap_summary.csv",
  pairwise_overlap_csv = "EBV_pairwise_sample_overlap_full.csv",
  fisher_tests_csv = "EBV_public_TRB_sharing_fisher_tests.csv",
  wilcox_tests_csv = "EBV_public_TRB_sharing_wilcox_tests.csv",
  pairwise_tetramer_tests_csv = "EBV_pairwise_public_TRB_burden_unique_tests.csv",
  pairwise_heatmap_values_csv = "EBV_pairwise_public_TRB_burden_lower_triangle_heatmap_values.csv",
  
  public_fraction_plot = "EBV_public_cell_fraction_by_tetramer.png",
  public_prevalence_plot = "EBV_public_clone_prevalence_by_tetramer.png",
  cross_cohort_fraction_plot = "EBV_cross_cohort_public_cell_fraction_by_tetramer.png",
  cross_cohort_prevalence_plot = "EBV_cross_cohort_public_clone_prevalence_by_tetramer.png",
  pairwise_heatmap_plot = "EBV_pairwise_public_TRB_burden_lower_triangle_heatmap.png",
  
  workbook_file = "EBV_tetramer_sample_level_public_TRB_clonotype_sharing_stats.xlsx",
  session_info_file = "sessionInfo_public_TRB_clonotype_sharing_from_RDS.txt"
)


# ----------------------------- #
# 2. Helper functions
# ----------------------------- #

create_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
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
  str_trim(as.character(x))
}


clean_vdj <- function(x) {
  
  x <- clean_chr(x)
  x[x %in% c("", "NA", "None", "none", "nan", "NaN", "NULL", "null")] <- NA_character_
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


clean_tetramer <- function(x) {
  
  # Handles both "GLCT" and "A02_GLCT" style names.
  x %>%
    clean_chr() %>%
    str_remove("^[^_]+_")
}


safe_filename <- function(x) {
  
  x %>%
    str_replace_all("\\*", "STAR") %>%
    str_replace_all("[<>:\"/\\\\|\\?\\*]", "_")
}


safe_fold <- function(a, b) {
  
  if (is.na(a) || is.na(b)) {
    return(NA_real_)
  }
  
  if (b == 0 && a == 0) {
    return(NA_real_)
  }
  
  if (b == 0 && a > 0) {
    return(Inf)
  }
  
  a / b
}


adjust_bh_without_na <- function(p) {
  
  out <- rep(NA_real_, length(p))
  idx <- !is.na(p)
  
  if (any(idx)) {
    out[idx] <- p.adjust(p[idx], method = "BH")
  }
  
  out
}


binom_ci <- function(x, n) {
  
  if (is.na(x) || is.na(n) || n == 0) {
    return(tibble(ci_lower = NA_real_, ci_upper = NA_real_))
  }
  
  bt <- binom.test(x, n)
  
  tibble(
    ci_lower = unname(bt$conf.int[1]),
    ci_upper = unname(bt$conf.int[2])
  )
}


morisita_horn <- function(x, y) {
  
  if (sum(x) == 0 || sum(y) == 0) {
    return(NA_real_)
  }
  
  p <- x / sum(x)
  q <- y / sum(y)
  
  numerator <- 2 * sum(p * q)
  denominator <- sum(p^2) + sum(q^2)
  
  if (denominator == 0) {
    return(NA_real_)
  }
  
  numerator / denominator
}


# ----------------------------- #
# 3. Load RDS and build filtered metadata
# ----------------------------- #

create_dir(config$outdir)

message("Loading baseline EBV object:")
message(config$input_rds)

merged <- load_seurat_object(config$input_rds)

meta <- merged@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell_barcode")

required_cols <- c(
  config$sample_col,
  config$cohort_col,
  config$batch_col,
  config$tetramer_col,
  config$trb_cdr3_col,
  config$trb_v_col,
  config$trb_j_col
)

check_required_columns(
  df = meta,
  required_cols = required_cols,
  object_name = "baseline EBV metadata"
)

df <- meta %>%
  mutate(
    sample = clean_chr(.data[[config$sample_col]]),
    cohort = clean_cohort(
      .data[[config$cohort_col]],
      control_label = config$control_label,
      ms_label = config$ms_label
    ),
    batch = clean_chr(.data[[config$batch_col]]),
    tetramer = clean_tetramer(.data[[config$tetramer_col]]),
    
    TRB_cdr3 = clean_vdj(.data[[config$trb_cdr3_col]]),
    TRB_v_gene = clean_vdj(.data[[config$trb_v_col]]),
    TRB_j_gene = clean_vdj(.data[[config$trb_j_col]]),
    
    cell_id_for_sharing = paste(sample, cell_barcode, sep = "_"),
    
    clonotype_for_sharing = if_else(
      !is.na(TRB_cdr3) & !is.na(TRB_v_gene) & !is.na(TRB_j_gene),
      paste(TRB_cdr3, TRB_v_gene, TRB_j_gene, sep = "|"),
      NA_character_
    )
  ) %>%
  filter(
    str_detect(batch, config$gemebv_pattern),
    tetramer %in% config$target_tetramers,
    
    !is.na(sample), sample != "",
    !is.na(cohort), cohort != "",
    !is.na(clonotype_for_sharing), clonotype_for_sharing != ""
  ) %>%
  distinct(cell_id_for_sharing, .keep_all = TRUE) %>%
  mutate(
    tetramer = factor(tetramer, levels = config$target_tetramers)
  )

if (nrow(df) == 0) {
  stop("No cells remained after GEMEBV/tetramer/TRB filtering.")
}

filtered_metadata_path <- file.path(config$outdir, config$filtered_metadata_csv)

write_csv(df, filtered_metadata_path)

message("Cells per tetramer:")
print(df %>% count(tetramer, sort = FALSE))

message("Cells per batch:")
print(df %>% count(batch, sort = TRUE))

message("Cells per tetramer per sample:")
print(df %>% count(tetramer, sample, sort = TRUE))

missing_tetramers <- setdiff(
  config$target_tetramers,
  unique(as.character(df$tetramer))
)

if (length(missing_tetramers) > 0) {
  warning(
    "These target tetramers were not found after filtering: ",
    paste(missing_tetramers, collapse = ", ")
  )
}


# ----------------------------- #
# 4. Pairwise participant sharing per tetramer
# ----------------------------- #

make_pairwise_sharing_private <- function(df, tet_name, min_pair_cells_shared = 3) {
  
  df_tet <- df %>%
    filter(as.character(tetramer) == tet_name) %>%
    filter(
      !is.na(sample), sample != "",
      !is.na(clonotype_for_sharing), clonotype_for_sharing != ""
    )
  
  if (nrow(df_tet) == 0) {
    return(NULL)
  }
  
  clones_by_sample <- df_tet %>%
    distinct(sample, clonotype_for_sharing)
  
  total_clones <- clones_by_sample %>%
    count(sample, name = "total_clones")
  
  samples <- sort(unique(clones_by_sample$sample))
  
  if (length(samples) == 0) {
    return(NULL)
  }
  
  clone_cells_sample <- df_tet %>%
    count(sample, clonotype_for_sharing, name = "n_cells")
  
  sets_named <- clones_by_sample %>%
    group_by(sample) %>%
    summarise(
      clones = list(unique(clonotype_for_sharing)),
      .groups = "drop"
    ) %>%
    { setNames(.$clones, .$sample) }
  
  pairs <- expand_grid(
    sample_1 = samples,
    sample_2 = samples
  )
  
  pair_shared_n <- function(a, b) {
    
    if (a == b) {
      return(NA_integer_)
    }
    
    shared <- intersect(sets_named[[a]], sets_named[[b]])
    
    if (length(shared) == 0) {
      return(0L)
    }
    
    a_counts <- clone_cells_sample %>%
      filter(sample == a, clonotype_for_sharing %in% shared) %>%
      select(clonotype_for_sharing, n_cells) %>%
      rename(n_a = n_cells)
    
    b_counts <- clone_cells_sample %>%
      filter(sample == b, clonotype_for_sharing %in% shared) %>%
      select(clonotype_for_sharing, n_cells) %>%
      rename(n_b = n_cells)
    
    pair_counts <- inner_join(
      a_counts,
      b_counts,
      by = "clonotype_for_sharing"
    ) %>%
      mutate(n_pair = n_a + n_b)
    
    sum(pair_counts$n_pair >= min_pair_cells_shared)
  }
  
  pair_df <- pairs %>%
    mutate(
      shared_clones_ge_cutoff =
        map2_int(sample_1, sample_2, pair_shared_n)
    ) %>%
    left_join(
      total_clones,
      by = c("sample_1" = "sample")
    ) %>%
    rename(total_clones_sample1 = total_clones) %>%
    mutate(
      shared_prop = if_else(
        is.na(shared_clones_ge_cutoff) |
          is.na(total_clones_sample1) |
          total_clones_sample1 == 0,
        NA_real_,
        shared_clones_ge_cutoff / total_clones_sample1
      ),
      shared_percent = shared_prop * 100,
      log10_shared_percent = if_else(
        shared_percent > 0,
        log10(shared_percent),
        NA_real_
      ),
      tetramer = tet_name
    )
  
  private_df <- pair_df %>%
    filter(sample_1 != sample_2) %>%
    group_by(sample_1) %>%
    summarise(
      private_prop = 1 - sum(shared_prop, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      private_prop = pmax(private_prop, 0),
      private_percent = private_prop * 100,
      log10_private_percent = if_else(
        private_percent > 0,
        log10(private_percent),
        NA_real_
      )
    )
  
  pair_df2 <- pair_df %>%
    left_join(private_df, by = "sample_1") %>%
    mutate(
      final_prop = if_else(
        sample_1 == sample_2,
        private_prop,
        shared_prop
      ),
      final_percent = final_prop * 100,
      final_log10_percent = if_else(
        final_percent > 0,
        log10(final_percent),
        NA_real_
      )
    )
  
  sample_cohort <- df_tet %>%
    distinct(sample, cohort)
  
  pair_df2 <- pair_df2 %>%
    left_join(sample_cohort, by = c("sample_1" = "sample")) %>%
    rename(cohort_1 = cohort) %>%
    left_join(sample_cohort, by = c("sample_2" = "sample")) %>%
    rename(cohort_2 = cohort)
  
  list(
    full = pair_df2,
    three_col = pair_df2 %>%
      transmute(
        sample_1,
        sample_2,
        log10_value = final_log10_percent
      )
  )
}

pairwise_private_dir <- file.path(config$outdir, "pairwise_private_by_tetramer")
create_dir(pairwise_private_dir)

walk(config$target_tetramers, function(tet) {
  
  out <- make_pairwise_sharing_private(
    df = df,
    tet_name = tet,
    min_pair_cells_shared = config$min_pair_cells_shared
  )
  
  if (is.null(out) || nrow(out$full) == 0) {
    message("Skipping ", tet, " pairwise-private output: no data.")
    return(invisible(NULL))
  }
  
  safe_tet <- safe_filename(tet)
  
  write_csv(
    out$full,
    file.path(
      pairwise_private_dir,
      paste0(safe_tet, "_pairwise_sharing_PRIVATE_FULL.csv")
    )
  )
  
  write_csv(
    out$three_col,
    file.path(
      pairwise_private_dir,
      paste0(safe_tet, "_pairwise_sharing_PRIVATE_3col.csv")
    )
  )
  
  message("Wrote pairwise-private files for: ", safe_tet)
})


# ----------------------------- #
# 5. Count clonotypes per tetramer/sample
# ----------------------------- #

clone_counts <- df %>%
  count(
    tetramer,
    sample,
    cohort,
    clonotype_for_sharing,
    name = "n_cells"
  )

clone_counts_path <- file.path(config$outdir, config$clone_counts_csv)
write_csv(clone_counts, clone_counts_path)


# ----------------------------- #
# 6. Define public clonotypes
# ----------------------------- #

public_clone_lookup <- clone_counts %>%
  group_by(tetramer, clonotype_for_sharing) %>%
  summarise(
    n_samples_with_clone = n_distinct(sample),
    n_cohorts_with_clone = n_distinct(cohort),
    total_cells_in_clone = sum(n_cells),
    .groups = "drop"
  ) %>%
  mutate(
    passes_cell_cutoff =
      total_cells_in_clone >= config$min_public_clone_cells,
    
    is_public =
      n_samples_with_clone >= 2 & passes_cell_cutoff,
    
    is_cross_cohort_public =
      n_samples_with_clone >= 2 &
      n_cohorts_with_clone >= 2 &
      passes_cell_cutoff
  )

public_clone_lookup_path <- file.path(
  config$outdir,
  config$public_clone_lookup_csv
)

write_csv(public_clone_lookup, public_clone_lookup_path)

message("Public clonotypes per tetramer using cell cutoff:")
print(
  public_clone_lookup %>%
    filter(is_public) %>%
    count(tetramer, sort = FALSE)
)

message("Cross-cohort public clonotypes per tetramer using cell cutoff:")
print(
  public_clone_lookup %>%
    filter(is_cross_cohort_public) %>%
    count(tetramer, sort = FALSE)
)


# ----------------------------- #
# 7. Sample-level public sharing table
# ----------------------------- #

sample_level_sharing <- clone_counts %>%
  left_join(
    public_clone_lookup %>%
      select(
        tetramer,
        clonotype_for_sharing,
        is_public,
        is_cross_cohort_public,
        n_samples_with_clone,
        n_cohorts_with_clone
      ),
    by = c("tetramer", "clonotype_for_sharing")
  ) %>%
  group_by(tetramer, sample, cohort) %>%
  summarise(
    total_cells = sum(n_cells),
    total_clonotypes = n_distinct(clonotype_for_sharing),
    
    public_cells = sum(n_cells[is_public], na.rm = TRUE),
    public_clonotypes = n_distinct(clonotype_for_sharing[is_public]),
    
    cross_cohort_public_cells =
      sum(n_cells[is_cross_cohort_public], na.rm = TRUE),
    cross_cohort_public_clonotypes =
      n_distinct(clonotype_for_sharing[is_cross_cohort_public]),
    
    public_cell_fraction = public_cells / total_cells,
    public_clonotype_fraction = public_clonotypes / total_clonotypes,
    
    cross_cohort_public_cell_fraction =
      cross_cohort_public_cells / total_cells,
    cross_cohort_public_clonotype_fraction =
      cross_cohort_public_clonotypes / total_clonotypes,
    
    has_public_clone = public_clonotypes > 0,
    has_cross_cohort_public_clone = cross_cohort_public_clonotypes > 0,
    
    public_burden_10pct = public_cell_fraction >= 0.10,
    public_burden_20pct = public_cell_fraction >= 0.20,
    
    cross_cohort_public_burden_10pct =
      cross_cohort_public_cell_fraction >= 0.10,
    cross_cohort_public_burden_20pct =
      cross_cohort_public_cell_fraction >= 0.20,
    
    .groups = "drop"
  ) %>%
  mutate(
    tetramer = factor(tetramer, levels = config$target_tetramers)
  ) %>%
  arrange(tetramer, sample)

sample_level_sharing_path <- file.path(
  config$outdir,
  config$sample_level_sharing_csv
)

write_csv(sample_level_sharing, sample_level_sharing_path)


# ----------------------------- #
# 8. Tetramer-level summary
# ----------------------------- #

tetramer_summary_base <- sample_level_sharing %>%
  group_by(tetramer) %>%
  summarise(
    n_samples = n(),
    total_cells = sum(total_cells),
    total_clonotypes = sum(total_clonotypes),
    
    median_cells_per_sample = median(total_cells),
    q1_cells_per_sample = unname(quantile(total_cells, 0.25)),
    q3_cells_per_sample = unname(quantile(total_cells, 0.75)),
    
    samples_with_public_clone = sum(has_public_clone),
    prop_samples_with_public_clone =
      samples_with_public_clone / n_samples,
    
    samples_with_cross_cohort_public_clone =
      sum(has_cross_cohort_public_clone),
    prop_samples_with_cross_cohort_public_clone =
      samples_with_cross_cohort_public_clone / n_samples,
    
    samples_with_public_burden_10pct = sum(public_burden_10pct),
    prop_samples_with_public_burden_10pct =
      samples_with_public_burden_10pct / n_samples,
    
    samples_with_public_burden_20pct = sum(public_burden_20pct),
    prop_samples_with_public_burden_20pct =
      samples_with_public_burden_20pct / n_samples,
    
    samples_with_cross_cohort_public_burden_10pct =
      sum(cross_cohort_public_burden_10pct),
    prop_samples_with_cross_cohort_public_burden_10pct =
      samples_with_cross_cohort_public_burden_10pct / n_samples,
    
    samples_with_cross_cohort_public_burden_20pct =
      sum(cross_cohort_public_burden_20pct),
    prop_samples_with_cross_cohort_public_burden_20pct =
      samples_with_cross_cohort_public_burden_20pct / n_samples,
    
    median_public_cell_fraction =
      median(public_cell_fraction, na.rm = TRUE),
    q1_public_cell_fraction =
      unname(quantile(public_cell_fraction, 0.25, na.rm = TRUE)),
    q3_public_cell_fraction =
      unname(quantile(public_cell_fraction, 0.75, na.rm = TRUE)),
    
    median_public_clonotype_fraction =
      median(public_clonotype_fraction, na.rm = TRUE),
    q1_public_clonotype_fraction =
      unname(quantile(public_clonotype_fraction, 0.25, na.rm = TRUE)),
    q3_public_clonotype_fraction =
      unname(quantile(public_clonotype_fraction, 0.75, na.rm = TRUE)),
    
    median_cross_cohort_public_cell_fraction =
      median(cross_cohort_public_cell_fraction, na.rm = TRUE),
    q1_cross_cohort_public_cell_fraction =
      unname(quantile(cross_cohort_public_cell_fraction, 0.25, na.rm = TRUE)),
    q3_cross_cohort_public_cell_fraction =
      unname(quantile(cross_cohort_public_cell_fraction, 0.75, na.rm = TRUE)),
    
    .groups = "drop"
  )

public_clone_ci <- map2_dfr(
  tetramer_summary_base$samples_with_public_clone,
  tetramer_summary_base$n_samples,
  binom_ci
) %>%
  rename(
    public_clone_prevalence_ci_lower = ci_lower,
    public_clone_prevalence_ci_upper = ci_upper
  )

cross_cohort_public_clone_ci <- map2_dfr(
  tetramer_summary_base$samples_with_cross_cohort_public_clone,
  tetramer_summary_base$n_samples,
  binom_ci
) %>%
  rename(
    cross_cohort_public_clone_prevalence_ci_lower = ci_lower,
    cross_cohort_public_clone_prevalence_ci_upper = ci_upper
  )

tetramer_summary <- bind_cols(
  tetramer_summary_base,
  public_clone_ci,
  cross_cohort_public_clone_ci
) %>%
  mutate(
    substantial_public_sharing = case_when(
      samples_with_public_burden_10pct >= 3 &
        prop_samples_with_public_clone >= 0.30 ~ "Yes",
      TRUE ~ "No"
    ),
    
    substantial_cross_cohort_public_sharing = case_when(
      samples_with_cross_cohort_public_burden_10pct >= 3 &
        prop_samples_with_cross_cohort_public_clone >= 0.30 ~ "Yes",
      TRUE ~ "No"
    ),
    
    tetramer = factor(tetramer, levels = config$target_tetramers)
  ) %>%
  arrange(tetramer)


# ----------------------------- #
# 9. One-vs-rest Fisher tests
# ----------------------------- #

fisher_one_vs_rest <- function(dat, target_tetramer, outcome_col, test_name) {
  
  dat2 <- dat %>%
    mutate(
      is_target = as.character(tetramer) == target_tetramer,
      outcome = .data[[outcome_col]]
    )
  
  target_yes <- sum(dat2$is_target & dat2$outcome, na.rm = TRUE)
  target_no <- sum(dat2$is_target & !dat2$outcome, na.rm = TRUE)
  
  other_yes <- sum(!dat2$is_target & dat2$outcome, na.rm = TRUE)
  other_no <- sum(!dat2$is_target & !dat2$outcome, na.rm = TRUE)
  
  mat <- matrix(
    c(target_yes, target_no, other_yes, other_no),
    nrow = 2,
    byrow = TRUE
  )
  
  rownames(mat) <- c(target_tetramer, "Other")
  colnames(mat) <- c("Yes", "No")
  
  if (sum(mat) == 0) {
    return(tibble(
      target_tetramer = target_tetramer,
      test = test_name,
      target_yes = target_yes,
      target_no = target_no,
      other_yes = other_yes,
      other_no = other_no,
      odds_ratio = NA_real_,
      p_value = NA_real_
    ))
  }
  
  ft <- fisher.test(mat, alternative = "greater")
  
  tibble(
    target_tetramer = target_tetramer,
    test = test_name,
    target_yes = target_yes,
    target_no = target_no,
    other_yes = other_yes,
    other_no = other_no,
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value
  )
}

fisher_has_public_clone <- map_dfr(
  config$target_tetramers,
  ~ fisher_one_vs_rest(
    dat = sample_level_sharing,
    target_tetramer = .x,
    outcome_col = "has_public_clone",
    test_name = "Fisher_has_public_clone"
  )
) %>%
  mutate(FDR = adjust_bh_without_na(p_value))

fisher_public_burden_10pct <- map_dfr(
  config$target_tetramers,
  ~ fisher_one_vs_rest(
    dat = sample_level_sharing,
    target_tetramer = .x,
    outcome_col = "public_burden_10pct",
    test_name = "Fisher_public_burden_10pct"
  )
) %>%
  mutate(FDR = adjust_bh_without_na(p_value))

fisher_cross_cohort_has_public_clone <- map_dfr(
  config$target_tetramers,
  ~ fisher_one_vs_rest(
    dat = sample_level_sharing,
    target_tetramer = .x,
    outcome_col = "has_cross_cohort_public_clone",
    test_name = "Fisher_has_cross_cohort_public_clone"
  )
) %>%
  mutate(FDR = adjust_bh_without_na(p_value))

fisher_cross_cohort_public_burden_10pct <- map_dfr(
  config$target_tetramers,
  ~ fisher_one_vs_rest(
    dat = sample_level_sharing,
    target_tetramer = .x,
    outcome_col = "cross_cohort_public_burden_10pct",
    test_name = "Fisher_cross_cohort_public_burden_10pct"
  )
) %>%
  mutate(FDR = adjust_bh_without_na(p_value))


# ----------------------------- #
# 10. One-vs-rest Wilcoxon tests
# ----------------------------- #

wilcox_one_vs_rest <- function(dat, target_tetramer, metric_col, test_name) {
  
  target_values <- dat %>%
    filter(as.character(tetramer) == target_tetramer) %>%
    pull(.data[[metric_col]])
  
  other_values <- dat %>%
    filter(as.character(tetramer) != target_tetramer) %>%
    pull(.data[[metric_col]])
  
  target_values <- target_values[!is.na(target_values)]
  other_values <- other_values[!is.na(other_values)]
  
  target_median <- median(target_values, na.rm = TRUE)
  other_median <- median(other_values, na.rm = TRUE)
  
  if (length(target_values) < 2 || length(other_values) < 2) {
    return(tibble(
      target_tetramer = target_tetramer,
      test = test_name,
      n_target = length(target_values),
      n_other = length(other_values),
      target_median = target_median,
      other_median = other_median,
      fold_difference = safe_fold(target_median, other_median),
      p_value = NA_real_
    ))
  }
  
  wt <- wilcox.test(
    target_values,
    other_values,
    alternative = "greater",
    exact = FALSE
  )
  
  tibble(
    target_tetramer = target_tetramer,
    test = test_name,
    n_target = length(target_values),
    n_other = length(other_values),
    target_median = target_median,
    other_median = other_median,
    fold_difference = safe_fold(target_median, other_median),
    p_value = wt$p.value
  )
}

wilcox_public_cell_fraction <- map_dfr(
  config$target_tetramers,
  ~ wilcox_one_vs_rest(
    dat = sample_level_sharing,
    target_tetramer = .x,
    metric_col = "public_cell_fraction",
    test_name = "Wilcoxon_public_cell_fraction"
  )
) %>%
  mutate(FDR = adjust_bh_without_na(p_value))

wilcox_cross_cohort_public_cell_fraction <- map_dfr(
  config$target_tetramers,
  ~ wilcox_one_vs_rest(
    dat = sample_level_sharing,
    target_tetramer = .x,
    metric_col = "cross_cohort_public_cell_fraction",
    test_name = "Wilcoxon_cross_cohort_public_cell_fraction"
  )
) %>%
  mutate(FDR = adjust_bh_without_na(p_value))


# ----------------------------- #
# 11. Pairwise sample overlap metrics
# ----------------------------- #

calc_pairwise_overlap <- function(dat) {
  
  tet <- unique(as.character(dat$tetramer))
  
  counts <- dat %>%
    count(sample, clonotype_for_sharing, name = "n_cells")
  
  sample_info <- dat %>%
    distinct(sample, cohort)
  
  wide <- counts %>%
    pivot_wider(
      names_from = clonotype_for_sharing,
      values_from = n_cells,
      values_fill = 0
    )
  
  if (nrow(wide) < 2) {
    return(tibble(
      tetramer = tet,
      sample1 = NA_character_,
      sample2 = NA_character_,
      cohort1 = NA_character_,
      cohort2 = NA_character_,
      pair_type = NA_character_,
      shared_clonotypes = NA_real_,
      union_clonotypes = NA_real_,
      jaccard = NA_real_,
      overlap_coefficient = NA_real_,
      morisita_horn = NA_real_
    ))
  }
  
  mat <- wide %>%
    select(-sample) %>%
    as.matrix()
  
  rownames(mat) <- wide$sample
  
  pairs <- combn(wide$sample, 2, simplify = FALSE)
  
  map_dfr(pairs, function(pair) {
    
    s1 <- pair[1]
    s2 <- pair[2]
    
    x <- mat[s1, ]
    y <- mat[s2, ]
    
    shared <- sum(x > 0 & y > 0)
    union <- sum(x > 0 | y > 0)
    min_size <- min(sum(x > 0), sum(y > 0))
    
    cohort1 <- sample_info %>%
      filter(sample == s1) %>%
      pull(cohort) %>%
      unique()
    
    cohort2 <- sample_info %>%
      filter(sample == s2) %>%
      pull(cohort) %>%
      unique()
    
    cohort1 <- cohort1[1]
    cohort2 <- cohort2[1]
    
    tibble(
      tetramer = tet,
      sample1 = s1,
      sample2 = s2,
      cohort1 = cohort1,
      cohort2 = cohort2,
      pair_type = ifelse(cohort1 == cohort2, "within_cohort", "cross_cohort"),
      shared_clonotypes = shared,
      union_clonotypes = union,
      jaccard = ifelse(union == 0, NA_real_, shared / union),
      overlap_coefficient = ifelse(min_size == 0, NA_real_, shared / min_size),
      morisita_horn = morisita_horn(x, y)
    )
  })
}

pairwise_overlap <- df %>%
  group_split(tetramer) %>%
  map_dfr(calc_pairwise_overlap) %>%
  mutate(
    tetramer = factor(tetramer, levels = config$target_tetramers)
  )

pairwise_summary <- pairwise_overlap %>%
  group_by(tetramer) %>%
  summarise(
    n_sample_pairs = sum(!is.na(sample1)),
    
    median_shared_clonotypes = median(shared_clonotypes, na.rm = TRUE),
    mean_shared_clonotypes = mean(shared_clonotypes, na.rm = TRUE),
    
    median_jaccard = median(jaccard, na.rm = TRUE),
    mean_jaccard = mean(jaccard, na.rm = TRUE),
    
    median_overlap_coefficient = median(overlap_coefficient, na.rm = TRUE),
    mean_overlap_coefficient = mean(overlap_coefficient, na.rm = TRUE),
    
    median_morisita_horn = median(morisita_horn, na.rm = TRUE),
    mean_morisita_horn = mean(morisita_horn, na.rm = TRUE),
    
    median_cross_cohort_shared_clonotypes =
      median(shared_clonotypes[pair_type == "cross_cohort"], na.rm = TRUE),
    median_cross_cohort_jaccard =
      median(jaccard[pair_type == "cross_cohort"], na.rm = TRUE),
    median_cross_cohort_morisita_horn =
      median(morisita_horn[pair_type == "cross_cohort"], na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    tetramer = factor(tetramer, levels = config$target_tetramers)
  ) %>%
  arrange(tetramer)

tetramer_summary_final <- tetramer_summary %>%
  left_join(pairwise_summary, by = "tetramer") %>%
  arrange(tetramer)


# ----------------------------- #
# 12. Combine statistical outputs
# ----------------------------- #

fisher_tests_all <- bind_rows(
  fisher_has_public_clone,
  fisher_public_burden_10pct,
  fisher_cross_cohort_has_public_clone,
  fisher_cross_cohort_public_burden_10pct
) %>%
  mutate(
    target_tetramer = factor(target_tetramer, levels = config$target_tetramers)
  ) %>%
  arrange(test, target_tetramer)

wilcox_tests_all <- bind_rows(
  wilcox_public_cell_fraction,
  wilcox_cross_cohort_public_cell_fraction
) %>%
  mutate(
    target_tetramer = factor(target_tetramer, levels = config$target_tetramers)
  ) %>%
  arrange(test, target_tetramer)


# ----------------------------- #
# 13. Pairwise tetramer tests and lower-triangle heatmap
# ----------------------------- #

pairwise_tetramer_tests <- combn(config$target_tetramers, 2, simplify = FALSE) %>%
  map_dfr(function(x) {
    
    tet1 <- x[1]
    tet2 <- x[2]
    
    values1 <- sample_level_sharing %>%
      filter(as.character(tetramer) == tet1) %>%
      pull(public_cell_fraction)
    
    values2 <- sample_level_sharing %>%
      filter(as.character(tetramer) == tet2) %>%
      pull(public_cell_fraction)
    
    values1 <- values1[!is.na(values1)]
    values2 <- values2[!is.na(values2)]
    
    med1 <- median(values1, na.rm = TRUE)
    med2 <- median(values2, na.rm = TRUE)
    
    pval <- if (length(values1) < 2 || length(values2) < 2) {
      NA_real_
    } else {
      wilcox.test(
        values1,
        values2,
        alternative = "two.sided",
        exact = FALSE
      )$p.value
    }
    
    tibble(
      tetramer_1 = tet1,
      tetramer_2 = tet2,
      n_1 = length(values1),
      n_2 = length(values2),
      median_1 = med1,
      median_2 = med2,
      median_difference = med1 - med2,
      p_value = pval
    )
  }) %>%
  mutate(
    FDR = adjust_bh_without_na(p_value),
    significance = case_when(
      is.na(FDR) ~ "",
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ ""
    ),
    FDR_label = case_when(
      is.na(FDR) ~ "NA",
      FDR < 0.001 ~ "FDR<0.001",
      TRUE ~ paste0("FDR=", signif(FDR, 2))
    )
  )

tetramer_positions <- tibble(
  tetramer = config$target_tetramers,
  idx = seq_along(config$target_tetramers)
)

pairwise_heatmap_df <- pairwise_tetramer_tests %>%
  left_join(
    tetramer_positions %>%
      rename(tetramer_1 = tetramer, idx_1 = idx),
    by = "tetramer_1"
  ) %>%
  left_join(
    tetramer_positions %>%
      rename(tetramer_2 = tetramer, idx_2 = idx),
    by = "tetramer_2"
  ) %>%
  mutate(
    tetramer_row = ifelse(idx_1 > idx_2, tetramer_1, tetramer_2),
    tetramer_col = ifelse(idx_1 > idx_2, tetramer_2, tetramer_1),
    
    median_row = ifelse(idx_1 > idx_2, median_1, median_2),
    median_col = ifelse(idx_1 > idx_2, median_2, median_1),
    
    median_difference = median_row - median_col,
    insufficient_data = is.na(p_value),
    
    label = case_when(
      insufficient_data ~ "",
      TRUE ~ significance
    ),
    
    tetramer_row = factor(tetramer_row, levels = rev(config$target_tetramers)),
    tetramer_col = factor(tetramer_col, levels = config$target_tetramers)
  )

bad_self_tests <- pairwise_heatmap_df %>%
  filter(as.character(tetramer_row) == as.character(tetramer_col))

if (nrow(bad_self_tests) > 0) {
  stop("Self-comparisons are present in the lower-triangle heatmap data.")
}

message("Pairwise comparisons with insufficient data:")
print(
  pairwise_heatmap_df %>%
    filter(insufficient_data) %>%
    select(tetramer_row, tetramer_col, n_1, n_2, p_value, FDR)
)

heatmap_cols <- c(
  "#005b96",
  "#6497b1",
  "#b3cde0",
  "#ffefea",
  "#fbd9d3",
  "#ffb09c",
  "#fe5757",
  "#cb2424",
  "#900000"
)

max_abs_diff <- max(abs(pairwise_heatmap_df$median_difference), na.rm = TRUE)

if (!is.finite(max_abs_diff) || max_abs_diff == 0) {
  max_abs_diff <- 1
}

p_pairwise_lower_heatmap <- ggplot(
  pairwise_heatmap_df,
  aes(
    x = tetramer_col,
    y = tetramer_row,
    fill = median_difference
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = label),
    size = 5,
    fontface = "bold",
    colour = "black"
  ) +
  scale_fill_gradientn(
    colours = heatmap_cols,
    values = scales::rescale(
      seq(-max_abs_diff, max_abs_diff, length.out = length(heatmap_cols))
    ),
    limits = c(-max_abs_diff, max_abs_diff),
    na.value = "grey90",
    name = "Median difference\n(row - column)"
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = "Column tetramer",
    y = "Row tetramer",
    title = "Pairwise comparison of public TRB clonotype burden",
    subtitle = paste0(
      "Lower triangle only; fill shows median difference in sample-level ",
      "public cell fraction; stars show BH-adjusted Wilcoxon FDR."
    )
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(face = "bold")
  )


# ----------------------------- #
# 14. Plots
# ----------------------------- #

p_public_fraction <- ggplot(
  tetramer_summary_final,
  aes(x = tetramer, y = median_public_cell_fraction)
) +
  geom_col() +
  theme_classic(base_size = 12) +
  labs(
    x = "Tetramer specificity",
    y = "Median sample-level public cell fraction",
    title = "Sample-level public clonotype burden across EBV specificities"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_public_prevalence <- ggplot(
  tetramer_summary_final,
  aes(x = tetramer, y = prop_samples_with_public_clone)
) +
  geom_col() +
  theme_classic(base_size = 12) +
  labs(
    x = "Tetramer specificity",
    y = "Proportion of samples with public clonotypes",
    title = "Prevalence of public clonotype sharing across EBV specificities"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_cross_cohort_fraction <- ggplot(
  tetramer_summary_final,
  aes(x = tetramer, y = median_cross_cohort_public_cell_fraction)
) +
  geom_col() +
  theme_classic(base_size = 12) +
  labs(
    x = "Tetramer specificity",
    y = "Median sample-level cross-cohort public cell fraction",
    title = "Cross-cohort public clonotype burden across EBV specificities"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_cross_cohort_prevalence <- ggplot(
  tetramer_summary_final,
  aes(x = tetramer, y = prop_samples_with_cross_cohort_public_clone)
) +
  geom_col() +
  theme_classic(base_size = 12) +
  labs(
    x = "Tetramer specificity",
    y = "Proportion of samples with cross-cohort public clonotypes",
    title = "Cross-cohort public clonotype sharing across EBV specificities"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ----------------------------- #
# 15. Save CSV outputs and plots
# ----------------------------- #

tetramer_summary_path <- file.path(config$outdir, config$tetramer_summary_csv)
pairwise_summary_path <- file.path(config$outdir, config$pairwise_summary_csv)
pairwise_overlap_path <- file.path(config$outdir, config$pairwise_overlap_csv)
fisher_tests_path <- file.path(config$outdir, config$fisher_tests_csv)
wilcox_tests_path <- file.path(config$outdir, config$wilcox_tests_csv)
pairwise_tetramer_tests_path <- file.path(config$outdir, config$pairwise_tetramer_tests_csv)
pairwise_heatmap_values_path <- file.path(config$outdir, config$pairwise_heatmap_values_csv)

write_csv(tetramer_summary_final, tetramer_summary_path)
write_csv(pairwise_summary, pairwise_summary_path)
write_csv(pairwise_overlap, pairwise_overlap_path)
write_csv(fisher_tests_all, fisher_tests_path)
write_csv(wilcox_tests_all, wilcox_tests_path)
write_csv(pairwise_tetramer_tests, pairwise_tetramer_tests_path)
write_csv(pairwise_heatmap_df, pairwise_heatmap_values_path)

public_fraction_plot_path <- file.path(config$outdir, config$public_fraction_plot)
public_prevalence_plot_path <- file.path(config$outdir, config$public_prevalence_plot)
cross_cohort_fraction_plot_path <- file.path(config$outdir, config$cross_cohort_fraction_plot)
cross_cohort_prevalence_plot_path <- file.path(config$outdir, config$cross_cohort_prevalence_plot)
pairwise_heatmap_plot_path <- file.path(config$outdir, config$pairwise_heatmap_plot)

ggsave(
  filename = public_fraction_plot_path,
  plot = p_public_fraction,
  width = config$plot_width,
  height = config$plot_height,
  dpi = config$plot_dpi
)

ggsave(
  filename = public_prevalence_plot_path,
  plot = p_public_prevalence,
  width = config$plot_width,
  height = config$plot_height,
  dpi = config$plot_dpi
)

ggsave(
  filename = cross_cohort_fraction_plot_path,
  plot = p_cross_cohort_fraction,
  width = config$plot_width,
  height = config$plot_height,
  dpi = config$plot_dpi
)

ggsave(
  filename = cross_cohort_prevalence_plot_path,
  plot = p_cross_cohort_prevalence,
  width = config$plot_width,
  height = config$plot_height,
  dpi = config$plot_dpi
)

ggsave(
  filename = pairwise_heatmap_plot_path,
  plot = p_pairwise_lower_heatmap,
  width = config$heatmap_width,
  height = config$heatmap_height,
  dpi = config$heatmap_dpi
)


# ----------------------------- #
# 16. Export Excel workbook
# ----------------------------- #

workbook_path <- file.path(config$outdir, config$workbook_file)

wb <- createWorkbook()

addWorksheet(wb, "tetramer_summary")
writeData(wb, "tetramer_summary", tetramer_summary_final)

addWorksheet(wb, "sample_level_sharing")
writeData(wb, "sample_level_sharing", sample_level_sharing)

addWorksheet(wb, "public_clone_lookup")
writeData(wb, "public_clone_lookup", public_clone_lookup)

addWorksheet(wb, "clone_counts")
writeData(wb, "clone_counts", clone_counts)

addWorksheet(wb, "pairwise_summary")
writeData(wb, "pairwise_summary", pairwise_summary)

addWorksheet(wb, "pairwise_overlap")
writeData(wb, "pairwise_overlap", pairwise_overlap)

addWorksheet(wb, "fisher_tests_all")
writeData(wb, "fisher_tests_all", fisher_tests_all)

addWorksheet(wb, "wilcox_tests_all")
writeData(wb, "wilcox_tests_all", wilcox_tests_all)

addWorksheet(wb, "pairwise_tetramer_tests")
writeData(wb, "pairwise_tetramer_tests", pairwise_tetramer_tests)

addWorksheet(wb, "pairwise_heatmap_values")
writeData(wb, "pairwise_heatmap_values", pairwise_heatmap_df)

saveWorkbook(wb, workbook_path, overwrite = TRUE)


# ----------------------------- #
# 17. Save session information
# ----------------------------- #

session_info_path <- file.path(config$outdir, config$session_info_file)

writeLines(
  c(
    "Configuration:",
    capture.output(str(config)),
    "",
    "Input source:",
    "Seurat RDS metadata from merged@meta.data",
    "",
    "Clonotype definition:",
    "TRB_cdr3 + TRB_v_gene + TRB_j_gene",
    "",
    "Participant/sample column used:",
    config$sample_col,
    "",
    "Cells retained after filtering:",
    as.character(nrow(df)),
    "",
    "Output files:",
    filtered_metadata_path,
    clone_counts_path,
    public_clone_lookup_path,
    sample_level_sharing_path,
    tetramer_summary_path,
    pairwise_summary_path,
    pairwise_overlap_path,
    fisher_tests_path,
    wilcox_tests_path,
    pairwise_tetramer_tests_path,
    pairwise_heatmap_values_path,
    public_fraction_plot_path,
    public_prevalence_plot_path,
    cross_cohort_fraction_plot_path,
    cross_cohort_prevalence_plot_path,
    pairwise_heatmap_plot_path,
    workbook_path,
    "",
    "Session information:",
    capture.output(sessionInfo())
  ),
  session_info_path
)


# ----------------------------- #
# 18. Print key tables
# ----------------------------- #

message("\nPublic TRB clonotype sharing workflow complete.")
message("Filtered metadata: ", filtered_metadata_path)
message("Tetramer summary: ", tetramer_summary_path)
message("Sample-level sharing: ", sample_level_sharing_path)
message("Public clone lookup: ", public_clone_lookup_path)
message("Pairwise-private files: ", pairwise_private_dir)
message("Workbook: ", workbook_path)
message("Session info: ", session_info_path)

message("\nFinal tetramer summary:")
print(tetramer_summary_final)

message("\nFisher tests:")
print(fisher_tests_all)

message("\nWilcoxon tests:")
print(wilcox_tests_all)
