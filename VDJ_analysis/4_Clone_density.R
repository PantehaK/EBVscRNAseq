#!/usr/bin/env Rscript

# ==============================================================================
# Clone density quantification across EBV antigen groups and T cell clusters
# ==============================================================================
#
# Purpose:
#   This script quantifies donor-level clonal density across T cell clusters and
#   antigen groups in the final baseline EBV Seurat object.
#
#   It:
#     1. Loads the final baseline EBV Seurat object.
#     2. Defines TRB clonotypes using:
#          TRB_cdr3 + TRB_v_gene + TRB_j_gene
#     3. Assigns antigen groups:
#          Negative, Latent, Lytic
#     4. Calculates clone size once per donor/sample across all antigen groups.
#     5. Calculates donor-level clone density per:
#          donor x cohort x cluster x antigen group
#     6. Runs pairwise tests for boxplots:
#          - MS vs Control within each cluster and antigen group
#          - paired antigen-group comparisons within each cluster and cohort
#     7. Runs cluster-vs-cluster tests for heatmap enrichment:
#          - within each cohort x antigen group panel
#     8. Saves:
#          - donor-level clone density table
#          - pairwise statistics
#          - cluster enrichment statistics
#          - heatmap source table
#          - heatmap PNG/PDF
#          - one-page-per-cluster boxplot PDF
#          - session information
#
# Notes:
#   - Paths are intentionally generic for GitHub.
#   - Update the config block before running.
#   - Clone size is calculated across antigen groups within each donor, matching
#     the original analysis intent.
#   - The default donor/sample column is `sample`. If your manuscript workflow
#     uses `id` as the donor identifier, change config$donor_col to "id".
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
  library(scales)
  library(readr)
  library(broom)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  
  # Input RDS
  input_rds = "path/to/Publication_RDS/Baseline_EBV_integrated_seurat.rds",
  
  # Output directory
  output_dir = "path/to/Publication_data/EBV_baseline/clone_density_quantification",
  
  # Metadata columns
  donor_col = "sample",
  cohort_col = "cohort",
  age_col = "age",
  sex_col = "sex",
  celltype_col = "celltype_new",
  virus_col = "virus",
  lifecycle_col = "lifecycle",
  latency_col = "latency",
  
  # TRB columns used to define clonotype
  trb_cdr3_col = "TRB_cdr3",
  trb_v_col = "TRB_v_gene",
  trb_j_col = "TRB_j_gene",
  
  # Orders
  celltype_order = c(
    "Naive/early TCM",
    "TEM",
    "Late TEM",
    "CTL",
    "CD69+ early activated T",
    "innate-like T",
    "CD69+ TEM"
  ),
  
  antigen_order = c("Negative", "Latent", "Lytic"),
  cohort_order = c("Control", "MS"),
  
  # Testing thresholds
  min_cells_for_boxplot_test = 5,
  min_cells_for_heatmap_test = 5,
  min_paired_donors_for_paired_test = 2,
  min_group_donors_for_unpaired_test = 2,
  
  # Plot settings
  control_colour = "#2166AC",
  ms_colour = "#B2182B",
  heatmap_low_colour = "#2166AC",
  heatmap_high_colour = "#B2182B",
  
  heatmap_width = 8,
  heatmap_height = 6,
  heatmap_dpi = 600,
  
  boxplot_pdf_width = 6.5,
  boxplot_pdf_height = 6,
  
  # Output files
  clone_density_csv = "cluster_clone_density_per_sample.csv",
  boxplot_tests_csv = "boxplot_pairwise_clone_density_tests_FDR_within_cluster.csv",
  cluster_tests_csv = "cluster_vs_cluster_clone_density_tests_FDR_within_cohort_antigen.csv",
  cluster_hits_csv = "clusters_significantly_enriched_for_clone_density.csv",
  heatmap_table_csv = "heatmap_clone_density_cluster_enrichment_by_cohort.csv",
  
  heatmap_png = "heatmap_cluster_enrichment_clone_density_by_cohort.png",
  heatmap_pdf = "heatmap_cluster_enrichment_clone_density_by_cohort.pdf",
  boxplot_pdf = "boxplot_clone_density_antigen_groups_by_cluster_one_page_each_MS_vs_Control_FDR_within_cluster.pdf",
  
  metadata_with_clone_density_csv = "baseline_metadata_with_TRB_clone_density.csv",
  session_info_file = "sessionInfo_clone_density_quantification.txt"
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


clean_vdj <- function(x) {
  
  x <- as.character(x)
  x <- str_trim(x)
  x[x %in% c("", "NA", "None", "none", "nan", "NaN", "NULL", "null")] <- NA_character_
  x
}


clean_cohort <- function(x, cohort_order = c("Control", "MS")) {
  
  x_raw <- as.character(x)
  x_clean <- str_to_lower(str_trim(x_raw))
  
  out <- case_when(
    x_clean %in% c("control", "ctrl", "nms", "non-ms", "nonms", "non_ms", "healthy control", "healthy controls", "hc") ~ "Control",
    x_clean %in% c("ms", "rrms", "pwms", "multiple sclerosis") ~ "MS",
    TRUE ~ x_raw
  )
  
  factor(out, levels = cohort_order)
}


mean_or_na <- function(x) {
  
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  
  mean(x, na.rm = TRUE)
}


median_or_na <- function(x) {
  
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  
  median(x, na.rm = TRUE)
}


safe_wilcox_formula <- function(formula, data) {
  
  tryCatch(
    suppressWarnings(
      wilcox.test(formula, data = data, exact = FALSE)$p.value
    ),
    error = function(e) NA_real_
  )
}


safe_wilcox_paired <- function(x, y) {
  
  tryCatch(
    suppressWarnings(
      wilcox.test(x, y, paired = TRUE, exact = FALSE)$p.value
    ),
    error = function(e) NA_real_
  )
}


adjust_bh_without_na <- function(p) {
  
  out <- rep(NA_real_, length(p))
  idx <- !is.na(p)
  
  if (any(idx)) {
    out[idx] <- p.adjust(p[idx], method = "BH")
  }
  
  out
}


format_raw_p <- function(p) {
  
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "p<0.001",
    TRUE ~ paste0("p=", signif(p, 2))
  )
}


format_fdr <- function(fdr) {
  
  case_when(
    is.na(fdr) ~ "",
    fdr < 0.001 ~ "FDR<0.001",
    TRUE ~ paste0("FDR=", signif(fdr, 2))
  )
}


fdr_to_stars <- function(fdr) {
  
  case_when(
    is.na(fdr) ~ "",
    fdr < 0.001 ~ "***",
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  )
}


assign_antigen_group <- function(virus, lifecycle, latency) {
  
  virus_lower <- str_to_lower(as.character(virus))
  lifecycle_lower <- str_to_lower(as.character(lifecycle))
  latency_lower <- str_to_lower(as.character(latency))
  
  case_when(
    virus_lower %in% c("negative", "neg", "non-ebv", "non_ebv") ~ "Negative",
    
    lifecycle_lower == "lytic" ~ "Lytic",
    latency_lower == "lytic" ~ "Lytic",
    
    virus_lower == "ebv" &
      (
        lifecycle_lower %in% c("latent", "latency") |
          as.character(latency) %in% c("Latency I/0", "Latency II", "Latency III", "Latent")
      ) ~ "Latent",
    
    TRUE ~ NA_character_
  )
}


get_x_position <- function(antigen, cohort, antigen_order) {
  
  antigen_index <- match(as.character(antigen), antigen_order)
  cohort_offset <- ifelse(as.character(cohort) == "Control", -0.1875, 0.1875)
  
  antigen_index + cohort_offset
}


make_plot_annotations <- function(test_tbl, plot_df, antigen_order) {
  
  sig_df <- test_tbl %>%
    filter(
      !is.na(FDR),
      FDR < 0.05
    )
  
  if (nrow(sig_df) == 0) {
    return(tibble())
  }
  
  y_range <- range(plot_df$mean_log_clone_size, na.rm = TRUE)
  y_span <- diff(y_range)
  
  if (!is.finite(y_span) || y_span == 0) {
    y_span <- 0.1
  }
  
  y_start <- max(plot_df$mean_log_clone_size, na.rm = TRUE) + y_span * 0.12
  y_step <- y_span * 0.12
  
  sig_df %>%
    arrange(
      comparison_type,
      group1_antigen,
      group2_antigen,
      group1_cohort,
      group2_cohort
    ) %>%
    mutate(
      x1 = get_x_position(group1_antigen, group1_cohort, antigen_order),
      x2 = get_x_position(group2_antigen, group2_cohort, antigen_order),
      xmid = (x1 + x2) / 2,
      y = y_start + (row_number() - 1) * y_step,
      y_tick = y - y_span * 0.03,
      label = sig_label_FDR
    )
}


# ----------------------------- #
# 3. Load object and prepare metadata
# ----------------------------- #

create_dir(config$output_dir)

message("Loading baseline EBV object:")
message(config$input_rds)

merged <- load_seurat_object(config$input_rds)

meta <- merged@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell")

required_cols <- c(
  config$donor_col,
  config$cohort_col,
  config$age_col,
  config$sex_col,
  config$celltype_col,
  config$virus_col,
  config$lifecycle_col,
  config$latency_col,
  config$trb_cdr3_col,
  config$trb_v_col,
  config$trb_j_col
)

check_required_columns(
  meta,
  required_cols = required_cols,
  object_name = "baseline EBV metadata"
)


# ----------------------------- #
# 4. Create TRB clonotype ID and antigen group
# ----------------------------- #

df <- meta %>%
  mutate(
    donor_id = as.character(.data[[config$donor_col]]),
    cohort = clean_cohort(.data[[config$cohort_col]], config$cohort_order),
    age = suppressWarnings(as.numeric(.data[[config$age_col]])),
    sex = as.character(.data[[config$sex_col]]),
    celltype_new = as.character(.data[[config$celltype_col]]),
    
    TRB_cdr3 = clean_vdj(.data[[config$trb_cdr3_col]]),
    TRB_v_gene = clean_vdj(.data[[config$trb_v_col]]),
    TRB_j_gene = clean_vdj(.data[[config$trb_j_col]]),
    
    clone_id = case_when(
      !is.na(TRB_cdr3) & !is.na(TRB_v_gene) & !is.na(TRB_j_gene) ~
        paste(TRB_cdr3, TRB_v_gene, TRB_j_gene, sep = "|"),
      TRUE ~ NA_character_
    ),
    
    antigen_group = assign_antigen_group(
      virus = .data[[config$virus_col]],
      lifecycle = .data[[config$lifecycle_col]],
      latency = .data[[config$latency_col]]
    ),
    
    cohort = factor(cohort, levels = config$cohort_order),
    antigen_group = factor(antigen_group, levels = config$antigen_order),
    celltype_new = factor(celltype_new, levels = config$celltype_order)
  ) %>%
  filter(
    !is.na(donor_id), donor_id != "",
    !is.na(cohort),
    !is.na(celltype_new),
    !is.na(clone_id),
    !is.na(antigen_group)
  )

message("\nAntigen group counts:")
print(table(df$antigen_group, useNA = "ifany"))

message("\nCohort by antigen group:")
print(table(df$cohort, df$antigen_group, useNA = "ifany"))

message("\nCelltype by antigen group:")
print(table(df$celltype_new, df$antigen_group, useNA = "ifany"))


# ----------------------------- #
# 5. Calculate clone size once per donor
# ----------------------------- #
#
# Important:
#   clone_size is calculated across all antigen groups within each donor.
#   Do not calculate clone size separately inside Negative / Latent / Lytic.
#
# ----------------------------- #

df <- df %>%
  group_by(donor_id, clone_id) %>%
  mutate(
    clone_size = n(),
    log_clone_size = log(clone_size + 1)
  ) %>%
  ungroup()

clone_meta <- df %>%
  select(cell, clone_id, clone_size, log_clone_size, antigen_group) %>%
  column_to_rownames("cell")

merged <- AddMetaData(
  merged,
  metadata = clone_meta
)

metadata_with_clone_density_path <- file.path(
  config$output_dir,
  config$metadata_with_clone_density_csv
)

write_csv(
  df,
  metadata_with_clone_density_path
)


# ----------------------------- #
# 6. Donor-level clonal density
# ----------------------------- #

cluster_clone_density <- df %>%
  group_by(donor_id, cohort, age, sex, celltype_new, antigen_group) %>%
  summarise(
    n_cells = n(),
    mean_log_clone_size = mean(log_clone_size, na.rm = TRUE),
    median_log_clone_size = median(log_clone_size, na.rm = TRUE),
    total_log_clone_density = sum(log_clone_size, na.rm = TRUE),
    .groups = "drop"
  )

clone_density_path <- file.path(
  config$output_dir,
  config$clone_density_csv
)

write_csv(
  cluster_clone_density,
  clone_density_path
)


# ----------------------------- #
# 7. Boxplot statistics
# ----------------------------- #
#
# These are for the boxplots.
# FDR family = within each T-cell cluster.
#
# ----------------------------- #

test_df <- cluster_clone_density %>%
  filter(n_cells >= config$min_cells_for_boxplot_test) %>%
  mutate(
    celltype_new = factor(celltype_new, levels = config$celltype_order),
    antigen_group = factor(antigen_group, levels = config$antigen_order),
    cohort = factor(cohort, levels = config$cohort_order)
  ) %>%
  filter(
    !is.na(celltype_new),
    !is.na(antigen_group),
    !is.na(cohort),
    !is.na(mean_log_clone_size)
  )

# 7A. MS vs Control within the same antigen group and cluster.
same_antigen_MS_vs_Control <- test_df %>%
  group_by(celltype_new, antigen_group) %>%
  group_modify(~{
    
    dat <- .x
    ag <- as.character(.y$antigen_group)
    
    control_vals <- dat$mean_log_clone_size[dat$cohort == "Control"]
    ms_vals <- dat$mean_log_clone_size[dat$cohort == "MS"]
    
    n_control <- sum(dat$cohort == "Control")
    n_ms <- sum(dat$cohort == "MS")
    
    raw_p <- if (
      n_control >= config$min_group_donors_for_unpaired_test &&
      n_ms >= config$min_group_donors_for_unpaired_test
    ) {
      safe_wilcox_formula(mean_log_clone_size ~ cohort, dat)
    } else {
      NA_real_
    }
    
    tibble(
      comparison_type = "MS_vs_Control_same_antigen",
      comparison_label = paste0(ag, ": MS vs Control"),
      
      group1_cohort = "Control",
      group1_antigen = ag,
      group2_cohort = "MS",
      group2_antigen = ag,
      
      n_group1 = n_control,
      n_group2 = n_ms,
      n_paired = NA_integer_,
      
      mean_group1 = mean_or_na(control_vals),
      mean_group2 = mean_or_na(ms_vals),
      median_group1 = median_or_na(control_vals),
      median_group2 = median_or_na(ms_vals),
      
      delta_group2_minus_group1 =
        mean_or_na(ms_vals) - mean_or_na(control_vals),
      raw_p_value = raw_p
    )
  }) %>%
  ungroup() %>%
  select(-antigen_group)

# 7B. Paired antigen-group comparisons within the same cohort and cluster.
antigen_pairs <- list(
  c("Negative", "Latent"),
  c("Negative", "Lytic"),
  c("Latent", "Lytic")
)

within_cohort_antigen_tests <- test_df %>%
  group_by(celltype_new, cohort) %>%
  group_modify(~{
    
    dat <- .x
    co <- as.character(.y$cohort)
    
    bind_rows(lapply(antigen_pairs, function(pair) {
      
      ag1 <- pair[1]
      ag2 <- pair[2]
      
      dat1 <- dat %>%
        filter(antigen_group == ag1) %>%
        select(donor_id, value1 = mean_log_clone_size)
      
      dat2 <- dat %>%
        filter(antigen_group == ag2) %>%
        select(donor_id, value2 = mean_log_clone_size)
      
      paired_dat <- inner_join(dat1, dat2, by = "donor_id")
      
      n_group1 <- n_distinct(dat1$donor_id)
      n_group2 <- n_distinct(dat2$donor_id)
      n_paired <- nrow(paired_dat)
      
      raw_p <- if (n_paired >= config$min_paired_donors_for_paired_test) {
        safe_wilcox_paired(paired_dat$value1, paired_dat$value2)
      } else {
        NA_real_
      }
      
      tibble(
        comparison_type = "Within_cohort_antigen_groups_paired",
        comparison_label = paste0(co, ": ", ag2, " vs ", ag1),
        
        group1_cohort = co,
        group1_antigen = ag1,
        group2_cohort = co,
        group2_antigen = ag2,
        
        n_group1 = n_group1,
        n_group2 = n_group2,
        n_paired = n_paired,
        
        mean_group1 = mean_or_na(dat1$value1),
        mean_group2 = mean_or_na(dat2$value2),
        median_group1 = median_or_na(dat1$value1),
        median_group2 = median_or_na(dat2$value2),
        
        delta_group2_minus_group1 =
          mean_or_na(dat2$value2) - mean_or_na(dat1$value1),
        raw_p_value = raw_p
      )
    }))
  }) %>%
  ungroup() %>%
  select(-cohort)

# 7C. Combine boxplot tests and adjust FDR within each cluster.
boxplot_pairwise_tests <- bind_rows(
  same_antigen_MS_vs_Control,
  within_cohort_antigen_tests
) %>%
  group_by(celltype_new) %>%
  mutate(
    FDR = adjust_bh_without_na(raw_p_value)
  ) %>%
  ungroup() %>%
  mutate(
    sig_label_FDR = fdr_to_stars(FDR),
    raw_p_label = format_raw_p(raw_p_value),
    FDR_label = format_fdr(FDR)
  )

boxplot_tests_path <- file.path(
  config$output_dir,
  config$boxplot_tests_csv
)

write_csv(
  boxplot_pairwise_tests,
  boxplot_tests_path
)

print(boxplot_pairwise_tests)


# ----------------------------- #
# 8. Heatmap statistics
# ----------------------------- #
#
# Within each cohort x antigen group, compare cluster vs cluster clone density.
# FDR family = within each cohort x antigen_group panel.
#
# ----------------------------- #

heatmap_input <- cluster_clone_density %>%
  filter(n_cells >= config$min_cells_for_heatmap_test) %>%
  mutate(
    celltype_new = factor(celltype_new, levels = config$celltype_order),
    antigen_group = factor(antigen_group, levels = config$antigen_order),
    cohort = factor(cohort, levels = config$cohort_order)
  ) %>%
  filter(
    !is.na(donor_id),
    !is.na(cohort),
    !is.na(antigen_group),
    !is.na(celltype_new),
    !is.na(mean_log_clone_size)
  )

cluster_cluster_tests <- heatmap_input %>%
  group_by(cohort, antigen_group) %>%
  group_modify(~{
    
    dat <- .x
    
    present_clusters <- dat %>%
      pull(celltype_new) %>%
      as.character() %>%
      unique()
    
    present_clusters <- config$celltype_order[
      config$celltype_order %in% present_clusters
    ]
    
    if (length(present_clusters) < 2) {
      return(tibble())
    }
    
    cluster_pairs <- combn(present_clusters, 2, simplify = FALSE)
    
    bind_rows(lapply(cluster_pairs, function(pair) {
      
      cl1 <- pair[1]
      cl2 <- pair[2]
      
      dat1 <- dat %>%
        filter(celltype_new == cl1) %>%
        select(donor_id, value1 = mean_log_clone_size)
      
      dat2 <- dat %>%
        filter(celltype_new == cl2) %>%
        select(donor_id, value2 = mean_log_clone_size)
      
      paired_dat <- inner_join(dat1, dat2, by = "donor_id")
      
      n_cl1 <- n_distinct(dat1$donor_id)
      n_cl2 <- n_distinct(dat2$donor_id)
      n_paired <- nrow(paired_dat)
      
      raw_p <- if (n_paired >= config$min_paired_donors_for_paired_test) {
        safe_wilcox_paired(paired_dat$value1, paired_dat$value2)
      } else {
        NA_real_
      }
      
      tibble(
        cluster_1 = cl1,
        cluster_2 = cl2,
        n_cluster_1 = n_cl1,
        n_cluster_2 = n_cl2,
        n_paired = n_paired,
        
        mean_cluster_1 = mean_or_na(dat1$value1),
        mean_cluster_2 = mean_or_na(dat2$value2),
        median_cluster_1 = median_or_na(dat1$value1),
        median_cluster_2 = median_or_na(dat2$value2),
        
        mean_paired_difference_cluster_2_minus_cluster_1 =
          mean_or_na(paired_dat$value2 - paired_dat$value1),
        
        raw_p_value = raw_p
      )
    }))
  }) %>%
  ungroup() %>%
  group_by(cohort, antigen_group) %>%
  mutate(
    FDR = adjust_bh_without_na(raw_p_value)
  ) %>%
  ungroup() %>%
  mutate(
    sig_label_FDR = fdr_to_stars(FDR),
    raw_p_label = format_raw_p(raw_p_value),
    FDR_label = format_fdr(FDR)
  )

cluster_tests_path <- file.path(
  config$output_dir,
  config$cluster_tests_csv
)

write_csv(
  cluster_cluster_tests,
  cluster_tests_path
)

print(cluster_cluster_tests)


# ----------------------------- #
# 9. Identify cluster enrichment hits
# ----------------------------- #

cluster_enrichment_hits <- cluster_cluster_tests %>%
  filter(
    !is.na(FDR),
    FDR < 0.05,
    !is.na(mean_paired_difference_cluster_2_minus_cluster_1),
    mean_paired_difference_cluster_2_minus_cluster_1 != 0
  ) %>%
  mutate(
    enriched_cluster = if_else(
      mean_paired_difference_cluster_2_minus_cluster_1 > 0,
      cluster_2,
      cluster_1
    )
  ) %>%
  group_by(cohort, antigen_group, enriched_cluster) %>%
  summarise(
    n_significant_pairwise_wins = n(),
    best_FDR = min(FDR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sig_label_enrichment = fdr_to_stars(best_FDR)
  ) %>%
  rename(celltype_new = enriched_cluster)

cluster_hits_path <- file.path(
  config$output_dir,
  config$cluster_hits_csv
)

write_csv(
  cluster_enrichment_hits,
  cluster_hits_path
)


# ----------------------------- #
# 10. Build and save heatmap table
# ----------------------------- #

cluster_density_heatmap <- heatmap_input %>%
  group_by(cohort, antigen_group, celltype_new) %>%
  summarise(
    n_donors = n_distinct(donor_id),
    mean_density = mean(mean_log_clone_size, na.rm = TRUE),
    median_density = median(mean_log_clone_size, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(cohort, antigen_group) %>%
  mutate(
    panel_mean_density = mean(mean_density, na.rm = TRUE),
    panel_sd_density = sd(mean_density, na.rm = TRUE),
    density_z_within_panel = case_when(
      is.na(panel_sd_density) ~ 0,
      panel_sd_density == 0 ~ 0,
      TRUE ~ (mean_density - panel_mean_density) / panel_sd_density
    )
  ) %>%
  ungroup() %>%
  left_join(
    cluster_enrichment_hits,
    by = c("cohort", "antigen_group", "celltype_new")
  ) %>%
  mutate(
    n_significant_pairwise_wins = replace_na(n_significant_pairwise_wins, 0L),
    sig_label_enrichment = replace_na(sig_label_enrichment, ""),
    mean_density_label = round(mean_density, 2),
    celltype_new = factor(celltype_new, levels = rev(config$celltype_order)),
    antigen_group = factor(antigen_group, levels = config$antigen_order),
    cohort = factor(cohort, levels = config$cohort_order)
  )

heatmap_table_path <- file.path(
  config$output_dir,
  config$heatmap_table_csv
)

write_csv(
  cluster_density_heatmap,
  heatmap_table_path
)


# ----------------------------- #
# 11. Plot heatmap
# ----------------------------- #

p_cluster_enrichment_heatmap <- ggplot(
  cluster_density_heatmap,
  aes(
    x = antigen_group,
    y = celltype_new,
    fill = density_z_within_panel
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.45
  ) +
  geom_text(
    aes(label = mean_density_label),
    size = 3.2,
    vjust = 1.35
  ) +
  geom_text(
    aes(label = sig_label_enrichment),
    size = 7,
    vjust = -0.15
  ) +
  facet_wrap(~ cohort, nrow = 1) +
  scale_fill_gradient2(
    low = config$heatmap_low_colour,
    mid = "white",
    high = config$heatmap_high_colour,
    midpoint = 0,
    na.value = "grey90",
    name = "Relative clone density\nz-score within cohort/antigen"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Cluster enrichment of expanded clones",
    subtitle = paste0(
      "Tile values show mean donor-level log(clone size + 1). ",
      "Colour shows relative density within each cohort-antigen group. ",
      "Asterisks indicate FDR < 0.05 for cluster-vs-cluster comparisons within that panel."
    )
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 11),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 13),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10),
    legend.title = element_text(size = 10)
  )

heatmap_png_path <- file.path(
  config$output_dir,
  config$heatmap_png
)

heatmap_pdf_path <- file.path(
  config$output_dir,
  config$heatmap_pdf
)

ggsave(
  filename = heatmap_png_path,
  plot = p_cluster_enrichment_heatmap,
  width = config$heatmap_width,
  height = config$heatmap_height,
  dpi = config$heatmap_dpi
)

ggsave(
  filename = heatmap_pdf_path,
  plot = p_cluster_enrichment_heatmap,
  width = config$heatmap_width,
  height = config$heatmap_height
)


# ----------------------------- #
# 12. Boxplot PDF: one page per cluster
# ----------------------------- #

cluster_clone_density <- cluster_clone_density %>%
  mutate(
    celltype_new = factor(celltype_new, levels = config$celltype_order),
    antigen_group = factor(antigen_group, levels = config$antigen_order),
    cohort = factor(cohort, levels = config$cohort_order)
  )

cluster_levels <- config$celltype_order[
  config$celltype_order %in% unique(as.character(cluster_clone_density$celltype_new))
]

boxplot_pdf_path <- file.path(
  config$output_dir,
  config$boxplot_pdf
)

pdf(
  boxplot_pdf_path,
  width = config$boxplot_pdf_width,
  height = config$boxplot_pdf_height,
  useDingbats = FALSE
)

for (cl in cluster_levels) {
  
  plot_df <- cluster_clone_density %>%
    filter(celltype_new == cl)
  
  stat_df <- boxplot_pairwise_tests %>%
    filter(celltype_new == cl)
  
  ann_df <- make_plot_annotations(
    test_tbl = stat_df,
    plot_df = plot_df,
    antigen_order = config$antigen_order
  )
  
  p <- ggplot(
    plot_df,
    aes(
      x = antigen_group,
      y = mean_log_clone_size,
      fill = cohort,
      color = cohort
    )
  ) +
    geom_boxplot(
      outlier.shape = NA,
      width = 0.65,
      alpha = 0.35,
      position = position_dodge(width = 0.75)
    ) +
    geom_point(
      position = position_jitterdodge(
        jitter.width = 0.12,
        dodge.width = 0.75
      ),
      size = 2.2,
      alpha = 0.85
    ) +
    scale_fill_manual(
      values = c(
        "Control" = config$control_colour,
        "MS" = config$ms_colour
      ),
      drop = FALSE
    ) +
    scale_color_manual(
      values = c(
        "Control" = config$control_colour,
        "MS" = config$ms_colour
      ),
      drop = FALSE
    ) +
    labs(
      x = NULL,
      y = "Donor-level mean log(clone size + 1)",
      fill = NULL,
      color = NULL,
      title = cl,
      subtitle = "Asterisks show FDR < 0.05 corrected within this cluster"
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10),
      legend.position = "top"
    )
  
  if (nrow(ann_df) > 0) {
    
    y_range <- range(plot_df$mean_log_clone_size, na.rm = TRUE)
    y_span <- diff(y_range)
    
    if (!is.finite(y_span) || y_span == 0) {
      y_span <- 0.1
    }
    
    y_upper <- max(ann_df$y, na.rm = TRUE) + y_span * 0.15
    
    p <- p +
      geom_segment(
        data = ann_df,
        aes(x = x1, xend = x2, y = y, yend = y),
        inherit.aes = FALSE,
        linewidth = 0.45
      ) +
      geom_segment(
        data = ann_df,
        aes(x = x1, xend = x1, y = y_tick, yend = y),
        inherit.aes = FALSE,
        linewidth = 0.45
      ) +
      geom_segment(
        data = ann_df,
        aes(x = x2, xend = x2, y = y_tick, yend = y),
        inherit.aes = FALSE,
        linewidth = 0.45
      ) +
      geom_text(
        data = ann_df,
        aes(x = xmid, y = y + y_span * 0.025, label = label),
        inherit.aes = FALSE,
        size = 8
      ) +
      coord_cartesian(ylim = c(NA, y_upper))
  }
  
  print(p)
}

dev.off()


# ----------------------------- #
# 13. Save session information
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
    "Clone definition:",
    "TRB_cdr3 + TRB_v_gene + TRB_j_gene",
    "",
    "Clone size calculation:",
    "Calculated once across all antigen groups within each donor_id.",
    "",
    "Output files:",
    metadata_with_clone_density_path,
    clone_density_path,
    boxplot_tests_path,
    cluster_tests_path,
    cluster_hits_path,
    heatmap_table_path,
    heatmap_png_path,
    heatmap_pdf_path,
    boxplot_pdf_path,
    "",
    "Session information:",
    capture.output(sessionInfo())
  ),
  session_info_path
)


# ----------------------------- #
# 14. Completion messages
# ----------------------------- #

message("\nClone density quantification workflow complete.")
message("Metadata with clone density: ", metadata_with_clone_density_path)
message("Donor-level clone density: ", clone_density_path)
message("Boxplot statistics: ", boxplot_tests_path)
message("Cluster-vs-cluster statistics: ", cluster_tests_path)
message("Cluster enrichment hits: ", cluster_hits_path)
message("Heatmap source table: ", heatmap_table_path)
message("Heatmap PNG: ", heatmap_png_path)
message("Heatmap PDF: ", heatmap_pdf_path)
message("Boxplot PDF: ", boxplot_pdf_path)
message("Session info: ", session_info_path)
