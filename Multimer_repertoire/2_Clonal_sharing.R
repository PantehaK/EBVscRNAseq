
library(tidyverse)
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readr)

# --------------
#### In this script, we will create a table that determines shared clones between samples for each multimer specificity. This data will then be used in https://flourish.studio/ to create chord diagrams (circos plots) for each tetramer using the proportion of shared clones values per sample.

################## Clonal sharing between MS and Control samples per tetramer ###############
df <- read_xlsx(
  "L:/Lab_CoreyS/Smith Lab - Manuscripts/EBV MS manuscript/Data/Single cell RNA seq data/EBV T cell baseline/Baseline_EBV_TCR_data_module_scored.xlsx"
)

# -----------------------------
# 2. Define clone key (sample-specific)
# -----------------------------
df <- df %>%
  mutate(
    clone_key = paste(sample, TRB_cdr3, TRB_v_gene, TRB_j_gene, sep = "|")
  )

# -----------------------------
# 3. Assign numeric clone IDs
# -----------------------------
clone_map <- df %>%
  distinct(clone_key) %>%
  mutate(clone_id = row_number())

df <- df %>%
  left_join(clone_map, by = "clone_key")

# -----------------------------
# 4. Clone frequency WITHIN id x tetramer
# -----------------------------
df_tet <- df %>%
  filter(!is.na(tetramer), tetramer != "")

clone_freq <- df_tet %>%
  group_by(id, tetramer, clone_id) %>%
  summarise(clone_frequency = n(), .groups = "drop")

# total tetramer-specific cells per id (DENOMINATOR)
id_tet_totals <- df_tet %>%
  count(id, tetramer, name = "total_cells_id_tet")

# -----------------------------
# 5. Proportion per id x tetramer
# -----------------------------
clone_freq <- clone_freq %>%
  left_join(id_tet_totals, by = c("id", "tetramer")) %>%
  mutate(clone_proportion = clone_frequency / total_cells_id_tet)

# -----------------------------
# 6. Expansion category with min cells gate
#    If total_cells_id_tet < 3 => force Small
# -----------------------------
clone_freq <- clone_freq %>%
  mutate(
    clonal_expansion = case_when(
      clone_frequency == 1 ~ "Small expansion x < 10%",
      clone_proportion > 0.1 & clone_proportion <= 0.5 ~ "Medium expansion 10% < x < 50%",
      clone_proportion > 0.5 ~ "Hyperexpanded x > 50%",
      TRUE ~ "Small expansion x < 10%"
    )
  )

# -----------------------------
# 7. Join back (include tetramer in join keys)
# -----------------------------
df_final <- df %>%
  left_join(
    clone_freq %>%
      select(id, tetramer, clone_id, clone_frequency, clone_proportion, total_cells_id_tet, clonal_expansion),
    by = c("id", "tetramer", "clone_id")
  ) %>%
  select(-clone_key)

# -----------------------------
# Optional sanity check
# -----------------------------
# How many (id,tetramer) have <3 cells?
df_final %>%
  filter(!is.na(tetramer), tetramer != "") %>%
  distinct(id, tetramer, total_cells_id_tet) %>%
  summarise(n_lt3 = sum(total_cells_id_tet < 3), n_total = n())



library(tidyverse)

# --------------------------
# Inputs
# --------------------------
tet_order <- c(
  "GLCT", "RAKF",
  "RPPI", "FLRG", "QAKW",
  "AVFD*", "LLDF", "YLQQ",
  "CLGG", "FLYA",
  "RPQK*", "YNLR*"
)

antigens_keep <- c("BZLF1","BMLF1","EBNA3A","EBNA3C","LMP1","LMP2","EBNA1")

exp_levels <- c(
  "Hyperexpanded x > 50%",
  "Medium expansion 10% < x < 50%",
  "Small expansion x < 10%"
)

# Colors
fill_control <- c(
  "Small expansion x < 10%" = "#cfe6ff",
  "Medium expansion 10% < x < 50%" = "#3f9ef8",
  "Hyperexpanded x > 50%" = "#1f5fb8"
)

fill_ms <- c(
  "Small expansion x < 10%" = "#ffbaba",
  "Medium expansion 10% < x < 50%" = "#ff5252",
  "Hyperexpanded x > 50%" = "#a70000"
)

df_plot <- df_final %>%
  mutate(
    tetramer = factor(tetramer, levels = tet_order),
    antigen  = factor(antigen, levels = antigens_keep),
    clonal_expansion = factor(clonal_expansion, levels = exp_levels),
    cohort = factor(cohort, levels = c("Control","MS"))
  ) %>%
  filter(
    tetramer %in% tet_order,
    antigen %in% antigens_keep,
    !is.na(clonal_expansion),
    !is.na(cohort),
    str_starts(batch, "GEMEBV")          # <-- NEW: only GEMEBV batches
  )

# 1. Beta chain

# ============================================================
# Pairwise clonal sharing (per tetramer), clones >= 3 cells only
# Output: one CSV per tetramer (FULL + 3-column version)
# Sharing metric: Jaccard = |A ∩ B| / |A ∪ B|  (includes private clones)
# Log transform: log10(percent shared) by default for the 3-col file
# ============================================================
------------
  # Inputs
  # --------------------------
tet_order <- c(
  "GLCT", "RAKF",
  "RPPI", "FLRG", "QAKW",
  "AVFD*", "LLDF", "YLQQ",
  "CLGG", "FLYA",
  "RPQK*", "YNLR*"
)

# Use your already-filtered data (recommended)
# df_plot should already include your GEMEBV filtering etc.
df_in <- df_plot

out_dir <- "L:/Lab_CoreyS/Smith Lab - Manuscripts/EBV MS manuscript/Data/Single cell RNA seq data/EBV T cell baseline/circos_plots_beta/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

safe_filename <- function(x) {
  x %>%
    str_replace_all("\\*", "STAR") %>%
    str_replace_all("[<>:\"/\\\\|\\?\\*]", "_")
}

# --------------------------
# Core function for one tetramer
# --------------------------
make_pairwise_sharing_private <- function(df, tet_name, min_pair_cells_shared = 3) {
  
  df_tet <- df %>%
    filter(tetramer == tet_name) %>%
    filter(!is.na(sample), sample != "",
           !is.na(TRB_cdr3), TRB_cdr3 != "",
           !is.na(TRB_v_gene), TRB_v_gene != "",
           !is.na(TRB_j_gene), TRB_j_gene != "") %>%
    mutate(clone_key = paste(TRB_cdr3, TRB_v_gene, TRB_j_gene, sep = "|"))
  
  if (nrow(df_tet) == 0) return(NULL)
  
  # --- clone sets per sample (ALL clones; private included; no min-cells filter)
  clones_by_sample <- df_tet %>%
    distinct(sample, clone_key)
  
  total_clones <- clones_by_sample %>%
    count(sample, name = "total_clones")
  
  samples <- sort(unique(clones_by_sample$sample))
  if (length(samples) == 0) return(NULL)
  
  # --- cells per clone per sample (needed for the ≥3 cells across pair rule)
  clone_cells_sample <- df_tet %>%
    count(sample, clone_key, name = "n_cells")
  
  # named list of clone sets
  sets_named <- clones_by_sample %>%
    group_by(sample) %>%
    summarise(clones = list(unique(clone_key)), .groups = "drop") %>%
    { setNames(.$clones, .$sample) }
  
  # all directed pairs A->B (including diagonal; diagonal will be overwritten with private later)
  pairs <- expand_grid(sample_1 = samples, sample_2 = samples)
  
  # helper to compute eligible shared clones for a pair (A,B)
  pair_shared_n <- function(a, b) {
    if (a == b) return(NA_integer_)  # diagonal handled later
    
    shared <- intersect(sets_named[[a]], sets_named[[b]])
    if (length(shared) == 0) return(0L)
    
    # cells in A + cells in B for each shared clone
    a_counts <- clone_cells_sample %>%
      filter(sample == a, clone_key %in% shared) %>%
      select(clone_key, n_cells) %>%
      rename(n_a = n_cells)
    
    b_counts <- clone_cells_sample %>%
      filter(sample == b, clone_key %in% shared) %>%
      select(clone_key, n_cells) %>%
      rename(n_b = n_cells)
    
    pair_counts <- inner_join(a_counts, b_counts, by = "clone_key") %>%
      mutate(n_pair = n_a + n_b)
    
    sum(pair_counts$n_pair >= min_pair_cells_shared)
  }
  
  # compute eligible shared clone counts for each directed pair
  pair_df <- pairs %>%
    mutate(shared_clones_ge3 = map2_int(sample_1, sample_2, pair_shared_n)) %>%
    left_join(total_clones, by = c("sample_1" = "sample")) %>%
    rename(total_clones_sample1 = total_clones) %>%
    mutate(
      shared_prop = if_else(
        is.na(shared_clones_ge3) | is.na(total_clones_sample1) | total_clones_sample1 == 0,
        NA_real_,
        shared_clones_ge3 / total_clones_sample1
      ),
      shared_percent = shared_prop * 100,
      log10_shared_percent = if_else(shared_percent > 0, log10(shared_percent), NA_real_),
      tetramer = tet_name
    )
  
  # --- private proportion per sample_1 = 1 - sum(shared_prop to all OTHER samples)
  private_df <- pair_df %>%
    filter(sample_1 != sample_2) %>%
    group_by(sample_1) %>%
    summarise(
      private_prop = 1 - sum(shared_prop, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      private_prop = pmax(private_prop, 0),              # clamp to 0 if overlaps make it negative
      private_percent = private_prop * 100,
      log10_private_percent = if_else(private_percent > 0, log10(private_percent), NA_real_)
    )
  
  # replace diagonal rows with private values
  pair_df2 <- pair_df %>%
    left_join(private_df, by = c("sample_1")) %>%
    mutate(
      # for diagonal, use private; for off-diagonal keep shared
      final_prop = if_else(sample_1 == sample_2, private_prop, shared_prop),
      final_percent = final_prop * 100,
      final_log10_percent = if_else(final_percent > 0, log10(final_percent), NA_real_)
    )
  
  # add cohort columns if available
  if ("cohort" %in% names(df_tet)) {
    samp_cohort <- df_tet %>% distinct(sample, cohort)
    pair_df2 <- pair_df2 %>%
      left_join(samp_cohort, by = c("sample_1" = "sample")) %>%
      rename(cohort_1 = cohort) %>%
      left_join(samp_cohort, by = c("sample_2" = "sample")) %>%
      rename(cohort_2 = cohort)
  }
  
  list(
    full = pair_df2,
    three_col = pair_df2 %>% transmute(sample_1, sample_2, log10_value = final_log10_percent)
  )
}

# --------------------------
# Run per tetramer + save
# --------------------------
walk(tet_order, function(tet) {
  
  out <- make_pairwise_sharing_private(df_in, tet, min_pair_cells_shared = 3)
  
  if (is.null(out) || nrow(out$full) == 0) {
    message("Skipping ", tet, " (no data).")
    return(invisible(NULL))
  }
  
  safe_tet <- safe_filename(tet)
  
  write_csv(out$full,     file.path(out_dir, paste0(safe_tet, "_pairwise_sharing_PRIVATE_FULL.csv")))
  write_csv(out$three_col, file.path(out_dir, paste0(safe_tet, "_pairwise_sharing_PRIVATE_3col.csv")))
  
  message("Wrote: ", safe_tet)
})



# 2. Alpha chain

# ============================================================
# Pairwise clonal sharing (per tetramer) using TRA chains
# Shared clones between samples require >=3 CELLS across the pair
# Private clones are ALWAYS included in denominator (even singletons)
# Diagonal = private proportion = 1 - sum(shared_props to other samples)
# Output: one CSV per tetramer (FULL + 3-column version)
# ============================================================


# --------------------------
# Inputs
# --------------------------
tet_order <- c(
  "GLCT", "RAKF",
  "RPPI", "FLRG", "QAKW",
  "AVFD*", "LLDF", "YLQQ",
  "CLGG", "FLYA",
  "RPQK*", "YNLR*"
)

# Use your already-filtered data (recommended)
df_in <- df_plot

out_dir <- "L:/Lab_CoreyS/Smith Lab - Manuscripts/EBV MS manuscript/Data/Single cell RNA seq data/EBV T cell baseline/circos_plots_alpha/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

safe_filename <- function(x) {
  x %>%
    str_replace_all("\\*", "STAR") %>%
    str_replace_all("[<>:\"/\\\\|\\?\\*]", "_")
}

# --------------------------
# Core function for one tetramer (TRA-based)
# --------------------------
make_pairwise_sharing_private_TRA <- function(df, tet_name, min_pair_cells_shared = 3) {
  
  df_tet <- df %>%
    filter(tetramer == tet_name) %>%
    filter(!is.na(sample), sample != "",
           !is.na(TRA_cdr3), TRA_cdr3 != "",
           !is.na(TRA_v_gene), TRA_v_gene != "",
           !is.na(TRA_j_gene), TRA_j_gene != "") %>%
    mutate(clone_key = paste(TRA_cdr3, TRA_v_gene, TRA_j_gene, sep = "|"))
  
  if (nrow(df_tet) == 0) return(NULL)
  
  # --- clone sets per sample (ALL clones; private included; no min-cells filter)
  clones_by_sample <- df_tet %>%
    distinct(sample, clone_key)
  
  total_clones <- clones_by_sample %>%
    count(sample, name = "total_clones")
  
  samples <- sort(unique(clones_by_sample$sample))
  if (length(samples) == 0) return(NULL)
  
  # --- cells per clone per sample (needed for >=3 cells across pair rule)
  clone_cells_sample <- df_tet %>%
    count(sample, clone_key, name = "n_cells")
  
  # named list: sample -> vector of clone_keys
  sets_named <- clones_by_sample %>%
    group_by(sample) %>%
    summarise(clones = list(unique(clone_key)), .groups = "drop") %>%
    { setNames(.$clones, .$sample) }
  
  # all directed pairs A->B (including diagonal)
  pairs <- expand_grid(sample_1 = samples, sample_2 = samples)
  
  # eligible shared clones for a directed pair (A,B)
  pair_shared_n <- function(a, b) {
    if (a == b) return(NA_integer_)  # diagonal handled later
    
    shared <- intersect(sets_named[[a]], sets_named[[b]])
    if (length(shared) == 0) return(0L)
    
    a_counts <- clone_cells_sample %>%
      filter(sample == a, clone_key %in% shared) %>%
      select(clone_key, n_cells) %>%
      rename(n_a = n_cells)
    
    b_counts <- clone_cells_sample %>%
      filter(sample == b, clone_key %in% shared) %>%
      select(clone_key, n_cells) %>%
      rename(n_b = n_cells)
    
    pair_counts <- inner_join(a_counts, b_counts, by = "clone_key") %>%
      mutate(n_pair = n_a + n_b)
    
    sum(pair_counts$n_pair >= min_pair_cells_shared)
  }
  
  # compute eligible shared clone counts for each directed pair
  pair_df <- pairs %>%
    mutate(shared_clones_ge3 = map2_int(sample_1, sample_2, pair_shared_n)) %>%
    left_join(total_clones, by = c("sample_1" = "sample")) %>%
    rename(total_clones_sample1 = total_clones) %>%
    mutate(
      shared_prop = if_else(
        is.na(shared_clones_ge3) | is.na(total_clones_sample1) | total_clones_sample1 == 0,
        NA_real_,
        shared_clones_ge3 / total_clones_sample1
      ),
      shared_percent = shared_prop * 100,
      log10_shared_percent = if_else(shared_percent > 0, log10(shared_percent), NA_real_),
      tetramer = tet_name,
      chain = "TRA"
    )
  
  # private proportion per sample_1 = 1 - sum(shared_prop to all OTHER samples)
  private_df <- pair_df %>%
    filter(sample_1 != sample_2) %>%
    group_by(sample_1) %>%
    summarise(private_prop = 1 - sum(shared_prop, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      private_prop = pmax(private_prop, 0),  # clamp to 0 if overlaps make it negative
      private_percent = private_prop * 100,
      log10_private_percent = if_else(private_percent > 0, log10(private_percent), NA_real_)
    )
  
  # replace diagonal with private
  pair_df2 <- pair_df %>%
    left_join(private_df, by = "sample_1") %>%
    mutate(
      final_prop = if_else(sample_1 == sample_2, private_prop, shared_prop),
      final_percent = final_prop * 100,
      final_log10_percent = if_else(final_percent > 0, log10(final_percent), NA_real_)
    )
  
  # add cohort columns if present
  if ("cohort" %in% names(df_tet)) {
    samp_cohort <- df_tet %>% distinct(sample, cohort)
    pair_df2 <- pair_df2 %>%
      left_join(samp_cohort, by = c("sample_1" = "sample")) %>%
      rename(cohort_1 = cohort) %>%
      left_join(samp_cohort, by = c("sample_2" = "sample")) %>%
      rename(cohort_2 = cohort)
  }
  
  list(
    full = pair_df2,
    three_col = pair_df2 %>% transmute(sample_1, sample_2, log10_value = final_log10_percent)
  )
}

# --------------------------
# Run per tetramer + save
# --------------------------
walk(tet_order, function(tet) {
  
  out <- make_pairwise_sharing_private_TRA(df_in, tet, min_pair_cells_shared = 3)
  
  if (is.null(out) || nrow(out$full) == 0) {
    message("Skipping ", tet, " (no data).")
    return(invisible(NULL))
  }
  
  safe_tet <- safe_filename(tet)
  
  write_csv(out$full,      file.path(out_dir, paste0(safe_tet, "_TRA_pairwise_sharing_PRIVATE_FULL.csv")))
  write_csv(out$three_col, file.path(out_dir, paste0(safe_tet, "_TRA_pairwise_sharing_PRIVATE_3col.csv")))
  
  message("Wrote: ", safe_tet)
})





