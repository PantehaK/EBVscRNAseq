#!/usr/bin/env Rscript

# ==============================================================================
# Clone size, cytotoxicity score models, and clone-size UMAP overlays
# ==============================================================================
#
# Purpose:
#   This script:
#     1. Loads the final baseline EBV Seurat object.
#     2. Defines sample-level TRB clonotypes using:
#          TRB_cdr3 + TRB_v_gene + TRB_j_gene
#     3. Calculates clone size within each participant/sample ID.
#     4. Calculates log2(clone size + 1).
#     5. Tests the association between clone size and cytotoxicity score in
#        latent EBV-specific T cells using mixed models.
#     6. Fits cohort-specific models for MS and Control participants.
#     7. Saves cohort-specific clone size vs cytotoxicity plots.
#     8. Saves UMAP overlays for:
#          - latent MS
#          - latent Control
#          - lytic MS
#          - lytic Control
#          - negative MS
#          - negative Control
#
# Notes:
#   - Paths are intentionally generic for GitHub.
#   - Update the config block before running.
#   - Clone size is calculated within `id_col` using TRB-only clone identity.
#   - The mixed model uses `random_effect_col`, which defaults to `sample`.
#     Change this to `id` if that is the correct donor/sample identifier in your
#     object for random effects.
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
  library(readr)
  library(lme4)
  library(lmerTest)
  library(patchwork)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  
  # Input RDS
  input_rds = "path/to/Publication_RDS/Baseline_EBV_integrated_seurat.rds",
  
  # Output directory
  out_dir = "path/to/Publication_data/EBV_baseline/Clonal_expansion",
  
  # Seurat reduction for UMAP overlays
  reduction = "umap.harmony",
  
  # Metadata columns
  id_col = "id",
  random_effect_col = "sample",
  cohort_col = "cohort",
  lifecycle_col = "lifecycle",
  age_col = "age",
  cytotoxicity_col = "Cytotoxicity_score",
  
  # TRB columns used to define clone identity
  trb_cdr3_col = "TRB_cdr3",
  trb_v_col = "TRB_v_gene",
  trb_j_col = "TRB_j_gene",
  
  # Analysis subset
  latent_label = "Latent",
  lytic_label = "Lytic",
  negative_label = "negative",
  control_label = "Control",
  ms_label = "MS",
  
  # Plot settings
  ms_colour = "#e31a1c",
  control_colour = "#1f77b4",
  background_colour = "grey80",
  
  scatter_width = 10,
  scatter_height = 5,
  scatter_dpi = 600,
  
  umap_width = 7,
  umap_height = 6,
  umap_dpi = 600,
  
  combined_umap_width = 12,
  combined_umap_height = 16,
  combined_umap_dpi = 600,
  
  # Output files
  metadata_with_clone_size_csv = "baseline_metadata_with_TRB_clone_size.csv",
  latent_model_input_csv = "latent_clone_size_cytotoxicity_model_input.csv",
  model_summary_txt = "clone_size_cytotoxicity_lmer_summaries.txt",
  model_coefficients_csv = "clone_size_cytotoxicity_model_coefficients.csv",
  
  ms_scatter_png = "clone_size_vs_cytotoxicity_score_MS_Latent.png",
  control_scatter_png = "clone_size_vs_cytotoxicity_score_Control_Latent.png",
  
  latent_ms_umap_png = "umap_clone_size_latent_MS.png",
  latent_control_umap_png = "umap_clone_size_latent_Control.png",
  lytic_ms_umap_png = "umap_clone_size_lytic_MS.png",
  lytic_control_umap_png = "umap_clone_size_lytic_Control.png",
  negative_ms_umap_png = "umap_clone_size_negative_MS.png",
  negative_control_umap_png = "umap_clone_size_negative_Control.png",
  combined_umap_png = "umap_clone_size_lifecycle_cohort_combined.png",
  
  session_info_file = "sessionInfo_clone_size_cytotoxicity_umap.txt"
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
  str_trim(as.character(x))
}


clean_cohort <- function(x, control_label = "Control", ms_label = "MS") {
  
  x_raw <- as.character(x)
  x_clean <- str_to_lower(str_trim(x_raw))
  
  case_when(
    x_clean %in% c(
      "control", "controls", "healthy control", "healthy controls",
      "hc", "nms", "non-ms", "non_ms"
    ) ~ control_label,
    
    x_clean %in% c(
      "ms", "pwms", "multiple sclerosis", "rrms", "ppms", "spms"
    ) ~ ms_label,
    
    TRUE ~ x_raw
  )
}


extract_model_term <- function(model, model_name, term = "log2_clone_size") {
  
  coef_df <- as.data.frame(summary(model)$coefficients) %>%
    rownames_to_column("term")
  
  coef_df %>%
    filter(term == !!term) %>%
    mutate(model = model_name, .before = 1)
}


extract_all_model_coefficients <- function(model, model_name) {
  
  as.data.frame(summary(model)$coefficients) %>%
    rownames_to_column("term") %>%
    mutate(model = model_name, .before = 1)
}


make_prediction_df <- function(model, data, cfg) {
  
  x_seq <- seq(
    min(data$log2_clone_size, na.rm = TRUE),
    max(data$log2_clone_size, na.rm = TRUE),
    length.out = 100
  )
  
  pred_df <- tibble(
    log2_clone_size = x_seq,
    age_model = mean(data$age_model, na.rm = TRUE)
  )
  
  # lme4 predict requires the random-effect column to exist in newdata,
  # even when re.form = NA ignores random effects.
  pred_df[[cfg$random_effect_col]] <- factor(
    levels(data[[cfg$random_effect_col]])[1],
    levels = levels(data[[cfg$random_effect_col]])
  )
  
  pred_df$pred <- predict(
    model,
    newdata = pred_df,
    re.form = NA,
    allow.new.levels = TRUE
  )
  
  pred_df
}


make_scatter_plot <- function(data, pred_df, label_text, cohort_label, colour) {
  
  ggplot(
    data,
    aes(x = log2_clone_size, y = Cytotoxicity_score)
  ) +
    geom_point(
      colour = colour,
      alpha = 0.3,
      size = 1
    ) +
    geom_line(
      data = pred_df,
      aes(x = log2_clone_size, y = pred),
      colour = "black",
      linewidth = 1
    ) +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = label_text,
      hjust = -0.05,
      vjust = 1.1,
      size = 4,
      fontface = "bold"
    ) +
    labs(
      title = paste(cohort_label, "Latent EBV-specific T cells"),
      x = expression(log[2]~"(clone size + 1)"),
      y = "Cytotoxicity score"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text = element_text(size = 16),
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
}


plot_clone_overlay <- function(
    plot_df,
    lifecycle_to_plot,
    cohort_to_plot,
    title_text,
    high_colour,
    background_colour = "grey80"
) {
  
  overlay_df <- plot_df %>%
    filter(
      lifecycle == lifecycle_to_plot,
      cohort == cohort_to_plot
    )
  
  # Background = all cells except the exact lifecycle/cohort being plotted.
  bg_df <- plot_df %>%
    filter(
      !(lifecycle == lifecycle_to_plot & cohort == cohort_to_plot)
    )
  
  ggplot() +
    geom_point(
      data = bg_df,
      aes(umap_1, umap_2),
      colour = background_colour,
      size = 0.5,
      alpha = 0.6
    ) +
    geom_point(
      data = overlay_df,
      aes(
        umap_1,
        umap_2,
        fill = log2_clone_size,
        size = log2_clone_size
      ),
      shape = 21,
      colour = "black",
      stroke = 0.2,
      alpha = 0.9
    ) +
    scale_fill_gradient(
      low = "white",
      high = high_colour,
      na.value = "transparent",
      name = "log2(Clone size + 1)"
    ) +
    scale_size_continuous(
      range = c(0.1, 4),
      name = "log2(Clone size + 1)"
    ) +
    coord_equal() +
    theme_classic(base_size = 13) +
    labs(
      title = title_text,
      x = "UMAP 1",
      y = "UMAP 2"
    )
}


# ----------------------------- #
# 3. Load object and calculate TRB clone size
# ----------------------------- #

create_dir(config$out_dir)

message("Loading baseline EBV object:")
message(config$input_rds)

merged <- load_seurat_object(config$input_rds)

meta <- merged@meta.data %>%
  as.data.frame()

required_cols <- c(
  config$id_col,
  config$random_effect_col,
  config$cohort_col,
  config$lifecycle_col,
  config$age_col,
  config$cytotoxicity_col,
  config$trb_cdr3_col,
  config$trb_v_col,
  config$trb_j_col
)

check_required_columns(
  df = meta,
  required_cols = required_cols,
  object_name = "baseline EBV metadata"
)

meta <- meta %>%
  rownames_to_column("cell_barcode") %>%
  mutate(
    id = clean_chr(.data[[config$id_col]]),
    cohort = clean_cohort(
      .data[[config$cohort_col]],
      control_label = config$control_label,
      ms_label = config$ms_label
    ),
    lifecycle = clean_chr(.data[[config$lifecycle_col]]),
    age_model = suppressWarnings(as.numeric(.data[[config$age_col]])),
    Cytotoxicity_score = suppressWarnings(as.numeric(.data[[config$cytotoxicity_col]])),
    
    TRB_cdr3 = clean_chr(.data[[config$trb_cdr3_col]]),
    TRB_v_gene = clean_chr(.data[[config$trb_v_col]]),
    TRB_j_gene = clean_chr(.data[[config$trb_j_col]]),
    
    clone_id_TRB = if_else(
      !is.na(TRB_cdr3) & TRB_cdr3 != "" &
        !is.na(TRB_v_gene) & TRB_v_gene != "" &
        !is.na(TRB_j_gene) & TRB_j_gene != "",
      paste(TRB_cdr3, TRB_v_gene, TRB_j_gene, sep = "|"),
      NA_character_
    )
  )

clone_sizes <- meta %>%
  filter(
    !is.na(id), id != "",
    !is.na(clone_id_TRB), clone_id_TRB != ""
  ) %>%
  count(id, clone_id_TRB, name = "clone_size")

meta <- meta %>%
  left_join(
    clone_sizes,
    by = c("id", "clone_id_TRB")
  ) %>%
  mutate(
    clone_size = replace_na(clone_size, 1L),
    log2_clone_size = log2(clone_size + 1)
  ) %>%
  column_to_rownames("cell_barcode")

merged@meta.data <- meta

metadata_output_path <- file.path(config$out_dir, config$metadata_with_clone_size_csv)
write_csv(
  meta %>% rownames_to_column("cell_barcode"),
  metadata_output_path
)

message("Saved metadata with TRB clone size: ", metadata_output_path)


# ----------------------------- #
# 4. Prepare latent model input
# ----------------------------- #

df_latent <- merged@meta.data %>%
  as.data.frame() %>%
  filter(lifecycle == config$latent_label) %>%
  filter(
    !is.na(Cytotoxicity_score),
    !is.na(log2_clone_size),
    !is.na(age_model),
    !is.na(.data[[config$random_effect_col]]),
    !is.na(cohort),
    cohort %in% c(config$control_label, config$ms_label)
  ) %>%
  mutate(
    cohort = factor(
      cohort,
      levels = c(config$control_label, config$ms_label)
    ),
    across(all_of(config$random_effect_col), as.factor)
  )

if (nrow(df_latent) == 0) {
  stop("No latent cells remained for clone size/cytotoxicity modelling.")
}

if (length(unique(df_latent$cohort)) < 2) {
  stop("Both Control and MS cohorts are required for the interaction model.")
}

latent_input_path <- file.path(config$out_dir, config$latent_model_input_csv)
write_csv(df_latent, latent_input_path)

message("Saved latent model input: ", latent_input_path)


# ----------------------------- #
# 5. Fit mixed models
# ----------------------------- #

message("Fitting interaction model...")

interaction_formula <- as.formula(
  paste0(
    "Cytotoxicity_score ~ log2_clone_size * cohort + age_model + (1 | ",
    config$random_effect_col,
    ")"
  )
)

m_int <- lmer(
  interaction_formula,
  data = df_latent
)

message("Fitting cohort-specific models...")

cohort_formula <- as.formula(
  paste0(
    "Cytotoxicity_score ~ log2_clone_size + age_model + (1 | ",
    config$random_effect_col,
    ")"
  )
)

df_ctl <- df_latent %>%
  filter(cohort == config$control_label)

df_ms <- df_latent %>%
  filter(cohort == config$ms_label)

if (nrow(df_ctl) == 0 || nrow(df_ms) == 0) {
  stop("Cohort-specific models require both Control and MS latent cells.")
}

m_ctl <- lmer(
  cohort_formula,
  data = df_ctl
)

m_ms <- lmer(
  cohort_formula,
  data = df_ms
)


# ----------------------------- #
# 6. Save model summaries and coefficients
# ----------------------------- #

model_summary_path <- file.path(config$out_dir, config$model_summary_txt)

writeLines(
  c(
    "Interaction model:",
    capture.output(summary(m_int)),
    "",
    "Control-only model:",
    capture.output(summary(m_ctl)),
    "",
    "MS-only model:",
    capture.output(summary(m_ms))
  ),
  model_summary_path
)

model_coefficients <- bind_rows(
  extract_all_model_coefficients(m_int, "interaction_model"),
  extract_all_model_coefficients(m_ctl, "control_model"),
  extract_all_model_coefficients(m_ms, "ms_model")
)

model_coefficients_path <- file.path(config$out_dir, config$model_coefficients_csv)
write_csv(model_coefficients, model_coefficients_path)

message("Saved model summaries: ", model_summary_path)
message("Saved model coefficients: ", model_coefficients_path)


# ----------------------------- #
# 7. Generate clone size/cytotoxicity plots
# ----------------------------- #

coef_ctl <- extract_model_term(m_ctl, "control_model", "log2_clone_size")
coef_ms <- extract_model_term(m_ms, "ms_model", "log2_clone_size")

beta_ctl <- coef_ctl$Estimate[1]
p_ctl <- coef_ctl$`Pr(>|t|)`[1]

beta_ms <- coef_ms$Estimate[1]
p_ms <- coef_ms$`Pr(>|t|)`[1]

label_ctl <- paste0(
  "β = ",
  round(beta_ctl, 3),
  "\nP ",
  format.pval(p_ctl, digits = 2, eps = 1e-300)
)

label_ms <- paste0(
  "β = ",
  round(beta_ms, 3),
  "\nP ",
  format.pval(p_ms, digits = 2, eps = 1e-300)
)

new_ctl <- make_prediction_df(m_ctl, df_ctl, config)
new_ms <- make_prediction_df(m_ms, df_ms, config)

p_ctl <- make_scatter_plot(
  data = df_ctl,
  pred_df = new_ctl,
  label_text = label_ctl,
  cohort_label = config$control_label,
  colour = config$control_colour
)

p_ms <- make_scatter_plot(
  data = df_ms,
  pred_df = new_ms,
  label_text = label_ms,
  cohort_label = config$ms_label,
  colour = config$ms_colour
)

ms_scatter_path <- file.path(config$out_dir, config$ms_scatter_png)
control_scatter_path <- file.path(config$out_dir, config$control_scatter_png)

ggsave(
  filename = ms_scatter_path,
  plot = p_ms,
  bg = "transparent",
  width = config$scatter_width,
  height = config$scatter_height,
  dpi = config$scatter_dpi
)

ggsave(
  filename = control_scatter_path,
  plot = p_ctl,
  bg = "transparent",
  width = config$scatter_width,
  height = config$scatter_height,
  dpi = config$scatter_dpi
)

message("Saved MS scatter plot: ", ms_scatter_path)
message("Saved Control scatter plot: ", control_scatter_path)


# ----------------------------- #
# 8. Generate clone-size UMAP overlays
# ----------------------------- #

if (!config$reduction %in% Reductions(merged)) {
  stop(
    "Reduction not found: ",
    config$reduction,
    ". Available reductions: ",
    paste(Reductions(merged), collapse = ", ")
  )
}

umap <- as.data.frame(Embeddings(merged, reduction = config$reduction)) %>%
  rownames_to_column("cell_barcode")

colnames(umap)[2:3] <- c("umap_1", "umap_2")

plot_df <- merged@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell_barcode") %>%
  left_join(umap, by = "cell_barcode") %>%
  mutate(
    lifecycle = clean_chr(lifecycle),
    cohort = clean_cohort(
      cohort,
      control_label = config$control_label,
      ms_label = config$ms_label
    ),
    log2_clone_size = as.numeric(log2_clone_size)
  )

required_plot_cols <- c(
  "lifecycle",
  "cohort",
  "log2_clone_size",
  "umap_1",
  "umap_2"
)

check_required_columns(
  df = plot_df,
  required_cols = required_plot_cols,
  object_name = "UMAP plotting dataframe"
)

p_latent_MS <- plot_clone_overlay(
  plot_df = plot_df,
  lifecycle_to_plot = config$latent_label,
  cohort_to_plot = config$ms_label,
  title_text = "MS latent cells coloured by clonal expansion",
  high_colour = config$ms_colour,
  background_colour = config$background_colour
)

p_latent_Control <- plot_clone_overlay(
  plot_df = plot_df,
  lifecycle_to_plot = config$latent_label,
  cohort_to_plot = config$control_label,
  title_text = "Control latent cells coloured by clonal expansion",
  high_colour = config$control_colour,
  background_colour = config$background_colour
)

p_lytic_MS <- plot_clone_overlay(
  plot_df = plot_df,
  lifecycle_to_plot = config$lytic_label,
  cohort_to_plot = config$ms_label,
  title_text = "MS lytic cells coloured by clonal expansion",
  high_colour = config$ms_colour,
  background_colour = config$background_colour
)

p_lytic_Control <- plot_clone_overlay(
  plot_df = plot_df,
  lifecycle_to_plot = config$lytic_label,
  cohort_to_plot = config$control_label,
  title_text = "Control lytic cells coloured by clonal expansion",
  high_colour = config$control_colour,
  background_colour = config$background_colour
)

p_negative_MS <- plot_clone_overlay(
  plot_df = plot_df,
  lifecycle_to_plot = config$negative_label,
  cohort_to_plot = config$ms_label,
  title_text = "MS negative cells coloured by clonal expansion",
  high_colour = config$ms_colour,
  background_colour = config$background_colour
)

p_negative_Control <- plot_clone_overlay(
  plot_df = plot_df,
  lifecycle_to_plot = config$negative_label,
  cohort_to_plot = config$control_label,
  title_text = "Control negative cells coloured by clonal expansion",
  high_colour = config$control_colour,
  background_colour = config$background_colour
)

umap_paths <- list(
  latent_ms = file.path(config$out_dir, config$latent_ms_umap_png),
  latent_control = file.path(config$out_dir, config$latent_control_umap_png),
  lytic_ms = file.path(config$out_dir, config$lytic_ms_umap_png),
  lytic_control = file.path(config$out_dir, config$lytic_control_umap_png),
  negative_ms = file.path(config$out_dir, config$negative_ms_umap_png),
  negative_control = file.path(config$out_dir, config$negative_control_umap_png)
)

ggsave(umap_paths$latent_ms, p_latent_MS, width = config$umap_width, height = config$umap_height, dpi = config$umap_dpi, bg = "transparent")
ggsave(umap_paths$latent_control, p_latent_Control, width = config$umap_width, height = config$umap_height, dpi = config$umap_dpi, bg = "transparent")
ggsave(umap_paths$lytic_ms, p_lytic_MS, width = config$umap_width, height = config$umap_height, dpi = config$umap_dpi, bg = "transparent")
ggsave(umap_paths$lytic_control, p_lytic_Control, width = config$umap_width, height = config$umap_height, dpi = config$umap_dpi, bg = "transparent")
ggsave(umap_paths$negative_ms, p_negative_MS, width = config$umap_width, height = config$umap_height, dpi = config$umap_dpi, bg = "transparent")
ggsave(umap_paths$negative_control, p_negative_Control, width = config$umap_width, height = config$umap_height, dpi = config$umap_dpi, bg = "transparent")

combined_clone_umap <- (
  p_latent_MS + p_latent_Control +
    p_lytic_MS + p_lytic_Control +
    p_negative_MS + p_negative_Control
) +
  plot_layout(ncol = 2)

combined_umap_path <- file.path(config$out_dir, config$combined_umap_png)

ggsave(
  filename = combined_umap_path,
  plot = combined_clone_umap,
  width = config$combined_umap_width,
  height = config$combined_umap_height,
  dpi = config$combined_umap_dpi,
  bg = "transparent"
)


# ----------------------------- #
# 9. Save session information
# ----------------------------- #

session_info_path <- file.path(config$out_dir, config$session_info_file)

writeLines(
  c(
    "Configuration:",
    capture.output(str(config)),
    "",
    "Clone size definition:",
    "TRB_cdr3 + TRB_v_gene + TRB_j_gene, calculated within id.",
    "",
    "Latent model input rows:",
    as.character(nrow(df_latent)),
    "",
    "Output files:",
    metadata_output_path,
    latent_input_path,
    model_summary_path,
    model_coefficients_path,
    ms_scatter_path,
    control_scatter_path,
    unlist(umap_paths),
    combined_umap_path,
    "",
    "Session information:",
    capture.output(sessionInfo())
  ),
  session_info_path
)

message("\nClone size/cytotoxicity/UMAP workflow complete.")
message("Session info: ", session_info_path)
