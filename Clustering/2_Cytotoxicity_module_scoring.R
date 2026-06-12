#!/usr/bin/env Rscript

# ==============================================================================
# Cytotoxicity module scoring and clone-size association analysis
# ==============================================================================
#
# Purpose:
#   This script adds a literature-supported cytotoxicity module score to a
#   clustered EBV baseline CD8+ T cell Seurat object, visualises the score on UMAP,
#   calculates clone sizes, and tests whether cytotoxicity score is associated
#   with clonal expansion in latent EBV-specific T cells.
#
#   It saves:
#     1. a module-scored Seurat object,
#     2. cytotoxicity UMAP,
#     3. clone-size versus cytotoxicity plots for MS and Control,
#     4. mixed-model summaries,
#     5. model coefficient table,
#     6. session information.
#
# Expected input:
#   A clustered baseline EBV CD8+ T cell Seurat object, e.g.
#     14_baseline_EBV_clustered.rds
#
# Notes:
#   - Clone size is calculated within sample/id using TRA_cdr3 + TRB_cdr3.
#   - The mixed model is restricted to latent EBV-specific cells.
#   - Models adjust for age and include sample as a random intercept.
#
# ==============================================================================


# ----------------------------- #
# 0. Load packages
# ----------------------------- #

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(readr)
  library(stringr)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  dataset_id = "baseline_ebv",
  
  input_rds = "path/to/input/14_baseline_EBV_clustered.rds",
  output_rds = "path/to/output/15_baseline_EBV_module_scored.rds",
  
  plot_dir = "path/to/output/Publication_data/EBV_baseline/Module_scores",
  model_dir = "path/to/output/Publication_data/EBV_baseline/Module_scores/models",
  
  cytotoxicity_umap_png = "path/to/output/Publication_data/EBV_baseline/Module_scores/umap_cytotoxicity.png",
  clone_cyto_ms_png = "path/to/output/Publication_data/EBV_baseline/Module_scores/clone_size_vs_cytotoxicity_score_MS_Latent.png",
  clone_cyto_control_png = "path/to/output/Publication_data/EBV_baseline/Module_scores/clone_size_vs_cytotoxicity_score_Control_Latent.png",
  
  model_coefficients_csv = "path/to/output/Publication_data/EBV_baseline/Module_scores/clone_size_cytotoxicity_model_coefficients.csv",
  latent_model_data_csv = "path/to/output/Publication_data/EBV_baseline/Module_scores/latent_clone_size_cytotoxicity_model_data.csv",
  model_summary_txt = "path/to/output/Publication_data/EBV_baseline/Module_scores/models/lmer_model_summaries.txt",
  session_info_file = "path/to/output/Publication_data/EBV_baseline/Module_scores/sessionInfo_cytotoxicity_module_score.txt",
  
  assay = "SCT",
  reduction = "umap.harmony",
  
  score_name = "Cytotoxicity",
  score_column = "Cytotoxicity_score",
  
  cytotoxicity_genes = c(
    # Granzymes
    "GZMB", "GZMH", "GZMK", "GZMA", "GZMM",
    
    # Cytotoxic machinery
    "PRF1", "GNLY", "NKG7", "FGFBP2", "CTSW",
    "KLRD1", "KLRB1", "KLRK1",
    
    # Effector transcription
    "TBX21", "EOMES", "ZEB2", "ID2",
    
    # Interferon / effector support
    "IFNG", "IFITM1", "IFITM2", "IFITM3",
    
    # Terminal differentiation
    "CX3CR1", "S100A4", "S100A10"
  ),
  
  # Clone-size settings.
  clone_group_col = "id",
  tra_col = "TRA_cdr3",
  trb_col = "TRB_cdr3",
  
  # Model settings.
  lifecycle_col = "lifecycle",
  latent_label = "Latent",
  cohort_col = "cohort",
  sample_col = "sample",
  age_col = "age",
  cohort_levels = c("Control", "MS"),
  
  # Plot colours.
  cytotoxicity_umap_cols = c(
    "#E6F2FF", "#B3D9FF", "#7FB3FF",
    "#9B6EDB", "#D65DB1", "#F768A1",
    "#F89540", "#D7301F"
  ),
  ms_colour = "#e31a1c",
  control_colour = "#1f77b4",
  
  verbose = TRUE
)


# ----------------------------- #
# 2. Helper functions
# ----------------------------- #

create_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}


create_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}


setup_dirs <- function(cfg) {
  create_dir(cfg$plot_dir)
  create_dir(cfg$model_dir)
  
  create_parent_dir(cfg$output_rds)
  create_parent_dir(cfg$cytotoxicity_umap_png)
  create_parent_dir(cfg$clone_cyto_ms_png)
  create_parent_dir(cfg$clone_cyto_control_png)
  create_parent_dir(cfg$model_coefficients_csv)
  create_parent_dir(cfg$latent_model_data_csv)
  create_parent_dir(cfg$model_summary_txt)
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


remove_existing_module_score_cols <- function(obj, cfg) {
  cols_to_remove <- grep(
    paste0("^", cfg$score_name),
    colnames(obj@meta.data),
    value = TRUE
  )
  
  cols_to_remove <- union(cols_to_remove, cfg$score_column)
  cols_to_remove <- intersect(cols_to_remove, colnames(obj@meta.data))
  
  if (length(cols_to_remove) > 0) {
    obj@meta.data <- obj@meta.data |>
      select(-all_of(cols_to_remove))
  }
  
  obj
}


add_cytotoxicity_score <- function(obj, cfg) {
  if (!cfg$assay %in% Assays(obj)) {
    stop("Assay not found in object: ", cfg$assay)
  }
  
  DefaultAssay(obj) <- cfg$assay
  
  genes_present <- intersect(cfg$cytotoxicity_genes, rownames(obj))
  genes_missing <- setdiff(cfg$cytotoxicity_genes, genes_present)
  
  if (length(genes_present) == 0) {
    stop("None of the cytotoxicity genes were found in the object.")
  }
  
  if (length(genes_missing) > 0) {
    message(
      "Missing cytotoxicity genes: ",
      paste(genes_missing, collapse = ", ")
    )
  }
  
  obj <- remove_existing_module_score_cols(obj, cfg)
  
  obj <- AddModuleScore(
    object = obj,
    assay = cfg$assay,
    features = list(genes_present),
    name = cfg$score_name
  )
  
  score_cols <- grep(
    paste0("^", cfg$score_name),
    colnames(obj@meta.data),
    value = TRUE
  )
  
  newest_score_col <- tail(score_cols, 1)
  
  obj[[cfg$score_column]] <- obj[[newest_score_col, drop = TRUE]]
  
  # Remove the temporary Seurat-generated score column if it is not the final name.
  if (newest_score_col != cfg$score_column) {
    obj@meta.data <- obj@meta.data |>
      select(-all_of(newest_score_col))
  }
  
  obj
}


save_cytotoxicity_umap <- function(obj, cfg) {
  if (!cfg$score_column %in% colnames(obj@meta.data)) {
    stop("Score column not found: ", cfg$score_column)
  }
  
  score_range <- range(obj@meta.data[[cfg$score_column]], na.rm = TRUE)
  
  p <- FeaturePlot(
    object = obj,
    features = cfg$score_column,
    reduction = cfg$reduction,
    order = FALSE,
    min.cutoff = NA,
    max.cutoff = NA
  ) +
    scale_color_gradientn(
      colours = cfg$cytotoxicity_umap_cols,
      limits = score_range,
      oob = scales::squish
    )
  
  ggsave(
    plot = p,
    filename = cfg$cytotoxicity_umap_png,
    height = 8,
    width = 9,
    dpi = 300,
    bg = "transparent"
  )
  
  invisible(p)
}


add_clone_size_metadata <- function(obj, cfg) {
  required_cols <- c(cfg$tra_col, cfg$trb_col, cfg$clone_group_col)
  missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
  
  if (length(missing_cols) > 0) {
    stop(
      "Cannot calculate clone size. Missing metadata column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  meta <- obj@meta.data |>
    mutate(
      clone_id = paste(
        .data[[cfg$tra_col]],
        .data[[cfg$trb_col]],
        sep = "_"
      )
    ) |>
    group_by(.data[[cfg$clone_group_col]], clone_id) |>
    mutate(
      clone_size = n()
    ) |>
    ungroup() |>
    mutate(
      log2_clone_size = log2(clone_size + 1)
    )
  
  obj@meta.data <- meta
  
  obj
}


make_latent_model_dataframe <- function(obj, cfg) {
  required_cols <- c(
    cfg$lifecycle_col,
    cfg$score_column,
    "log2_clone_size",
    cfg$age_col,
    cfg$sample_col,
    cfg$cohort_col
  )
  
  missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
  
  if (length(missing_cols) > 0) {
    stop(
      "Cannot build model dataframe. Missing metadata column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  obj@meta.data |>
    filter(.data[[cfg$lifecycle_col]] == cfg$latent_label) |>
    filter(
      !is.na(.data[[cfg$score_column]]),
      !is.na(log2_clone_size),
      !is.na(.data[[cfg$age_col]]),
      !is.na(.data[[cfg$sample_col]]),
      !is.na(.data[[cfg$cohort_col]])
    ) |>
    mutate(
      sample = factor(.data[[cfg$sample_col]]),
      cohort = factor(.data[[cfg$cohort_col]], levels = cfg$cohort_levels),
      age = as.numeric(.data[[cfg$age_col]]),
      Cytotoxicity_score = as.numeric(.data[[cfg$score_column]])
    )
}


safe_lmer <- function(formula, data) {
  tryCatch(
    lmer(formula, data = data),
    error = function(e) {
      warning("Model failed: ", e$message)
      NULL
    }
  )
}


extract_lmer_term <- function(model, term, model_name, cohort = NA_character_) {
  if (is.null(model)) {
    return(tibble(
      model = model_name,
      cohort = cohort,
      term = term,
      estimate = NA_real_,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_
    ))
  }
  
  coefs <- summary(model)$coefficients
  
  if (!term %in% rownames(coefs)) {
    return(tibble(
      model = model_name,
      cohort = cohort,
      term = term,
      estimate = NA_real_,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_
    ))
  }
  
  tibble(
    model = model_name,
    cohort = cohort,
    term = term,
    estimate = coefs[term, "Estimate"],
    std_error = coefs[term, "Std. Error"],
    statistic = coefs[term, "t value"],
    p_value = coefs[term, "Pr(>|t|)"]
  )
}


fit_clone_size_models <- function(df_latent, cfg) {
  if (nrow(df_latent) == 0) {
    stop("No latent cells available for model fitting.")
  }
  
  interaction_model <- safe_lmer(
    Cytotoxicity_score ~ log2_clone_size * cohort + age + (1 | sample),
    data = df_latent
  )
  
  cohort_models <- list()
  
  for (cohort_name in cfg$cohort_levels) {
    df_sub <- df_latent |>
      filter(cohort == cohort_name)
    
    if (nrow(df_sub) == 0 || n_distinct(df_sub$sample) < 2) {
      warning(
        "Skipping cohort-specific model for ",
        cohort_name,
        ": not enough cells/samples."
      )
      cohort_models[[cohort_name]] <- NULL
      next
    }
    
    cohort_models[[cohort_name]] <- safe_lmer(
      Cytotoxicity_score ~ log2_clone_size + age + (1 | sample),
      data = df_sub
    )
  }
  
  coefficient_table <- bind_rows(
    extract_lmer_term(
      interaction_model,
      "log2_clone_size",
      model_name = "interaction_model",
      cohort = "reference"
    ),
    bind_rows(lapply(names(cohort_models), function(cohort_name) {
      extract_lmer_term(
        cohort_models[[cohort_name]],
        "log2_clone_size",
        model_name = "cohort_specific_model",
        cohort = cohort_name
      )
    }))
  )
  
  list(
    interaction_model = interaction_model,
    cohort_models = cohort_models,
    coefficient_table = coefficient_table
  )
}


write_model_summaries <- function(model_results, cfg) {
  output <- c()
  
  output <- c(output, "Interaction model:")
  output <- c(output, "==================")
  output <- c(
    output,
    if (is.null(model_results$interaction_model)) {
      "Model failed."
    } else {
      capture.output(summary(model_results$interaction_model))
    }
  )
  
  for (cohort_name in names(model_results$cohort_models)) {
    output <- c(output, "")
    output <- c(output, paste0("Cohort-specific model: ", cohort_name))
    output <- c(output, "========================================")
    
    model <- model_results$cohort_models[[cohort_name]]
    
    output <- c(
      output,
      if (is.null(model)) {
        "Model failed or was skipped."
      } else {
        capture.output(summary(model))
      }
    )
  }
  
  writeLines(output, cfg$model_summary_txt)
}


make_prediction_dataframe <- function(model, df_latent, cohort_name) {
  if (is.null(model)) {
    return(tibble())
  }
  
  df_sub <- df_latent |>
    filter(cohort == cohort_name)
  
  x_seq <- seq(
    min(df_sub$log2_clone_size, na.rm = TRUE),
    max(df_sub$log2_clone_size, na.rm = TRUE),
    length.out = 100
  )
  
  new_df <- data.frame(
    log2_clone_size = x_seq,
    age = mean(df_sub$age, na.rm = TRUE),
    sample = factor(levels(df_latent$sample)[1], levels = levels(df_latent$sample))
  )
  
  new_df$pred <- predict(model, newdata = new_df, re.form = NA)
  
  as_tibble(new_df)
}


make_clone_cytotoxicity_plot <- function(
    df_latent,
    model,
    coefficient_table,
    cohort_name,
    point_colour,
    output_file
) {
  df_plot <- df_latent |>
    filter(cohort == cohort_name)
  
  if (nrow(df_plot) == 0 || is.null(model)) {
    warning("Skipping plot for ", cohort_name, ": missing data or model.")
    return(invisible(NULL))
  }
  
  pred_df <- make_prediction_dataframe(
    model = model,
    df_latent = df_latent,
    cohort_name = cohort_name
  )
  
  coef_row <- coefficient_table |>
    filter(
      model == "cohort_specific_model",
      cohort == cohort_name,
      term == "log2_clone_size"
    )
  
  beta <- coef_row$estimate[[1]]
  p_val <- coef_row$p_value[[1]]
  
  label_text <- paste0(
    "\u03b2 = ", round(beta, 3),
    "\nP ", format.pval(p_val, digits = 2, eps = 1e-300)
  )
  
  p <- ggplot(
    df_plot,
    aes(x = log2_clone_size, y = Cytotoxicity_score)
  ) +
    geom_point(
      colour = point_colour,
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
      title = paste(cohort_name, "Latent EBV-specific T cells"),
      x = expression(log[2]~"(clone size + 1)"),
      y = "Cytotoxicity score"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text = element_text(size = 16),
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
  
  ggsave(
    plot = p,
    filename = output_file,
    bg = "transparent",
    width = 10,
    height = 5,
    dpi = 600
  )
  
  invisible(p)
}


# ----------------------------- #
# 3. Run pipeline
# ----------------------------- #

setup_dirs(config)

message("Loading clustered baseline object...")
merged <- load_seurat_object(config$input_rds)

message("Adding cytotoxicity module score...")
merged <- add_cytotoxicity_score(merged, config)

message("Saving module-scored object...")
saveRDS(merged, config$output_rds)

message("Saving cytotoxicity UMAP...")
save_cytotoxicity_umap(merged, config)

message("Calculating clone sizes...")
merged <- add_clone_size_metadata(merged, config)

message("Preparing latent-cell model dataframe...")
df_latent <- make_latent_model_dataframe(merged, config)

write.csv(
  df_latent,
  config$latent_model_data_csv,
  row.names = FALSE
)

message("Fitting mixed models...")
model_results <- fit_clone_size_models(df_latent, config)

write.csv(
  model_results$coefficient_table,
  config$model_coefficients_csv,
  row.names = FALSE
)

write_model_summaries(model_results, config)

message("Saving clone-size association plots...")
make_clone_cytotoxicity_plot(
  df_latent = df_latent,
  model = model_results$cohort_models[["MS"]],
  coefficient_table = model_results$coefficient_table,
  cohort_name = "MS",
  point_colour = config$ms_colour,
  output_file = config$clone_cyto_ms_png
)

make_clone_cytotoxicity_plot(
  df_latent = df_latent,
  model = model_results$cohort_models[["Control"]],
  coefficient_table = model_results$coefficient_table,
  cohort_name = "Control",
  point_colour = config$control_colour,
  output_file = config$clone_cyto_control_png
)

message("Saving final object with clone-size metadata...")
saveRDS(merged, config$output_rds)

writeLines(
  capture.output(sessionInfo()),
  config$session_info_file
)

message("\nCytotoxicity module scoring and clone-size analysis complete.")
message("Saved module-scored object to: ", config$output_rds)
message("Saved cytotoxicity UMAP to: ", config$cytotoxicity_umap_png)
message("Saved model coefficients to: ", config$model_coefficients_csv)
message("Saved model summaries to: ", config$model_summary_txt)
message("Saved session info to: ", config$session_info_file)