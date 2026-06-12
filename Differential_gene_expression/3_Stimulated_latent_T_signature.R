#!/usr/bin/env Rscript

# ==============================================================================
# Stimulated versus baseline latent EBV-specific T cell DE analysis
# ==============================================================================
#
# Purpose:
#   This script compares stimulated latent EBV-specific T cells against baseline
#   latent EBV-specific T cells to identify a stimulation-associated latent
#   T cell transcriptional signature.
#
#   It performs:
#     1. latent group annotation,
#     2. Seurat FindMarkers negative-binomial DE with age as covariate,
#     3. optional glmmTMB negative-binomial GLMM validation with sample random effect,
#     4. volcano plot,
#     5. labelled volcano plot,
#     6. colour-coded waterfall/BioRender table,
#     7. summary tables and session information.
#
# Expected input:
#   Activated paired EBV CD8+ T cell object, for example:
#     15_activated_EBV_MS_stimLat_signature_scored.rds
#
# Notes:
#   - FindMarkers compares:
#       Stimulated Latent vs Baseline Latent
#   - GLMM validation uses raw RNA counts and includes:
#       latent_group + age + offset(log(library size)) + (1 | sample)
#   - log2FC is Stimulated Latent relative to Baseline Latent.
#
# ==============================================================================


# ----------------------------- #
# 0. Load packages
# ----------------------------- #

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(stringr)
  library(glmmTMB)
  library(Matrix)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  dataset_id = "activated_ebv",
  
  input_rds = "path/to/input/15_activated_EBV_MS_stimLat_signature_scored.rds",
  
  output_dir = "path/to/output/Publication_data/EBV_activated/DE_stimulated_vs_baseline_latent",
  
  findmarkers_csv = "path/to/output/Publication_data/EBV_activated/DE_stimulated_vs_baseline_latent/NegBinom_Stimulated_vs_Baseline_agecov.csv",
  glmm_validation_csv = "path/to/output/Publication_data/EBV_activated/DE_stimulated_vs_baseline_latent/GLMM_glmmTMB_Stimulated_vs_Baseline_agecov_randomSample.csv",
  volcano_png = "path/to/output/Publication_data/EBV_activated/DE_stimulated_vs_baseline_latent/volcano_plot_stim_genes.png",
  volcano_labelled_png = "path/to/output/Publication_data/EBV_activated/DE_stimulated_vs_baseline_latent/volcano_plot_stim_genes_labelled.png",
  waterfall_colour_csv = "path/to/output/Publication_data/EBV_activated/DE_stimulated_vs_baseline_latent/waterfall_plot_stim_genes_colour_coded.csv",
  de_summary_csv = "path/to/output/Publication_data/EBV_activated/DE_stimulated_vs_baseline_latent/DE_summary.csv",
  session_info_file = "path/to/output/Publication_data/EBV_activated/DE_stimulated_vs_baseline_latent/sessionInfo_stimulated_vs_baseline_latent_DE.txt",
  
  # Metadata columns.
  batch_col = "batch",
  lifecycle_col = "lifecycle",
  age_col = "age",
  donor_col = "sample",
  
  # Labels.
  latent_label = "Latent",
  non_latent_label = "Non-Latent",
  baseline_latent_label = "Baseline Latent",
  stimulated_latent_label = "Stimulated Latent",
  
  baseline_batch_pattern = "^GEMEBV",
  stimulated_batch_pattern = "^EBVLCL",
  
  # Assay settings.
  de_assay = "RNA",
  raw_counts_assay = "RNA",
  
  # FindMarkers settings.
  test_use = "negbinom",
  latent_vars = c("age"),
  min_pct = 0.25,
  logfc_threshold = 0.0,
  
  # Significance thresholds.
  padj_threshold = 0.05,
  logfc_threshold_for_sig = 2,
  
  # GLMM validation.
  run_glmm_validation = TRUE,
  
  # Use "all_sig" to validate all significant genes passing GLMM filters.
  # Or set an integer, e.g. 100.
  validate_top_n_genes = "all_sig",
  
  glmm_gene_padj_threshold = 0.05,
  glmm_gene_abs_logfc_threshold = 2,
  
  # Volcano labels.
  n_label_up = 15,
  n_label_down = 15,
  
  # Waterfall output.
  n_waterfall_up = 30,
  n_waterfall_down = 30,
  
  blue_palette = c(
    "#011f4b", "#03396c", "#005b96", "#6497b1", "#b3cde0"
  ),
  red_palette = c(
    "#ffefea", "#fbd9d3", "#ffb09c", "#fe5757", "#cb2424", "#900000"
  ),
  
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
  create_parent_dir(cfg$findmarkers_csv)
  create_parent_dir(cfg$glmm_validation_csv)
  create_parent_dir(cfg$volcano_png)
  create_parent_dir(cfg$volcano_labelled_png)
  create_parent_dir(cfg$waterfall_colour_csv)
  create_parent_dir(cfg$de_summary_csv)
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


get_assay_data_compat <- function(seurat_obj, assay, layer_or_slot = "counts") {
  tryCatch(
    GetAssayData(seurat_obj, assay = assay, layer = layer_or_slot),
    error = function(e) {
      GetAssayData(seurat_obj, assay = assay, slot = layer_or_slot)
    }
  )
}


add_latent_group <- function(obj, cfg) {
  required_cols <- c(cfg$batch_col, cfg$lifecycle_col)
  missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
  
  if (length(missing_cols) > 0) {
    stop(
      "Object metadata is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  md <- obj@meta.data
  
  latent_group <- case_when(
    !is.na(md[[cfg$batch_col]]) &
      grepl(cfg$baseline_batch_pattern, md[[cfg$batch_col]]) &
      md[[cfg$lifecycle_col]] == cfg$latent_label ~ cfg$baseline_latent_label,
    
    !is.na(md[[cfg$batch_col]]) &
      grepl(cfg$stimulated_batch_pattern, md[[cfg$batch_col]]) &
      md[[cfg$lifecycle_col]] == cfg$latent_label ~ cfg$stimulated_latent_label,
    
    TRUE ~ cfg$non_latent_label
  )
  
  obj$latent_group <- factor(
    latent_group,
    levels = c(
      cfg$non_latent_label,
      cfg$baseline_latent_label,
      cfg$stimulated_latent_label
    )
  )
  
  obj
}


prepare_de_object <- function(obj, cfg) {
  required_cols <- c("latent_group", cfg$age_col, cfg$donor_col)
  missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
  
  if (length(missing_cols) > 0) {
    stop(
      "Object metadata is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  md <- obj@meta.data
  
  keep_cells <- rownames(md)[
    md$latent_group %in% c(cfg$baseline_latent_label, cfg$stimulated_latent_label) &
      !is.na(md[[cfg$age_col]])
  ]
  
  if (length(keep_cells) == 0) {
    stop("No baseline/stimulated latent cells with non-missing age were found.")
  }
  
  obj <- subset(obj, cells = keep_cells)
  
  DefaultAssay(obj) <- cfg$de_assay
  
  obj$latent_group <- factor(
    as.character(obj$latent_group),
    levels = c(cfg$baseline_latent_label, cfg$stimulated_latent_label)
  )
  
  obj[[cfg$age_col]] <- as.numeric(obj[[cfg$age_col, drop = TRUE]])
  
  Idents(obj) <- "latent_group"
  
  obj
}


run_findmarkers_negbinom <- function(obj, cfg) {
  markers <- FindMarkers(
    object = obj,
    ident.1 = cfg$stimulated_latent_label,
    ident.2 = cfg$baseline_latent_label,
    test.use = cfg$test_use,
    latent.vars = cfg$latent_vars,
    min.pct = cfg$min_pct,
    logfc.threshold = cfg$logfc_threshold
  )
  
  markers |>
    rownames_to_column("gene") |>
    arrange(p_val_adj)
}


prepare_glmm_inputs <- function(obj, cfg) {
  counts <- get_assay_data_compat(
    seurat_obj = obj,
    assay = cfg$raw_counts_assay,
    layer_or_slot = "counts"
  )
  
  meta <- obj@meta.data
  
  if (!identical(colnames(counts), rownames(meta))) {
    stop("RNA count matrix columns do not match metadata rownames.")
  }
  
  meta$latent_group <- factor(
    as.character(meta$latent_group),
    levels = c(cfg$baseline_latent_label, cfg$stimulated_latent_label)
  )
  
  meta$age <- as.numeric(meta[[cfg$age_col]])
  meta$sample <- factor(meta[[cfg$donor_col]])
  
  if ("nCount_RNA" %in% colnames(meta)) {
    meta$libsize <- meta$nCount_RNA
  } else {
    meta$libsize <- Matrix::colSums(counts)
  }
  
  list(
    counts = counts,
    meta = meta
  )
}


fit_one_gene_glmmTMB <- function(gene, counts, meta, cfg) {
  if (!gene %in% rownames(counts)) {
    return(tibble(
      gene = gene,
      beta = NA_real_,
      se = NA_real_,
      z = NA_real_,
      p = NA_real_,
      status = "gene_missing"
    ))
  }
  
  df <- meta
  df$y <- as.numeric(counts[gene, ])
  
  df <- df |>
    filter(
      libsize > 0,
      !is.na(age),
      !is.na(latent_group),
      !is.na(sample)
    )
  
  fit <- tryCatch(
    glmmTMB(
      y ~ latent_group + age + offset(log(libsize)) + (1 | sample),
      data = df,
      family = glmmTMB::nbinom2
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(tibble(
      gene = gene,
      beta = NA_real_,
      se = NA_real_,
      z = NA_real_,
      p = NA_real_,
      status = "model_failed"
    ))
  }
  
  sm <- summary(fit)$coefficients$cond
  
  coef_name <- paste0("latent_group", cfg$stimulated_latent_label)
  
  if (!coef_name %in% rownames(sm)) {
    return(tibble(
      gene = gene,
      beta = NA_real_,
      se = NA_real_,
      z = NA_real_,
      p = NA_real_,
      status = "coef_missing"
    ))
  }
  
  tibble(
    gene = gene,
    beta = sm[coef_name, "Estimate"],
    se = sm[coef_name, "Std. Error"],
    z = sm[coef_name, "z value"],
    p = sm[coef_name, "Pr(>|z|)"],
    status = "processed"
  )
}


select_genes_for_glmm <- function(markers, cfg) {
  sig_genes <- markers |>
    filter(
      !is.na(p_val_adj),
      p_val_adj < cfg$glmm_gene_padj_threshold,
      !is.na(avg_log2FC),
      abs(avg_log2FC) > cfg$glmm_gene_abs_logfc_threshold
    ) |>
    pull(gene) |>
    as.character() |>
    unique()
  
  if (identical(cfg$validate_top_n_genes, "all_sig")) {
    return(sig_genes)
  }
  
  head(sig_genes, as.integer(cfg$validate_top_n_genes))
}


run_glmm_validation <- function(obj, markers, cfg) {
  if (!isTRUE(cfg$run_glmm_validation)) {
    return(tibble())
  }
  
  genes_to_validate <- select_genes_for_glmm(markers, cfg)
  
  if (length(genes_to_validate) == 0) {
    warning("No genes passed thresholds for GLMM validation.")
    return(tibble())
  }
  
  inputs <- prepare_glmm_inputs(obj, cfg)
  
  glmm_res <- bind_rows(lapply(
    genes_to_validate,
    fit_one_gene_glmmTMB,
    counts = inputs$counts,
    meta = inputs$meta,
    cfg = cfg
  )) |>
    mutate(
      p_adj = p.adjust(p, method = "BH")
    ) |>
    arrange(p_adj)
  
  glmm_res
}


prepare_volcano_dataframe <- function(markers, cfg) {
  markers |>
    mutate(
      avg_log2FC = as.numeric(avg_log2FC),
      p_val_adj = as.numeric(p_val_adj)
    ) |>
    filter(
      !is.na(avg_log2FC),
      !is.na(p_val_adj)
    ) |>
    mutate(
      padj_plot = pmax(p_val_adj, 1e-300),
      neglog10 = -log10(padj_plot),
      sig = case_when(
        p_val_adj < cfg$padj_threshold &
          avg_log2FC > cfg$logfc_threshold_for_sig ~ "Up",
        p_val_adj < cfg$padj_threshold &
          avg_log2FC < -cfg$logfc_threshold_for_sig ~ "Down",
        TRUE ~ "NS"
      )
    )
}


make_volcano_plot <- function(df, cfg, labelled = FALSE) {
  p <- ggplot(df, aes(x = avg_log2FC, y = neglog10)) +
    geom_point(aes(color = sig), alpha = 0.75, size = 1.6) +
    scale_color_manual(
      values = c(
        Up = "red3",
        Down = "dodgerblue3",
        NS = "grey80"
      )
    ) +
    geom_vline(
      xintercept = c(
        -cfg$logfc_threshold_for_sig,
        cfg$logfc_threshold_for_sig
      ),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(cfg$padj_threshold),
      linetype = "dashed"
    ) +
    theme_classic() +
    theme(
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15)
    ) +
    labs(
      x = "log2 fold change",
      y = "-log10(adjusted p-value)",
      color = NULL
    )
  
  if (isTRUE(labelled)) {
    top_up <- df |>
      filter(sig == "Up") |>
      arrange(desc(avg_log2FC)) |>
      slice_head(n = cfg$n_label_up)
    
    top_down <- df |>
      filter(sig == "Down") |>
      arrange(avg_log2FC) |>
      slice_head(n = cfg$n_label_down)
    
    label_df <- bind_rows(top_up, top_down)
    
    p <- p +
      geom_label_repel(
        data = label_df,
        aes(label = gene),
        size = 3,
        max.overlaps = Inf,
        box.padding = 0.35,
        point.padding = 0.2,
        fill = "white",
        color = "black",
        label.size = 0.2,
        segment.color = "grey40"
      )
  }
  
  p
}


pick_bin <- function(x, n_bins) {
  pmin(n_bins, pmax(1, floor(x * n_bins) + 1))
}


make_waterfall_colour_table <- function(df, cfg) {
  df2 <- df |>
    mutate(
      padj_plot = pmax(as.numeric(p_val_adj), 1e-300),
      avg_log2FC = as.numeric(avg_log2FC),
      sig_score = -log10(padj_plot)
    )
  
  score_min <- min(df2$sig_score, na.rm = TRUE)
  score_max <- max(df2$sig_score, na.rm = TRUE)
  
  if (score_max == score_min) {
    df2$sig_norm <- 1
  } else {
    df2 <- df2 |>
      mutate(
        sig_norm = (sig_score - score_min) / (score_max - score_min)
      )
  }
  
  top_up <- df2 |>
    filter(sig == "Up") |>
    arrange(desc(avg_log2FC)) |>
    slice_head(n = cfg$n_waterfall_up)
  
  top_down <- df2 |>
    filter(sig == "Down") |>
    arrange(avg_log2FC) |>
    slice_head(n = cfg$n_waterfall_down)
  
  top_genes <- bind_rows(top_up, top_down) |>
    mutate(
      pal_bin = as.integer(ifelse(
        sig == "Up",
        pick_bin(sig_norm, length(cfg$red_palette)),
        pick_bin(sig_norm, length(cfg$blue_palette))
      ))
    )
  
  top_genes$color_hex <- ifelse(
    top_genes$sig == "Up",
    cfg$red_palette[top_genes$pal_bin],
    cfg$blue_palette[length(cfg$blue_palette) - top_genes$pal_bin + 1]
  )
  
  top_genes$color_hex[!(top_genes$sig %in% c("Up", "Down"))] <- "grey80"
  
  top_genes |>
    select(
      gene,
      sig,
      avg_log2FC,
      p_val_adj,
      padj_plot,
      sig_score,
      sig_norm,
      color_hex
    ) |>
    arrange(sig, desc(abs(avg_log2FC)))
}


make_de_summary <- function(markers, volcano_df, glmm_res, cfg) {
  tibble(
    metric = c(
      "n_cells_in_DE_object",
      "n_genes_tested",
      "n_significant_adj_p",
      "n_up",
      "n_down",
      "n_glmm_genes_attempted",
      "n_glmm_processed"
    ),
    value = c(
      NA,
      nrow(markers),
      sum(markers$p_val_adj < cfg$padj_threshold, na.rm = TRUE),
      sum(volcano_df$sig == "Up", na.rm = TRUE),
      sum(volcano_df$sig == "Down", na.rm = TRUE),
      nrow(glmm_res),
      if ("status" %in% colnames(glmm_res)) {
        sum(glmm_res$status == "processed", na.rm = TRUE)
      } else {
        0
      }
    )
  )
}


# ----------------------------- #
# 3. Run pipeline
# ----------------------------- #

setup_dirs(config)

message("Loading activated EBV object...")
merged <- load_seurat_object(config$input_rds)

message("Defining latent_group...")
merged <- add_latent_group(merged, config)

message("Preparing stimulated versus baseline latent DE object...")
obj <- prepare_de_object(merged, config)

message("Running Seurat FindMarkers negative-binomial DE...")
markers_nb <- run_findmarkers_negbinom(obj, config)

write.csv(
  markers_nb,
  config$findmarkers_csv,
  row.names = FALSE
)

message("Running optional glmmTMB validation...")
glmm_res <- run_glmm_validation(obj, markers_nb, config)

write.csv(
  glmm_res,
  config$glmm_validation_csv,
  row.names = FALSE
)

message("Creating volcano plots...")
volcano_df <- prepare_volcano_dataframe(markers_nb, config)

p_volcano <- make_volcano_plot(
  df = volcano_df,
  cfg = config,
  labelled = FALSE
)

ggsave(
  plot = p_volcano,
  filename = config$volcano_png,
  width = 8,
  height = 6,
  dpi = 400,
  bg = "transparent"
)

p_labelled <- make_volcano_plot(
  df = volcano_df,
  cfg = config,
  labelled = TRUE
)

ggsave(
  plot = p_labelled,
  filename = config$volcano_labelled_png,
  width = 13,
  height = 10,
  dpi = 400,
  bg = "transparent"
)

message("Creating colour-coded waterfall table...")
waterfall_table <- make_waterfall_colour_table(
  df = volcano_df,
  cfg = config
)

write.csv(
  waterfall_table,
  config$waterfall_colour_csv,
  row.names = FALSE
)

message("Saving DE summary...")
de_summary <- make_de_summary(
  markers = markers_nb,
  volcano_df = volcano_df,
  glmm_res = glmm_res,
  cfg = config
)

de_summary$value[de_summary$metric == "n_cells_in_DE_object"] <- ncol(obj)

write.csv(
  de_summary,
  config$de_summary_csv,
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  config$session_info_file
)

message("\nStimulated versus baseline latent DE analysis complete.")
message("Saved FindMarkers results to: ", config$findmarkers_csv)
message("Saved GLMM validation to: ", config$glmm_validation_csv)
message("Saved volcano plot to: ", config$volcano_png)
message("Saved labelled volcano plot to: ", config$volcano_labelled_png)
message("Saved waterfall colour table to: ", config$waterfall_colour_csv)
message("Saved DE summary to: ", config$de_summary_csv)
message("Saved session info to: ", config$session_info_file)