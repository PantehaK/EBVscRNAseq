#!/usr/bin/env Rscript

# ==============================================================================
# EBV-specific versus non-EBV T cell differential expression analysis
# ==============================================================================
#
# Purpose:
#   This script compares EBV-specific T cells against non-EBV/background T cells
#   in the baseline EBV CD8+ T cell dataset to identify an EBV-specific
#   transcriptomic signature.
#
#   It performs:
#     1. Seurat FindMarkers negative binomial differential expression,
#     2. age adjustment using latent.vars,
#     3. optional gene-level negative-binomial GLMM validation with donor random effect,
#     4. volcano plot generation,
#     5. labelled volcano plot generation,
#     6. colour-coded top gene table for waterfall/BioRender-style plotting.
#
# Expected input:
#   A module-scored baseline EBV Seurat object, for example:
#     15_baseline_EBV_module_scored.rds
#
# Notes:
#   - FindMarkers compares virus == "EBV" against virus == "negative".
#   - GLMM validation uses raw RNA counts and includes:
#       virus + age + offset(log(library size)) + (1 | sample)
#   - GLMM can be slow, so validate_top_n_genes controls how many significant
#     genes are tested.
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
  library(lme4)
  library(Matrix)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  dataset_id = "baseline_ebv",
  
  input_rds = "path/to/input/15_baseline_EBV_module_scored.rds",
  
  output_dir = "path/to/output/Publication_data/EBV_baseline/DE_EBV_vs_negative",
  
  findmarkers_csv = "path/to/output/Publication_data/EBV_baseline/DE_EBV_vs_negative/NegBinom_EBV_vs_negative_agecov.csv",
  glmm_validation_csv = "path/to/output/Publication_data/EBV_baseline/DE_EBV_vs_negative/GLMM_glmerNB_validate_EBV_vs_negative_agecov_randomDonor.csv",
  volcano_png = "path/to/output/Publication_data/EBV_baseline/DE_EBV_vs_negative/volcano_plot_EBV_genes.png",
  volcano_labelled_png = "path/to/output/Publication_data/EBV_baseline/DE_EBV_vs_negative/volcano_plot_EBV_genes_labelled.png",
  waterfall_colour_csv = "path/to/output/Publication_data/EBV_baseline/DE_EBV_vs_negative/waterfall_plot_EBV_genes_colour_coded.csv",
  de_summary_csv = "path/to/output/Publication_data/EBV_baseline/DE_EBV_vs_negative/DE_summary.csv",
  session_info_file = "path/to/output/Publication_data/EBV_baseline/DE_EBV_vs_negative/sessionInfo_EBV_vs_negative_DE.txt",
  
  # Metadata columns.
  virus_col = "virus",
  ebv_label = "EBV",
  negative_label = "negative",
  age_col = "age",
  donor_col = "sample",
  
  # Assays.
  de_assay = "SCT",
  raw_counts_assay = "RNA",
  
  # FindMarkers settings.
  test_use = "negbinom",
  latent_vars = c("age"),
  min_pct = 0.25,
  logfc_threshold = 0.5,
  
  # Significance thresholds.
  padj_threshold = 0.05,
  logfc_threshold_for_sig = 0.5,
  
  # Optional GLMM validation.
  run_glmm_validation = TRUE,
  validate_top_n_genes = 271,
  
  # Volcano label settings.
  n_label_up = 15,
  n_label_down = 15,
  
  # Waterfall output settings.
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


prepare_de_object <- function(obj, cfg) {
  required_cols <- c(cfg$virus_col, cfg$age_col, cfg$donor_col)
  missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
  
  if (length(missing_cols) > 0) {
    stop(
      "Object metadata is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  obj <- subset(
    obj,
    subset = .data[[cfg$virus_col]] %in% c(cfg$negative_label, cfg$ebv_label) &
      !is.na(.data[[cfg$age_col]])
  )
  
  obj[[cfg$virus_col]] <- factor(
    obj[[cfg$virus_col, drop = TRUE]],
    levels = c(cfg$negative_label, cfg$ebv_label)
  )
  
  obj[[cfg$age_col]] <- as.numeric(obj[[cfg$age_col, drop = TRUE]])
  
  DefaultAssay(obj) <- cfg$de_assay
  
  obj <- PrepSCTFindMarkers(obj, assay = cfg$de_assay)
  
  Idents(obj) <- cfg$virus_col
  
  obj
}


run_findmarkers_negbinom <- function(obj, cfg) {
  markers <- FindMarkers(
    object = obj,
    ident.1 = cfg$ebv_label,
    ident.2 = cfg$negative_label,
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
  
  meta$virus <- factor(
    as.character(meta[[cfg$virus_col]]),
    levels = c(cfg$negative_label, cfg$ebv_label)
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


fit_one_gene_glmm <- function(gene, counts, meta, cfg) {
  y <- as.numeric(counts[gene, ])
  
  df <- meta
  df$y <- y
  
  df <- df |>
    filter(
      libsize > 0,
      !is.na(age),
      !is.na(virus),
      !is.na(sample)
    )
  
  fit <- tryCatch(
    glmer.nb(
      y ~ virus + age + offset(log(libsize)) + (1 | sample),
      data = df,
      control = glmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    ),
    error = function(e) NULL,
    warning = function(w) {
      suppressWarnings(
        tryCatch(
          glmer.nb(
            y ~ virus + age + offset(log(libsize)) + (1 | sample),
            data = df,
            control = glmerControl(
              optimizer = "bobyqa",
              optCtrl = list(maxfun = 2e5)
            )
          ),
          error = function(e) NULL
        )
      )
    }
  )
  
  if (is.null(fit)) {
    return(tibble(
      gene = gene,
      beta_virusEBV = NA_real_,
      se = NA_real_,
      z = NA_real_,
      p = NA_real_,
      status = "model_failed"
    ))
  }
  
  sm <- summary(fit)$coefficients
  
  if (!"virusEBV" %in% rownames(sm)) {
    return(tibble(
      gene = gene,
      beta_virusEBV = NA_real_,
      se = NA_real_,
      z = NA_real_,
      p = NA_real_,
      status = "virus_term_missing"
    ))
  }
  
  tibble(
    gene = gene,
    beta_virusEBV = sm["virusEBV", "Estimate"],
    se = sm["virusEBV", "Std. Error"],
    z = sm["virusEBV", "z value"],
    p = sm["virusEBV", "Pr(>|z|)"],
    status = "processed"
  )
}


run_glmm_validation <- function(obj, markers, cfg) {
  if (!isTRUE(cfg$run_glmm_validation)) {
    return(tibble())
  }
  
  sig_genes <- markers |>
    filter(p_val_adj < cfg$padj_threshold) |>
    pull(gene)
  
  genes_to_validate <- head(sig_genes, cfg$validate_top_n_genes)
  
  if (length(genes_to_validate) == 0) {
    warning("No significant genes available for GLMM validation.")
    return(tibble())
  }
  
  glmm_inputs <- prepare_glmm_inputs(obj, cfg)
  
  glmm_res <- bind_rows(lapply(
    genes_to_validate,
    fit_one_gene_glmm,
    counts = glmm_inputs$counts,
    meta = glmm_inputs$meta,
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
      "n_genes_tested",
      "n_significant_adj_p",
      "n_up",
      "n_down",
      "n_glmm_genes_attempted",
      "n_glmm_processed"
    ),
    value = c(
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

message("Loading baseline module-scored object...")
merged <- load_seurat_object(config$input_rds)

message("Preparing EBV versus negative DE object...")
obj <- prepare_de_object(merged, config)

message("Running Seurat FindMarkers negative binomial DE...")
markers_nb <- run_findmarkers_negbinom(obj, config)

write.csv(
  markers_nb,
  config$findmarkers_csv,
  row.names = FALSE
)

message("Running optional GLMM validation...")
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
waterfall_table <- make_waterfall_colour_table(volcano_df, config)

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

write.csv(
  de_summary,
  config$de_summary_csv,
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  config$session_info_file
)

message("\nEBV versus negative DE analysis complete.")
message("Saved FindMarkers results to: ", config$findmarkers_csv)
message("Saved GLMM validation to: ", config$glmm_validation_csv)
message("Saved volcano plot to: ", config$volcano_png)
message("Saved labelled volcano plot to: ", config$volcano_labelled_png)
message("Saved waterfall colour table to: ", config$waterfall_colour_csv)
message("Saved DE summary to: ", config$de_summary_csv)
message("Saved session info to: ", config$session_info_file)