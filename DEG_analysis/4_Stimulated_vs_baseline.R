#!/usr/bin/env Rscript

# ==============================================================================
# Stimulated versus baseline latent EBV-specific T cell transcriptional signature
# ==============================================================================
#
# Purpose:
#   This script identifies transcriptional differences between stimulated latent
#   EBV-specific T cells and baseline latent EBV-specific T cells using the
#   activated/baseline integrated Seurat object.
#
#   It performs two analyses:
#
#   1. Seurat negative binomial differential expression:
#        Stimulated Latent vs Baseline Latent, adjusting for age.
#
#   2. Gene-level negative binomial GLMM validation:
#        y ~ latent_group + age + offset(log(library size)) + (1 | donor)
#
#      This validation is run on selected significant genes from the Seurat DE
#      output and uses raw RNA counts with a donor-level random intercept.
#
# Input:
#   A Seurat RDS object containing both GEMEBV baseline and EBVLCL stimulated
#   cells, for example:
#     3_activated_EBV_module_scored_reannotated.rds
#
# Required metadata columns:
#   - batch
#   - lifecycle
#   - age
#   - sample or donor ID column
#
# Latent group definition:
#   Baseline Latent:
#     batch starts with GEMEBV AND lifecycle == Latent
#
#   Stimulated Latent:
#     batch starts with EBVLCL AND lifecycle == Latent
#
#   Non-Latent:
#     all other cells
#
# Notes:
#   - Paths are intentionally generic for GitHub.
#   - Update the config block before running.
#   - The GLMM step can be slow. Use config$glmm_top_n_genes to limit the number
#     of genes validated, or set it to Inf to validate all selected significant
#     genes.
#
# ==============================================================================


# ----------------------------- #
# 0. Load packages
# ----------------------------- #

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(Matrix)
  library(glmmTMB)
  library(readr)
  library(stringr)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  
  # Input RDS
  input_rds = "path/to/Publication_RDS/Activated_EBV_integrated_seurat.rds",
  
  # Output directory
  output_dir = "path/to/EBV_activated/DEG",
  
  # Metadata columns
  batch_col = "batch",
  lifecycle_col = "lifecycle",
  age_col = "age",
  donor_col = "sample",
  
  # Group assignment settings
  baseline_batch_pattern = "^GEMEBV",
  stimulated_batch_pattern = "^EBVLCL",
  latent_lifecycle_label = "Latent",
  baseline_label = "Baseline Latent",
  stimulated_label = "Stimulated Latent",
  non_latent_label = "Non-Latent",
  
  # Seurat DE settings
  seurat_de_assay = "RNA",
  seurat_de_test = "negbinom",
  latent_vars = c("age_for_de"),
  min_pct = 0.25,
  logfc_threshold = 0.0,
  only_pos = FALSE,
  
  # GLMM settings
  glmm_assay = "RNA",
  glmm_counts_layer = "counts",
  
  # Genes selected for GLMM validation from Seurat output
  seurat_fdr_threshold_for_glmm = 0.05,
  seurat_abs_logfc_threshold_for_glmm = 2,
  
  # Validate top N selected genes.
  # Use Inf to validate all selected genes.
  glmm_top_n_genes = Inf,
  
  # Output files
  filtered_metadata_csv = "Stimulated_vs_Baseline_Latent_filtered_metadata.csv",
  seurat_de_csv = "NegBinom_Stimulated_vs_Baseline_agecov.csv",
  seurat_sig_genes_csv = "NegBinom_Stimulated_vs_Baseline_agecov_FDR_logFC_significant_genes.csv",
  glmm_validation_csv = "GLMM_glmmTMB_Stimulated_vs_Baseline_agecov_randomDonor.csv",
  glmm_failed_genes_csv = "GLMM_failed_genes_Stimulated_vs_Baseline_agecov_randomDonor.csv",
  session_info_file = "sessionInfo_Stimulated_vs_Baseline_Latent_signature.txt"
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


get_assay_data_compat <- function(object, assay, data_layer = "counts") {
  
  # Compatible with Seurat v5 and Seurat v4.
  tryCatch(
    {
      GetAssayData(object, assay = assay, layer = data_layer)
    },
    error = function(e) {
      GetAssayData(object, assay = assay, slot = data_layer)
    }
  )
}


clean_chr <- function(x) {
  
  x <- str_trim(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL", "null", "None", "none")] <- NA_character_
  x
}


format_p <- function(p) {
  
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}


get_logfc_column <- function(markers_df) {
  
  possible_cols <- c("avg_log2FC", "avg_logFC", "logFC")
  hit <- intersect(possible_cols, colnames(markers_df))
  
  if (length(hit) == 0) {
    stop(
      "Could not find a log fold-change column. Tried: ",
      paste(possible_cols, collapse = ", ")
    )
  }
  
  hit[1]
}


get_latent_group_coef <- function(coef_table, stimulated_label = "Stimulated Latent") {
  
  # glmmTMB usually names the coefficient:
  #   latent_groupStimulated Latent
  # but this helper avoids hard-coding if R changes formatting.
  exact_name <- paste0("latent_group", stimulated_label)
  
  if (exact_name %in% rownames(coef_table)) {
    return(exact_name)
  }
  
  candidates <- rownames(coef_table)[
    str_detect(rownames(coef_table), "^latent_group")
  ]
  
  candidates <- candidates[
    str_detect(candidates, fixed(stimulated_label))
  ]
  
  if (length(candidates) == 1) {
    return(candidates)
  }
  
  NA_character_
}


fit_one_gene_glmm <- function(gene, counts, meta, stimulated_label) {
  
  if (!(gene %in% rownames(counts))) {
    return(tibble(
      gene = gene,
      beta = NA_real_,
      se = NA_real_,
      z = NA_real_,
      p = NA_real_,
      status = "failed_gene_not_in_counts"
    ))
  }
  
  y <- as.numeric(counts[gene, ])
  
  fit_df <- meta
  fit_df$y <- y
  
  fit_df <- fit_df %>%
    filter(
      !is.na(y),
      !is.na(latent_group),
      !is.na(age),
      !is.na(donor),
      !is.na(libsize),
      libsize > 0
    )
  
  if (nrow(fit_df) == 0 || length(unique(fit_df$latent_group)) < 2) {
    return(tibble(
      gene = gene,
      beta = NA_real_,
      se = NA_real_,
      z = NA_real_,
      p = NA_real_,
      status = "failed_no_usable_data"
    ))
  }
  
  if (length(unique(fit_df$y)) < 2) {
    return(tibble(
      gene = gene,
      beta = NA_real_,
      se = NA_real_,
      z = NA_real_,
      p = NA_real_,
      status = "failed_no_count_variation"
    ))
  }
  
  fit <- tryCatch(
    {
      glmmTMB(
        y ~ latent_group + age + offset(log(libsize)) + (1 | donor),
        data = fit_df,
        family = nbinom2
      )
    },
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    return(tibble(
      gene = gene,
      beta = NA_real_,
      se = NA_real_,
      z = NA_real_,
      p = NA_real_,
      status = paste0("failed_model_error: ", conditionMessage(fit))
    ))
  }
  
  sm <- tryCatch(
    {
      summary(fit)$coefficients$cond
    },
    error = function(e) NULL
  )
  
  if (is.null(sm)) {
    return(tibble(
      gene = gene,
      beta = NA_real_,
      se = NA_real_,
      z = NA_real_,
      p = NA_real_,
      status = "failed_no_coefficient_table"
    ))
  }
  
  coef_name <- get_latent_group_coef(
    coef_table = sm,
    stimulated_label = stimulated_label
  )
  
  if (is.na(coef_name) || !(coef_name %in% rownames(sm))) {
    return(tibble(
      gene = gene,
      beta = NA_real_,
      se = NA_real_,
      z = NA_real_,
      p = NA_real_,
      status = "failed_missing_stimulated_latent_coefficient"
    ))
  }
  
  tibble(
    gene = gene,
    beta = sm[coef_name, "Estimate"],
    se = sm[coef_name, "Std. Error"],
    z = sm[coef_name, "z value"],
    p = sm[coef_name, "Pr(>|z|)"],
    status = "ok"
  )
}


# ----------------------------- #
# 3. Load object and prepare metadata
# ----------------------------- #

create_dir(config$output_dir)

message("Loading activated EBV object:")
message(config$input_rds)

merged <- load_seurat_object(config$input_rds)

required_metadata_cols <- c(
  config$batch_col,
  config$lifecycle_col,
  config$age_col,
  config$donor_col
)

check_required_columns(
  df = merged@meta.data,
  required_cols = required_metadata_cols,
  object_name = "activated EBV metadata"
)

meta_clean <- merged@meta.data %>%
  as.data.frame() %>%
  mutate(
    batch_for_de = clean_chr(.data[[config$batch_col]]),
    lifecycle_for_de = clean_chr(.data[[config$lifecycle_col]]),
    age_for_de = suppressWarnings(as.numeric(.data[[config$age_col]])),
    donor_for_de = clean_chr(.data[[config$donor_col]]),
    
    latent_group = case_when(
      !is.na(batch_for_de) &
        grepl(config$baseline_batch_pattern, batch_for_de) &
        lifecycle_for_de == config$latent_lifecycle_label ~ config$baseline_label,
      
      !is.na(batch_for_de) &
        grepl(config$stimulated_batch_pattern, batch_for_de) &
        lifecycle_for_de == config$latent_lifecycle_label ~ config$stimulated_label,
      
      TRUE ~ config$non_latent_label
    ),
    
    latent_group = factor(
      latent_group,
      levels = c(
        config$non_latent_label,
        config$baseline_label,
        config$stimulated_label
      )
    )
  )

rownames(meta_clean) <- colnames(merged)
merged@meta.data <- meta_clean

message("Latent group counts:")
print(table(merged$latent_group, useNA = "ifany"))


# ----------------------------- #
# 4. Subset baseline and stimulated latent cells
# ----------------------------- #

obj <- subset(
  merged,
  subset = latent_group %in% c(config$stimulated_label, config$baseline_label) &
    !is.na(age_for_de) &
    !is.na(donor_for_de)
)

if (ncol(obj) == 0) {
  stop("No cells remained after baseline/stimulated latent filtering.")
}

obj$latent_group <- factor(
  as.character(obj$latent_group),
  levels = c(config$baseline_label, config$stimulated_label)
)

obj$age_for_de <- as.numeric(obj$age_for_de)
obj$donor_for_de <- factor(obj$donor_for_de)

message("Cells retained for stimulated vs baseline latent DE:")
print(table(obj$latent_group, useNA = "ifany"))

message("Donors retained for stimulated vs baseline latent DE:")
print(length(unique(obj$donor_for_de)))

filtered_metadata_path <- file.path(
  config$output_dir,
  config$filtered_metadata_csv
)

write_csv(
  obj@meta.data %>%
    as.data.frame() %>%
    rownames_to_column("cell_barcode"),
  filtered_metadata_path
)


# ----------------------------- #
# 5. Seurat negative binomial DE
# ----------------------------- #

if (!config$seurat_de_assay %in% Assays(obj)) {
  stop(
    "Seurat DE assay not found: ",
    config$seurat_de_assay,
    ". Available assays: ",
    paste(Assays(obj), collapse = ", ")
  )
}

DefaultAssay(obj) <- config$seurat_de_assay
Idents(obj) <- "latent_group"

message("Running Seurat FindMarkers...")
markers_nb <- FindMarkers(
  object = obj,
  ident.1 = config$stimulated_label,
  ident.2 = config$baseline_label,
  test.use = config$seurat_de_test,
  latent.vars = config$latent_vars,
  min.pct = config$min_pct,
  logfc.threshold = config$logfc_threshold,
  only.pos = config$only_pos
)

markers_nb <- markers_nb %>%
  rownames_to_column("gene") %>%
  arrange(p_val_adj, p_val)

logfc_col <- get_logfc_column(markers_nb)

seurat_de_path <- file.path(
  config$output_dir,
  config$seurat_de_csv
)

write_csv(markers_nb, seurat_de_path)

sig_genes_tbl <- markers_nb %>%
  filter(
    !is.na(p_val_adj),
    p_val_adj < config$seurat_fdr_threshold_for_glmm,
    !is.na(.data[[logfc_col]]),
    abs(.data[[logfc_col]]) > config$seurat_abs_logfc_threshold_for_glmm
  )

sig_genes <- sig_genes_tbl %>%
  pull(gene) %>%
  as.character() %>%
  unique()

sig_genes_path <- file.path(
  config$output_dir,
  config$seurat_sig_genes_csv
)

write_csv(sig_genes_tbl, sig_genes_path)

message("Number of selected Seurat genes for GLMM: ", length(sig_genes))

if (length(sig_genes) > 0) {
  message("Top selected genes:")
  print(head(sig_genes, 20))
}


# ----------------------------- #
# 6. GLMM validation using raw RNA counts
# ----------------------------- #

if (!config$glmm_assay %in% Assays(obj)) {
  stop(
    "GLMM assay not found: ",
    config$glmm_assay,
    ". Available assays: ",
    paste(Assays(obj), collapse = ", ")
  )
}

counts <- get_assay_data_compat(
  object = obj,
  assay = config$glmm_assay,
  data_layer = config$glmm_counts_layer
)

meta_glmm <- obj@meta.data %>%
  as.data.frame()

if (!identical(colnames(counts), rownames(meta_glmm))) {
  stop("Counts columns and metadata rownames are not aligned.")
}

meta_glmm <- meta_glmm %>%
  mutate(
    latent_group = factor(
      as.character(latent_group),
      levels = c(config$baseline_label, config$stimulated_label)
    ),
    age = as.numeric(age_for_de),
    donor = factor(donor_for_de)
  )

if ("nCount_RNA" %in% colnames(meta_glmm)) {
  meta_glmm$libsize <- meta_glmm$nCount_RNA
} else {
  meta_glmm$libsize <- Matrix::colSums(counts)
}

genes_available_for_glmm <- intersect(sig_genes, rownames(counts))

if (length(genes_available_for_glmm) == 0) {
  warning("No selected Seurat genes were available in the GLMM counts matrix.")
  genes_to_validate <- character(0)
} else {
  if (is.infinite(config$glmm_top_n_genes)) {
    genes_to_validate <- genes_available_for_glmm
  } else {
    genes_to_validate <- head(
      genes_available_for_glmm,
      config$glmm_top_n_genes
    )
  }
}

message("Genes selected for GLMM validation: ", length(genes_to_validate))

if (length(genes_to_validate) > 0) {
  
  glmm_res <- bind_rows(
    lapply(
      genes_to_validate,
      fit_one_gene_glmm,
      counts = counts,
      meta = meta_glmm,
      stimulated_label = config$stimulated_label
    )
  ) %>%
    mutate(
      p_adj = p.adjust(p, method = "BH"),
      p_label = format_p(p),
      p_adj_label = format_p(p_adj)
    ) %>%
    arrange(p_adj, p)
  
} else {
  
  glmm_res <- tibble(
    gene = character(),
    beta = numeric(),
    se = numeric(),
    z = numeric(),
    p = numeric(),
    status = character(),
    p_adj = numeric(),
    p_label = character(),
    p_adj_label = character()
  )
}

glmm_validation_path <- file.path(
  config$output_dir,
  config$glmm_validation_csv
)

write_csv(glmm_res, glmm_validation_path)

failed_genes <- glmm_res %>%
  filter(status != "ok")

failed_genes_path <- file.path(
  config$output_dir,
  config$glmm_failed_genes_csv
)

write_csv(failed_genes, failed_genes_path)


# ----------------------------- #
# 7. Save session information
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
    "Input source:",
    "Seurat RDS metadata and assays from activated/baseline integrated object",
    "",
    "Latent group counts:",
    capture.output(print(table(merged$latent_group, useNA = "ifany"))),
    "",
    "Cells retained for DE:",
    capture.output(print(table(obj$latent_group, useNA = "ifany"))),
    "",
    "Number of donors retained:",
    as.character(length(unique(obj$donor_for_de))),
    "",
    "Seurat logFC column used:",
    logfc_col,
    "",
    "Number of selected Seurat genes for GLMM:",
    as.character(length(sig_genes)),
    "",
    "Number of genes selected for GLMM validation:",
    as.character(length(genes_to_validate)),
    "",
    "Output files:",
    filtered_metadata_path,
    seurat_de_path,
    sig_genes_path,
    glmm_validation_path,
    failed_genes_path,
    "",
    "Session information:",
    capture.output(sessionInfo())
  ),
  session_info_path
)


# ----------------------------- #
# 8. Completion messages
# ----------------------------- #

message("\nStimulated vs baseline latent signature workflow complete.")
message("Filtered metadata: ", filtered_metadata_path)
message("Seurat DE results: ", seurat_de_path)
message("Seurat selected significant genes: ", sig_genes_path)
message("GLMM validation results: ", glmm_validation_path)
message("GLMM failed genes: ", failed_genes_path)
message("Session info: ", session_info_path)
