#!/usr/bin/env Rscript

# ==============================================================================
# Overall EBV transcriptional signature versus non-EBV T cells
# ==============================================================================
#
# Purpose:
#   This script identifies an overall EBV-associated transcriptional signature by
#   comparing EBV-specific T cells against non-EBV/negative T cells in the final
#   baseline EBV Seurat object.
#
#   It performs two analyses:
#
#   1. Seurat negative binomial differential expression:
#        EBV vs negative, adjusting for age.
#
#   2. Gene-level negative binomial GLMM validation:
#        y ~ virus + age + offset(log(library size)) + (1 | donor)
#
#      This validation is run on selected significant genes from the Seurat DE
#      output and uses raw RNA counts with a donor-level random intercept.
#
# Input:
#   A Seurat RDS object:
#     5_baseline_EBV_cluster_reannotated.rds
#
# Required metadata columns:
#   - virus
#   - age
#   - sample or donor ID column
#
# Notes:
#   - Paths are intentionally generic for GitHub.
#   - Update the config block before running.
#   - The GLMM step can be slow. Use config$glmm_top_n_genes to limit the number
#     of genes validated, or set it to Inf to validate all significant genes.
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
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  
  # Input RDS
  input_rds = "path/to/Publication_RDS/Baseline_EBV_integrated_seurat.rds",
  
  # Output directory
  output_dir = "path/to/Publication_data/EBV_baseline/DEG",
  
  # Metadata columns
  virus_col = "virus",
  age_col = "age",
  donor_col = "sample",
  
  # Virus labels in metadata
  ebv_label = "EBV",
  negative_label = "negative",
  
  # Seurat DE settings
  seurat_de_assay = "SCT",
  seurat_de_test = "negbinom",
  latent_vars = c("age_for_de"),
  min_pct = 0.25,
  logfc_threshold = 0.0,
  only_pos = FALSE,
  
  # GLMM settings
  glmm_assay = "RNA",
  glmm_counts_layer = "counts",
  
  # Validate the top N significant Seurat genes.
  # Use Inf to validate all FDR-significant genes.
  glmm_top_n_genes = 271,
  
  # Significance threshold used to select genes for GLMM validation
  seurat_fdr_threshold_for_glmm = 0.05,
  
  # Output files
  seurat_de_csv = "NegBinom_EBV_vs_negative_agecov.csv",
  seurat_sig_genes_csv = "NegBinom_EBV_vs_negative_agecov_FDR_significant_genes.csv",
  glmm_validation_csv = "GLMM_glmmTMB_validate_EBV_vs_negative_agecov_randomDonor.csv",
  glmm_failed_genes_csv = "GLMM_failed_genes_EBV_vs_negative_agecov_randomDonor.csv",
  filtered_metadata_csv = "EBV_vs_negative_DE_filtered_metadata.csv",
  session_info_file = "sessionInfo_EBV_vs_negative_signature.txt"
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


clean_virus <- function(x, ebv_label = "EBV", negative_label = "negative") {
  
  x_raw <- as.character(x)
  x_clean <- tolower(trimws(x_raw))
  
  dplyr::case_when(
    x_clean %in% c("ebv", "gemebv", "positive", "pos") ~ ebv_label,
    x_clean %in% c("negative", "neg", "non-ebv", "non_ebv", "none") ~ negative_label,
    TRUE ~ x_raw
  )
}


clean_donor <- function(x) {
  
  x <- trimws(as.character(x))
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


fit_one_gene_glmm <- function(gene, counts, meta) {
  
  y <- as.numeric(counts[gene, ])
  
  fit_df <- meta
  fit_df$y <- y
  
  fit_df <- fit_df %>%
    filter(
      !is.na(y),
      !is.na(virus),
      !is.na(age),
      !is.na(donor),
      !is.na(libsize),
      libsize > 0
    )
  
  if (nrow(fit_df) == 0 || length(unique(fit_df$virus)) < 2) {
    return(tibble(
      gene = gene,
      beta_virusEBV = NA_real_,
      se = NA_real_,
      z = NA_real_,
      p = NA_real_,
      status = "failed_no_usable_data"
    ))
  }
  
  if (length(unique(fit_df$y)) < 2) {
    return(tibble(
      gene = gene,
      beta_virusEBV = NA_real_,
      se = NA_real_,
      z = NA_real_,
      p = NA_real_,
      status = "failed_no_count_variation"
    ))
  }
  
  fit <- tryCatch(
    {
      glmmTMB(
        y ~ virus + age + offset(log(libsize)) + (1 | donor),
        data = fit_df,
        family = nbinom2
      )
    },
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    return(tibble(
      gene = gene,
      beta_virusEBV = NA_real_,
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
  
  if (is.null(sm) || !("virusEBV" %in% rownames(sm))) {
    return(tibble(
      gene = gene,
      beta_virusEBV = NA_real_,
      se = NA_real_,
      z = NA_real_,
      p = NA_real_,
      status = "failed_missing_virusEBV_coefficient"
    ))
  }
  
  tibble(
    gene = gene,
    beta_virusEBV = sm["virusEBV", "Estimate"],
    se = sm["virusEBV", "Std. Error"],
    z = sm["virusEBV", "z value"],
    p = sm["virusEBV", "Pr(>|z|)"],
    status = "ok"
  )
}


# ----------------------------- #
# 3. Load object and prepare metadata
# ----------------------------- #

create_dir(config$output_dir)

message("Loading baseline EBV object:")
message(config$input_rds)

merged <- load_seurat_object(config$input_rds)

required_metadata_cols <- c(
  config$virus_col,
  config$age_col,
  config$donor_col
)

check_required_columns(
  df = merged@meta.data,
  required_cols = required_metadata_cols,
  object_name = "baseline EBV metadata"
)

meta_clean <- merged@meta.data %>%
  as.data.frame() %>%
  mutate(
    virus_for_de = clean_virus(
      .data[[config$virus_col]],
      ebv_label = config$ebv_label,
      negative_label = config$negative_label
    ),
    age_for_de = suppressWarnings(as.numeric(.data[[config$age_col]])),
    donor_for_de = clean_donor(.data[[config$donor_col]])
  )

rownames(meta_clean) <- colnames(merged)
merged@meta.data <- meta_clean


# ----------------------------- #
# 4. Subset EBV and negative cells
# ----------------------------- #

obj <- subset(
  merged,
  subset = virus_for_de %in% c(config$ebv_label, config$negative_label) &
    !is.na(age_for_de) &
    !is.na(donor_for_de)
)

if (ncol(obj) == 0) {
  stop("No cells remained after EBV/negative/age/donor filtering.")
}

obj$virus_for_de <- factor(
  obj$virus_for_de,
  levels = c(config$negative_label, config$ebv_label)
)

obj$age_for_de <- as.numeric(obj$age_for_de)
obj$donor_for_de <- factor(obj$donor_for_de)

message("Cells retained for EBV vs negative DE:")
print(table(obj$virus_for_de, useNA = "ifany"))

message("Donors retained for EBV vs negative DE:")
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

# PrepSCTFindMarkers is only needed for SCT. Skip safely for other assays.
if (config$seurat_de_assay == "SCT") {
  message("Running PrepSCTFindMarkers...")
  obj <- PrepSCTFindMarkers(obj)
}

Idents(obj) <- "virus_for_de"

message("Running Seurat FindMarkers...")
markers_nb <- FindMarkers(
  object = obj,
  ident.1 = config$ebv_label,
  ident.2 = config$negative_label,
  test.use = config$seurat_de_test,
  latent.vars = config$latent_vars,
  min.pct = config$min_pct,
  logfc.threshold = config$logfc_threshold,
  only.pos = config$only_pos
)

markers_nb <- markers_nb %>%
  rownames_to_column("gene") %>%
  arrange(p_val_adj, p_val)

seurat_de_path <- file.path(
  config$output_dir,
  config$seurat_de_csv
)

write_csv(markers_nb, seurat_de_path)

sig_genes_tbl <- markers_nb %>%
  filter(!is.na(p_val_adj), p_val_adj < config$seurat_fdr_threshold_for_glmm)

sig_genes <- sig_genes_tbl %>%
  pull(gene)

sig_genes_path <- file.path(
  config$output_dir,
  config$seurat_sig_genes_csv
)

write_csv(sig_genes_tbl, sig_genes_path)

message("Number of FDR-significant Seurat genes: ", length(sig_genes))

if (length(sig_genes) > 0) {
  message("Top significant genes:")
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
    virus = factor(
      as.character(virus_for_de),
      levels = c(config$negative_label, config$ebv_label)
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
  warning("No FDR-significant Seurat genes were available in the GLMM counts matrix.")
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
      meta = meta_glmm
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
    beta_virusEBV = numeric(),
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
    "Seurat RDS metadata and assays from 5_baseline_EBV_cluster_reannotated.rds",
    "",
    "Cells retained for DE:",
    capture.output(print(table(obj$virus_for_de, useNA = "ifany"))),
    "",
    "Number of donors retained:",
    as.character(length(unique(obj$donor_for_de))),
    "",
    "Number of Seurat FDR-significant genes:",
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

message("\nOverall EBV vs negative signature workflow complete.")
message("Filtered metadata: ", filtered_metadata_path)
message("Seurat DE results: ", seurat_de_path)
message("Seurat significant genes: ", sig_genes_path)
message("GLMM validation results: ", glmm_validation_path)
message("GLMM failed genes: ", failed_genes_path)
message("Session info: ", session_info_path)
