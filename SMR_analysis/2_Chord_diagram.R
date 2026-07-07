#!/usr/bin/env Rscript

# ==============================================================================
# SMR overlap analysis with stimulated latent T-cell DE genes
# ==============================================================================
#
# Purpose:
#   This script integrates SMR-prioritised genes from multiple eQTL resources with
#   pseudobulk differential expression results from stimulated latent EBV-specific
#   T cells, then generates:
#
#     1. Filtered SMR gene tables.
#     2. Overlap tables between SMR genes and stimulated latent DE genes.
#     3. Chord-diagram edge tables linking:
#          SMR source -> gene -> stimulated latent T-cell cluster.
#     4. A chord diagram PDF.
#     5. Summary tables for manuscript/supplementary reporting.
#
# SMR data availability:
#   The SMR analysis data used as input for this workflow are available on:
#   10.5281/zenodo.21231045
#
# Pseudobulk input:
#   The stimulated latent MS-vs-Control pseudobulk result table should contain
#   at least:
#     - gene
#     - matched_cluster or cluster
#     - p_val_adj
#
# SMR input:
#   Each SMR input table should contain:
#     - Gene
#     - p_FDR
#     - a chromosome column
#     - a base-pair/position column
#     - a HEIDI p-value column
#
# Filtering:
#   SMR associations are retained if:
#     - p_FDR < config$smr_fdr_cutoff
#     - p_HEIDI > config$heidi_cutoff
#
#   Pseudobulk genes are retained if:
#     - p_val_adj < config$pseudobulk_fdr_cutoff
#
# Notes:
#   - Paths are intentionally generic for GitHub.
#   - Update the config block before running.
#   - The chord diagram is skipped gracefully if no overlap genes are detected.
#
# ==============================================================================


# ----------------------------- #
# 0. Load packages
# ----------------------------- #

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readxl)
  library(readr)
  library(stringr)
  library(circlize)
  library(openxlsx)
})


# ----------------------------- #
# 1. User configuration
# ----------------------------- #

config <- list(
  
  # Input directory containing SMR files.
  # SMR analysis data are available on Zenodo.org.
  smr_input_dir = "path/to/SMR_analysis",
  
  # Output directory
  output_dir = "outputs/SMR_overlap_stimulated_latent",
  
  # SMR input files
  smr_files = list(
    eqtlgen = list(
      file = "eQTLgen.fullSMR.FDR",
      source_label = "Blood eQTL (eQTLGen)"
    ),
    blood_scrnaseq = list(
      file = "blood_scRNAseq.fullSMR.FDR",
      source_label = "Immune cell eQTL (1k1k)"
    ),
    brain_scrnaseq = list(
      file = "brain_scRNAseq.fullSMR.FDR",
      source_label = "Brain cell eQTL (snRNA-seq)"
    ),
    brain_meta = list(
      file = "BrainMetaV2.fullSMR.FDR",
      source_label = "Brain eQTL (bulk cortex)"
    )
  ),
  
  # Stimulated latent pseudobulk DE table.
  # This table may be the same as Table S10 in the Supplementary Tables.
  pseudobulk_file = "path/to/Stimulated_MS_vs_control.xlsx",
  pseudobulk_sheet = 1,
  
  # Pseudobulk columns
  pseudobulk_gene_col = "gene",
  pseudobulk_cluster_col_candidates = c("matched_cluster", "cluster", "new_cluster"),
  pseudobulk_fdr_col = "p_val_adj",
  
  # Filtering thresholds
  smr_fdr_cutoff = 0.05,
  heidi_cutoff = 0.01,
  pseudobulk_fdr_cutoff = 0.05,
  
  # Clusters included in chord diagram
  clusters_to_keep = c("1", "0", "2", "6", "7"),
  
  # MHC region annotation
  mhc_chr = 6,
  mhc_start = 24000000,
  mhc_end = 32000000,
  
  # Chord diagram settings
  chord_pdf_width = 14,
  chord_pdf_height = 14,
  
  # Output files
  smr_all_filtered_csv = "SMR_all_sources_FDR_HEIDI_filtered_MHC.csv",
  smr_gene_summary_csv = "SMR_genes_passing_FDR_HEIDI_all_sources.csv",
  smr_summary_counts_csv = "SMR_FDR_HEIDI_summary_counts.csv",
  smr_supplementary_csv = "Supplementary_Table_SMR_FDR_HEIDI_passing_associations.csv",
  
  pseudobulk_filtered_csv = "Stimulated_latent_pseudobulk_DE_genes_FDR_filtered.csv",
  overlap_genes_csv = "combined_overlap_genes_pseudobulkFDR_SMR_FDR_HEIDI_MHC.csv",
  chord_edges_csv = "combined_chord_edges_pseudobulkFDR_SMR_FDR_HEIDI_MHC.csv",
  gene_source_cluster_summary_csv = "gene_source_cluster_summary_pseudobulkFDR_SMR_FDR_HEIDI_MHC.csv",
  gene_order_csv = "gene_order_for_chord_diagram.csv",
  
  chord_pdf = "Chord_combined_sources_genes_clusters_pseudobulkFDR_SMR_FDR_HEIDI_MHC.pdf",
  workbook_xlsx = "SMR_overlap_stimulated_latent_summary_tables.xlsx",
  session_info_file = "sessionInfo_SMR_overlap_stimulated_latent.txt"
)


# ----------------------------- #
# 2. Helper functions
# ----------------------------- #

create_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}


read_delimited_auto <- function(file_path) {
  
  if (!file.exists(file_path)) {
    stop("Input file does not exist: ", file_path)
  }
  
  # The original files were comma-separated despite their extension.
  # readr::read_delim() with delim = "," keeps behaviour explicit.
  readr::read_delim(
    file_path,
    delim = ",",
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    as.data.frame()
}


read_pseudobulk_table <- function(file_path, sheet = 1) {
  
  if (!file.exists(file_path)) {
    stop("Pseudobulk file does not exist: ", file_path)
  }
  
  ext <- tolower(tools::file_ext(file_path))
  
  if (ext %in% c("xlsx", "xls")) {
    readxl::read_excel(file_path, sheet = sheet) %>%
      as.data.frame()
  } else if (ext %in% c("csv")) {
    readr::read_csv(file_path, show_col_types = FALSE) %>%
      as.data.frame()
  } else if (ext %in% c("tsv", "txt")) {
    readr::read_tsv(file_path, show_col_types = FALSE) %>%
      as.data.frame()
  } else {
    stop(
      "Unsupported pseudobulk file extension: ",
      ext,
      ". Use .xlsx, .xls, .csv, .tsv or .txt."
    )
  }
}


find_col_any <- function(df, possible_names, label = "column") {
  
  hits <- names(df)[tolower(names(df)) %in% tolower(possible_names)]
  
  if (length(hits) == 0) {
    stop(
      "Missing ",
      label,
      ". Tried: ",
      paste(possible_names, collapse = ", ")
    )
  }
  
  hits[1]
}


to_numeric_clean <- function(x) {
  
  z <- trimws(as.character(x))
  z[z %in% c("", "NA", "N/A", "na", "n/a", "NaN", "nan", "null", "NULL")] <- NA
  z <- gsub(",", "", z)
  z <- gsub("<", "", z)
  z <- gsub(">", "", z)
  
  suppressWarnings(as.numeric(z))
}


safe_file_path <- function(dir, filename) {
  file.path(dir, filename)
}


# ----------------------------- #
# 3. Set output paths
# ----------------------------- #

create_dir(config$output_dir)

out_smr_all_filtered <- safe_file_path(config$output_dir, config$smr_all_filtered_csv)
out_smr_gene_summary <- safe_file_path(config$output_dir, config$smr_gene_summary_csv)
out_smr_summary_counts <- safe_file_path(config$output_dir, config$smr_summary_counts_csv)
out_smr_supplementary <- safe_file_path(config$output_dir, config$smr_supplementary_csv)

out_pseudobulk_filtered <- safe_file_path(config$output_dir, config$pseudobulk_filtered_csv)
out_overlap_genes <- safe_file_path(config$output_dir, config$overlap_genes_csv)
out_chord_edges <- safe_file_path(config$output_dir, config$chord_edges_csv)
out_gene_summary <- safe_file_path(config$output_dir, config$gene_source_cluster_summary_csv)
out_gene_order <- safe_file_path(config$output_dir, config$gene_order_csv)

out_chord_pdf <- safe_file_path(config$output_dir, config$chord_pdf)
out_workbook <- safe_file_path(config$output_dir, config$workbook_xlsx)
out_session_info <- safe_file_path(config$output_dir, config$session_info_file)


# ----------------------------- #
# 4. Load input data
# ----------------------------- #

message("Loading SMR input files from: ", config$smr_input_dir)

smr_input_list <- lapply(config$smr_files, function(x) {
  
  file_path <- file.path(config$smr_input_dir, x$file)
  
  list(
    data = read_delimited_auto(file_path),
    source_label = x$source_label,
    file_path = file_path
  )
})

message("Loading pseudobulk DE table: ", config$pseudobulk_file)

sc <- read_pseudobulk_table(
  file_path = config$pseudobulk_file,
  sheet = config$pseudobulk_sheet
)

names(sc) <- trimws(names(sc))

cluster_col <- find_col_any(
  df = sc,
  possible_names = config$pseudobulk_cluster_col_candidates,
  label = "pseudobulk cluster column"
)

if (!config$pseudobulk_gene_col %in% names(sc)) {
  stop("Pseudobulk gene column not found: ", config$pseudobulk_gene_col)
}

if (!config$pseudobulk_fdr_col %in% names(sc)) {
  stop("Pseudobulk FDR column not found: ", config$pseudobulk_fdr_col)
}

sc <- sc %>%
  rename(
    gene = all_of(config$pseudobulk_gene_col),
    cluster = all_of(cluster_col),
    p_val_adj = all_of(config$pseudobulk_fdr_col)
  ) %>%
  mutate(
    gene = as.character(gene),
    cluster = as.character(cluster),
    p_val_adj = to_numeric_clean(p_val_adj)
  )


# ----------------------------- #
# 5. Filter SMR results
# ----------------------------- #

filter_smr <- function(df, source_label) {
  
  names(df) <- trimws(names(df))
  
  gene_col <- find_col_any(
    df,
    possible_names = c("Gene", "gene", "GENE", "Probe", "probe"),
    label = paste0("gene column for ", source_label)
  )
  
  fdr_col <- find_col_any(
    df,
    possible_names = c("p_FDR", "P_FDR", "FDR", "p_adj", "padj"),
    label = paste0("SMR FDR column for ", source_label)
  )
  
  chr_col <- find_col_any(
    df,
    possible_names = c("Chr", "CHR", "chr", "SNPChr", "GeneChr", "ProbeChr"),
    label = paste0("chromosome column for ", source_label)
  )
  
  bp_col <- find_col_any(
    df,
    possible_names = c("BP", "bp", "SNPBP", "GeneBP", "Probe_bp", "ProbeBP", "probe_bp", "Position", "pos"),
    label = paste0("base-pair position column for ", source_label)
  )
  
  heidi_col <- find_col_any(
    df,
    possible_names = c("p_HEIDI", "P_HEIDI", "pheidi", "PHEIDI", "p_heidi", "P_heidi"),
    label = paste0("HEIDI p-value column for ", source_label)
  )
  
  df %>%
    mutate(
      gene_use = as.character(.data[[gene_col]]),
      p_fdr_use = to_numeric_clean(.data[[fdr_col]]),
      chr_use = to_numeric_clean(.data[[chr_col]]),
      bp_use = to_numeric_clean(.data[[bp_col]]),
      pheidi_use = to_numeric_clean(.data[[heidi_col]]),
      is_MHC = chr_use == config$mhc_chr &
        bp_use >= config$mhc_start &
        bp_use <= config$mhc_end
    ) %>%
    filter(
      !is.na(gene_use),
      gene_use != "",
      !is.na(p_fdr_use),
      p_fdr_use < config$smr_fdr_cutoff,
      !is.na(pheidi_use),
      pheidi_use > config$heidi_cutoff
    ) %>%
    transmute(
      gene = gene_use,
      source = source_label,
      p_FDR = p_fdr_use,
      p_HEIDI = pheidi_use,
      chr = chr_use,
      bp = bp_use,
      region = ifelse(is_MHC, "MHC_HEIDI_pass", "nonMHC_HEIDI_pass")
    )
}

smr_all_sig <- purrr::map_dfr(
  smr_input_list,
  ~ filter_smr(.x$data, .x$source_label)
) %>%
  distinct()

write_csv(smr_all_sig, out_smr_all_filtered)

message("SMR rows passing FDR + HEIDI rules: ", nrow(smr_all_sig))
message("SMR genes passing FDR + HEIDI rules: ", length(unique(smr_all_sig$gene)))
message("Region breakdown:")
print(table(smr_all_sig$region))


# ----------------------------- #
# 6. Filter pseudobulk DE results
# ----------------------------- #

sc_sig <- sc %>%
  filter(
    !is.na(gene),
    gene != "",
    !is.na(cluster),
    !is.na(p_val_adj),
    p_val_adj < config$pseudobulk_fdr_cutoff
  )

write_csv(sc_sig, out_pseudobulk_filtered)

message("Pseudobulk rows passing FDR < ", config$pseudobulk_fdr_cutoff, ": ", nrow(sc_sig))
message("Cluster distribution:")
print(table(sc_sig$cluster))


# ----------------------------- #
# 7. Overlap genes and chord edges
# ----------------------------- #

overlap_genes <- intersect(
  unique(sc_sig$gene),
  unique(smr_all_sig$gene)
)

message("Number of overlapping genes: ", length(overlap_genes))

write_csv(
  tibble(gene = sort(overlap_genes)),
  out_overlap_genes
)

source_links <- smr_all_sig %>%
  filter(gene %in% overlap_genes) %>%
  distinct(source, gene, region) %>%
  mutate(value = 1)

cluster_links <- sc_sig %>%
  filter(
    gene %in% overlap_genes,
    cluster %in% config$clusters_to_keep
  ) %>%
  mutate(cluster_label = paste0("Cluster ", cluster)) %>%
  distinct(cluster_label, gene) %>%
  mutate(value = 1)

genes_with_source <- unique(source_links$gene)
genes_with_cluster <- unique(cluster_links$gene)

genes_to_plot <- intersect(
  genes_with_source,
  genes_with_cluster
)

source_links <- source_links %>%
  filter(gene %in% genes_to_plot)

cluster_links <- cluster_links %>%
  filter(gene %in% genes_to_plot)

message("Genes with both SMR source and cluster links: ", length(genes_to_plot))

df_source <- source_links %>%
  transmute(
    from = source,
    to = gene,
    value = value
  )

df_cluster <- cluster_links %>%
  transmute(
    from = gene,
    to = cluster_label,
    value = value
  )

plot_df <- bind_rows(
  df_source,
  df_cluster
)

write_csv(plot_df, out_chord_edges)


# ----------------------------- #
# 8. Ordering and colours
# ----------------------------- #

source_order <- c(
  "Blood eQTL (eQTLGen)",
  "Immune cell eQTL (1k1k)",
  "Brain cell eQTL (snRNA-seq)",
  "Brain eQTL (bulk cortex)"
)

cluster_order <- paste0("Cluster ", config$clusters_to_keep)

gene_order_df <- tibble(gene = genes_to_plot) %>%
  left_join(
    source_links %>%
      count(gene, name = "n_sources"),
    by = "gene"
  ) %>%
  left_join(
    cluster_links %>%
      count(gene, name = "n_clusters"),
    by = "gene"
  ) %>%
  mutate(
    n_sources = ifelse(is.na(n_sources), 0, n_sources),
    n_clusters = ifelse(is.na(n_clusters), 0, n_clusters)
  ) %>%
  arrange(
    desc(n_sources),
    desc(n_clusters),
    gene
  )

write_csv(gene_order_df, out_gene_order)

gene_order <- gene_order_df$gene
sector_order <- c(source_order, gene_order, cluster_order)

# UMAP-consistent cluster colours
umap_cluster_cols <- c(
  "Cluster 0" = "#D89000",
  "Cluster 1" = "#4DAFE8",
  "Cluster 2" = "#009A70",
  "Cluster 3" = "#E6D800",
  "Cluster 4" = "#006FB0",
  "Cluster 5" = "#D84A00",
  "Cluster 6" = "#C76AA3",
  "Cluster 7" = "#9acd84",
  "Cluster 8" = "#c599ef"
)

sector_cols <- c(
  "Blood eQTL (eQTLGen)" = "#4C78A8",
  "Immune cell eQTL (1k1k)" = "#F58518",
  "Brain cell eQTL (snRNA-seq)" = "#54A24B",
  "Brain eQTL (bulk cortex)" = "#B279A2",
  umap_cluster_cols[cluster_order]
)

gene_cols <- rep("#9E9E9E", length(gene_order))
names(gene_cols) <- gene_order

grid_col <- c(sector_cols, gene_cols)

source_colour_map <- c(
  "Blood eQTL (eQTLGen)" = "#4C78A880",
  "Immune cell eQTL (1k1k)" = "#F5851880",
  "Brain cell eQTL (snRNA-seq)" = "#54A24B80",
  "Brain eQTL (bulk cortex)" = "#B279A280"
)

cluster_colour_map <- c(
  "Cluster 0" = "#D8900080",
  "Cluster 1" = "#4DAFE880",
  "Cluster 2" = "#009A7080",
  "Cluster 3" = "#E6D80080",
  "Cluster 4" = "#006FB080",
  "Cluster 5" = "#D84A0080",
  "Cluster 6" = "#C76AA380",
  "Cluster 7" = "#9acd8480",
  "Cluster 8" = "#c599ef80"
)

link_cols <- c(
  source_colour_map[df_source$from],
  cluster_colour_map[df_cluster$to]
)


# ----------------------------- #
# 9. Chord diagram
# ----------------------------- #

if (nrow(plot_df) > 0 && length(gene_order) > 0) {
  
  pdf(
    out_chord_pdf,
    height = config$chord_pdf_height,
    width = config$chord_pdf_width
  )
  
  circos.clear()
  
  par(mar = c(4, 4, 4, 4))
  
  circos.par(
    start.degree = 90,
    gap.degree = 4,
    track.margin = c(0.01, 0.01),
    cell.padding = c(0, 0, 0, 0)
  )
  
  chordDiagram(
    x = plot_df,
    order = sector_order,
    grid.col = grid_col,
    col = link_cols,
    directional = 0,
    transparency = 0.15,
    big.gap = 12,
    small.gap = 2,
    annotationTrack = "grid",
    preAllocateTracks = list(
      track.height = max(strwidth(sector_order)) * 1.4
    )
  )
  
  circos.track(
    track.index = 1,
    panel.fun = function(x, y) {
      
      sector_name <- CELL_META$sector.index
      
      label <- if (sector_name == "Blood eQTL (eQTLGen)") {
        "Blood eQTL\n(eQTLGen)"
      } else if (sector_name == "Immune cell eQTL (1k1k)") {
        "Immune cell eQTL\n(1k1k)"
      } else if (sector_name == "Brain cell eQTL (snRNA-seq)") {
        "Brain cell eQTL\n(snRNA-seq)"
      } else if (sector_name == "Brain eQTL (bulk cortex)") {
        "Brain eQTL\n(bulk cortex)"
      } else {
        sector_name
      }
      
      circos.text(
        x = CELL_META$xcenter,
        y = CELL_META$ylim[1],
        labels = label,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = 1.1
      )
    },
    bg.border = NA
  )
  
  dev.off()
  circos.clear()
  
} else {
  
  warning(
    "No chord diagram was generated because no genes had both SMR source links ",
    "and cluster links after filtering."
  )
}


# ----------------------------- #
# 10. Summary tables
# ----------------------------- #

gene_source_summary <- source_links %>%
  count(gene, name = "n_sources") %>%
  left_join(
    source_links %>%
      group_by(gene) %>%
      summarise(
        sources = paste(sort(unique(source)), collapse = ", "),
        smr_regions = paste(sort(unique(region)), collapse = ", "),
        .groups = "drop"
      ),
    by = "gene"
  ) %>%
  left_join(
    cluster_links %>%
      group_by(gene) %>%
      summarise(
        clusters = paste(sort(unique(cluster_label)), collapse = ", "),
        .groups = "drop"
      ),
    by = "gene"
  ) %>%
  arrange(desc(n_sources), gene)

write_csv(gene_source_summary, out_gene_summary)

smr_gene_summary <- smr_all_sig %>%
  distinct(
    gene,
    source,
    p_FDR,
    p_HEIDI,
    chr,
    bp,
    region
  ) %>%
  arrange(
    region,
    gene,
    p_FDR
  )

write_csv(smr_gene_summary, out_smr_gene_summary)

smr_summary_counts <- smr_all_sig %>%
  summarise(
    n_smr_associations = n(),
    n_unique_genes = n_distinct(gene),
    n_mhc_associations = sum(region == "MHC_HEIDI_pass"),
    n_non_mhc_associations = sum(region == "nonMHC_HEIDI_pass"),
    n_mhc_genes = n_distinct(gene[region == "MHC_HEIDI_pass"]),
    n_non_mhc_genes = n_distinct(gene[region == "nonMHC_HEIDI_pass"])
  )

write_csv(smr_summary_counts, out_smr_summary_counts)

smr_supp_table <- smr_all_sig %>%
  arrange(
    region,
    source,
    gene,
    p_FDR
  ) %>%
  select(
    gene,
    source,
    chr,
    bp,
    region,
    p_FDR,
    p_HEIDI
  )

write_csv(smr_supp_table, out_smr_supplementary)

print(smr_summary_counts)


# ----------------------------- #
# 11. Excel workbook
# ----------------------------- #

wb <- createWorkbook()

addWorksheet(wb, "SMR_filtered")
writeData(wb, "SMR_filtered", smr_all_sig)

addWorksheet(wb, "Pseudobulk_filtered")
writeData(wb, "Pseudobulk_filtered", sc_sig)

addWorksheet(wb, "Overlap_genes")
writeData(wb, "Overlap_genes", tibble(gene = sort(overlap_genes)))

addWorksheet(wb, "Chord_edges")
writeData(wb, "Chord_edges", plot_df)

addWorksheet(wb, "Gene_order")
writeData(wb, "Gene_order", gene_order_df)

addWorksheet(wb, "Gene_source_cluster")
writeData(wb, "Gene_source_cluster", gene_source_summary)

addWorksheet(wb, "SMR_gene_summary")
writeData(wb, "SMR_gene_summary", smr_gene_summary)

addWorksheet(wb, "SMR_summary_counts")
writeData(wb, "SMR_summary_counts", smr_summary_counts)

addWorksheet(wb, "SMR_supp_table")
writeData(wb, "SMR_supp_table", smr_supp_table)

addWorksheet(wb, "Notes")
writeData(
  wb,
  "Notes",
  data.frame(
    item = c(
      "SMR data availability",
      "SMR filtering",
      "HEIDI filtering",
      "MHC annotation",
      "Pseudobulk filtering",
      "Chord diagram edges",
      "SMR source labels"
    ),
    note = c(
      "The SMR analysis data used as input for this workflow are available on Zenodo.org.",
      paste0("SMR associations are retained when p_FDR < ", config$smr_fdr_cutoff, "."),
      paste0("SMR associations are retained when p_HEIDI > ", config$heidi_cutoff, "."),
      paste0("MHC is annotated as chr", config$mhc_chr, ":", config$mhc_start, "-", config$mhc_end, "."),
      paste0("Stimulated latent pseudobulk DE genes are retained when p_val_adj < ", config$pseudobulk_fdr_cutoff, "."),
      "Chord edges link SMR source -> overlapping gene -> stimulated latent T-cell cluster.",
      paste(source_order, collapse = "; ")
    )
  )
)

saveWorkbook(
  wb,
  out_workbook,
  overwrite = TRUE
)


# ----------------------------- #
# 12. Session information
# ----------------------------- #

writeLines(
  c(
    "Configuration:",
    capture.output(str(config)),
    "",
    "SMR data availability:",
    "The SMR analysis data used as input for this workflow are available on Zenodo.org.",
    "",
    "SMR files:",
    unlist(lapply(smr_input_list, function(x) x$file_path)),
    "",
    "Pseudobulk input:",
    config$pseudobulk_file,
    "",
    "Pseudobulk cluster column used:",
    cluster_col,
    "",
    "SMR rows passing FDR + HEIDI:",
    as.character(nrow(smr_all_sig)),
    "",
    "SMR unique genes passing FDR + HEIDI:",
    as.character(length(unique(smr_all_sig$gene))),
    "",
    "Pseudobulk rows passing FDR:",
    as.character(nrow(sc_sig)),
    "",
    "Overlapping genes:",
    as.character(length(overlap_genes)),
    "",
    "Genes in chord diagram:",
    as.character(length(genes_to_plot)),
    "",
    "SMR summary counts:",
    capture.output(print(smr_summary_counts)),
    "",
    "Output files:",
    out_smr_all_filtered,
    out_smr_gene_summary,
    out_smr_summary_counts,
    out_smr_supplementary,
    out_pseudobulk_filtered,
    out_overlap_genes,
    out_chord_edges,
    out_gene_summary,
    out_gene_order,
    out_chord_pdf,
    out_workbook,
    "",
    "Session information:",
    capture.output(sessionInfo())
  ),
  out_session_info
)


# ----------------------------- #
# 13. Completion messages
# ----------------------------- #

message("Done.")
message("Output directory: ", config$output_dir)
message("SMR filtered table: ", out_smr_all_filtered)
message("Overlap genes: ", out_overlap_genes)
message("Chord edges: ", out_chord_edges)
message("Gene summary: ", out_gene_summary)
message("SMR supplementary table: ", out_smr_supplementary)
message("Workbook: ", out_workbook)
message("Chord PDF: ", out_chord_pdf)
message("Session info: ", out_session_info)
