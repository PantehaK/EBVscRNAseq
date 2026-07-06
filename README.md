# Project Overview

This repository contains analysis scripts used to process and analyse single-cell RNA sequencing (scRNA-seq) data of EBV-specific CD8+ T cells from individuals with Multiple Sclerosis (MS) and non-MS controls.

The workflow processes 10x Genomics Cell Ranger outputs, performs quality control, integrates TCR and multimer metadata, performs clustering and transcriptional analysis, and analyses TCR repertoire features.

All analyses were performed in **R v4.5.0** using Seurat-based workflows.

# Repository Structure

The repository is organised into directories corresponding to major stages of the analysis pipeline.

```         
Scripts/
│
├── Quality_control/
├── Clustering/
├── Multimer_repertoire/
├── VDJ_analysis/
├── DEG_analysis/
├── MS_signature/
├── SMR_analysis/
│
├── Environment_container.R
├── renv/
├── renv.lock
└── sessionInfo.txt
```

Each directory contains scripts corresponding to a specific analysis stage.

# Reproducibility

The R environment used for this analysis is managed using **renv**.

The file `renv.lock` records the exact versions of R packages used in the project.

To recreate the computational environment:

``` r
install.packages("renv")
renv::restore()
```
# Analysis Workflow

The analysis is divided into four main stages:

1.  Quality control and metadata merging\
2.  Integration, Clustering and plotting\
3.  Multimer repertoire analysis with correlations\
4.  VDJ analysis\
5.  Differential gene expression analysis\
6.  MS signature phenotype\
7.  SMR analysis using relevant GWAS and eQTL datasets\

Scripts are organised into folders corresponding to each stage.

# 1. Quality Control

Directory: `Quality_control/`

These scripts process Cell Ranger outputs, perform QC filtering, detect doublets, and integrate VDJ and multimer metadata.

Scripts should be executed in the following order:

**1_Read_data.R**

-   Reads 10x filtered feature-barcode matrices
-   Creates initial Seurat objects

**2_QC_filtering.R**

-   Calculates QC metrics
-   Filters low-quality cells

**3_Doublet_removal.R**

-   Detects and removes doublets using DoubletFinder

**4_Cell_cycle_scoring.R**

-   Calculates cell cycle scores

**5_Tetramer_merging.R**

-   Integrates tetramer or multimer staining data

**6_VDJ_merge.R**

-   Integrates TCR VDJ sequencing data

**7_Merging_all_data.R**

-   Combines donor metadata with cell-level metadata

**8_Generate_TCR_combined_rds.R**

-   Creates integrated TCR and multimer metadata

**9_TCR_multimer_QC.R**

-   Performs final QC of TCR and multimer annotations

# 2. Clustering

Directory: `Clustering/`

These scripts perform integration, clustering, and module scoring.

**1_Harmony_integration_clustering.R**

-   Performs SCTransform normalisation
-   Runs Harmony integration
-   Performs dimensionality reduction and clustering
-   Performs module scoring for cytotoxicity-related genes

**2_cluster_composition.R**

-   Determines cluster composition
-   Performs relevant statistical analyses

**3_cluster_umaps.R**

-   Plots UMAPs for clusters and celltypes 

**4_feature_plots.R**

-   Plots feature plots of canonical cell type markers 

# 3. Multimer Repertoire Analysis

Directory: `Multimer_repertoire/`

These scripts analyse EBV-specific T cell frequency differences and correlates them with known EBV immunity markers.

**1_Estimate_tetramer_frequency.R**

-   Estimates frequencies of multimer-positive cells

**2_Spearman_correlation.R**

-   Performs spearman correlation test between estimated frequencies and markers of EBV burden

**3_HLA_stratified_response.R**

-   Plots estimated T cell frequencies based on protective/risk alleles in MS

**4_Cytotoxicity_score.R**

-   Generates median cytotoxicity scores per multimer-specific T cell per sample

# 4. T cell receptor repertoire analysis

Directory: `VDJ_analysis/`

These scripts analyse TCR diversity, expansion, clonality and relationships between antigen-specific T cells using relevant statistical testing.

**1_Expansion.R**

-   Estimates clonal expansion in antigen-specific T cells with statistical testing.

**2_D50_diversity.R**

-   Calculates D50 value as a measure of diversity and performs statistical testing.

**3_Clone_size.R**

-   Plots log-based clone size on UMAP of cells

**4_Clone_density.R**

-   Quantifies and compares clonal enrichment in different clusters.

**5_Clonal_sharing.R**

-   Quantifies clonal sharing in multimer-specific T cells between samples.

# 5. Differential Gene Expression Analysis

Directory: `DEG_analysis/`

These scripts define transcriptional signatures and perform differential expression analyses.

**1_Baseline_EBV_vs_other.R**

-   Defines baseline EBV-specific transcriptional signatures

**2_Baseline_lytic_vs_latent.R**

-   Defines baseline lytic and latent transcriptional signatures

**3_Baseline_MS_vs_control.R**

-   Performs pseudo-bulk differential expression analysis comparing MS and control samples in lytic, latent and non-EBV-specific T cells in each cluster.

**4_Stimulated_vs_baseline.R**

-   Defines transcriptional signatures for stimulated latent EBV-specific T cells

**5_Stimulated_MS_vs_control.R**

-   Performs pseudobulk differential expression analysis for stimulated latent-specific T cells in each cluster.

**6_hdWGCNA_pathway_enrichment.R**

-   Performs weighted gene co-expression network and pathway enrichment analysis on stimulated latent-specific T cells

# 6. MS Signature Phenotype

Directory: `MS_signature/`

These scripts define an MS-specific signature based on DEG analysis of stimulated latent-specific T cells.

**1_MS_signature_module_score.R**

-   Defines an MS signature and exports the mean signature of stimulated latent-specific T cells per sample with relevant statistical testing.

**2_MS_GWAS_gene_expression.R**

-   Exports expression levels of overlapping genes with MS GWAS (2019 IMSGC) enriched genes.

**3_Pheatmap_hierarchial_clustering.R**

-   Generates an unbiased hierarchial clustering heatmap visualising MS-associated gene signature across samples.

**4_GLM_CNS_injury.R**

-   Models the impact of MS signature and EBV-associated biomarkers on sGFAP and sNfL using GLM with covariates.
-   Plots forest plot with relevant statistical reporting.

# 7. SMR Analysis

Directory: `SMR_analysis/`

This script visualises the overlap between enriched genes from SMR-eQTL analyses of MS GWAS with differentially expressed genes post stimulated of latent-specific T cells using a chord diagram.

**1_SMR_analysis.R**

-   Performs SMR analysis on relevant GWAS and eQTL datasets.

**2_Chord_diagram.R**

-   Performs FDR and HEIDI filtering on SMR-anaylsis results
-   Plots a chord diagram

# Data Inputs

Raw sequencing data are not included in this repository due to file size constraints and confidentiality. FASTQ files as well as final metadata single-cell matrices will be available to download from EGA.

Required input data include:

-   10x Genomics filtered feature-barcode matrices
-   TCR VDJ contig annotations
-   supplementary_inputs_for_scripts available on Zenodo.org.

Input paths are defined within the scripts.

# Software

The analysis was performed in R v4.5.0 using various packages. Exact package versions are recorded in `renv.lock`.
