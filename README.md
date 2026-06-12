# Antigen-resolved remodeling of EBV-specific CD8 T-cell immunity in multiple sclerosis
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
├── Differential_gene_expression/
├── Multimer_repertoire/
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

1.  Quality control and metadata integration\
2.  Clustering and module scoring\
3.  Differential gene expression analysis\
4.  Multimer and TCR repertoire analysis

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

**6_Merging_VDJ_dataset.R**

-   Integrates TCR VDJ sequencing data

**7_Merging_all_metadata.R**

-   Combines donor metadata with cell-level metadata

**8_Generating_TCR_multimer_metadata.R**

-   Creates integrated TCR and multimer metadata

**9_QC_Multimer_TCR.R**

-   Performs final QC of TCR and multimer annotations

# 2. Clustering

Directory: `Clustering/`

These scripts perform integration, clustering, and module scoring.

**1_Harmony_integration_clustering.R**

-   Performs SCTransform normalisation
-   Runs Harmony integration
-   Performs dimensionality reduction and clustering
-   Relevant figures: 2A-F, 3A-D

**2_Cytotoxicity_module_scoring.R**

-   Calculates cytotoxicity-related transcriptional module scores
-   Relevant figures: 2G, S3

**3_MS_signature_module_score.R**

-   Calculates module scores for MS-associated transcriptional signatures
-   Relevant figures: 4D

# 3. Differential Gene Expression

Directory: `Differential_gene_expression/`

These scripts define transcriptional signatures and perform differential expression analyses.

**1_Baseline_EBV_signature.R**

-   Defines baseline EBV-specific transcriptional signatures
-   Relevant figures: S2A, S2B, Supplementary Table 2

**2_Pseudobulk_comparison_MS_vs_Control_baseline_EBV.R**

-   Performs pseudobulk differential expression analysis comparing MS and control samples
-   Relevant figures: Supplementary Table 3

**3_Stimulated_latent_T_signature.R**

-   Defines transcriptional signatures for stimulated latent EBV-specific T cells
-   Relevant figures: S4B, S5A, Supplementary Table 4

**4_Pseudobulk_comparison_MS_vs_Control_stimulated_latent.R**

-   Performs pseudobulk differential expression analysis for stimulated latent cells
-   Relevant figures: Supplementary Table 5

**5_Weighted_co_expression_analysis.R**

-   Performs weighted gene co-expression analysis
-   Relevant figures: Supplementary Table 6

**6_Pheatmap_MS_signature.R**

-   Generates heatmaps visualising MS-associated gene signatures
-   Relevant figures: 4C

# 4. Multimer and TCR Repertoire Analysis

Directory: `Multimer_repertoire/`

These scripts analyse EBV multimer binding and TCR repertoire characteristics.

**1_Estimating_multimer_frequency.R**

-   Estimates frequencies of multimer-positive cells
-   Relevant figures: 1B, 1C, 1D, Supplementary Table 10

**2_Clonal_sharing.R**

-   Identifies shared TCR clones
-   Relevant figures: 1E, S1A

**3_Clonal_diversity.R**

-   Calculates repertoire diversity metrics
-   Relevant figures: S1C

**4_Clonal_expansion.R**

-   Quantifies clonal expansion
-   Relevant figures: S1B

# Data Inputs

Raw sequencing data are not included in this repository due to file size constraints. FASTQ files as well as final metadata single-cell matrices will be available to download from https://zenodo.org/ 

Required input data include:

-   10x Genomics filtered feature-barcode matrices
-   TCR VDJ contig annotations
-   donor metadata
-   multimer metadata

Input paths are defined within the scripts.

# Software

The analysis was performed in R v4.5.0 using various packages. Exact package versions are recorded in `renv.lock`.
