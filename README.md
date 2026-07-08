# Project Overview

This repository contains analysis scripts used to process and analyse single-cell RNA sequencing (scRNA-seq) data of EBV-specific CD8+ T cells from individuals with multiple sclerosis (MS) and non-MS controls.

The workflow processes 10x Genomics Cell Ranger outputs, performs quality control, integrates TCR and multimer metadata, performs clustering and transcriptional analysis, and analyses TCR repertoire features.

All analyses were performed in **R v4.5.0** using Seurat-based workflows inside a reproducible Apptainer/Singularity container.

# Repository Structure

The repository is organised into directories corresponding to major stages of the analysis pipeline.

```text
│
├── Quality_control/
├── Clustering/
├── Multimer_repertoire/
├── VDJ_analysis/
├── DEG_analysis/
├── MS_signature/
├── SMR_analysis/
├── Container/
    ├── scRNAseq_v4.def
    ├── EBV_scRNAseq_docker.sh
    ├── container_test_commands.md

renv.lock
README.Rmd
README.md
sessionInfo.txt
```

Each directory contains scripts corresponding to a specific analysis stage.

# Reproducibility

The primary reproducible analysis environment is provided as an Apptainer/Singularity container image.

The final container image is:

```text
EBV_snRNAseq_v4.sif
```

The image is available from Docker Hub as an OCI artifact:

```bash
apptainer pull EBV_snRNAseq_v4.sif oras://registry-1.docker.io/pantehak/ebv_snrnaseq_v4:v4
```

The container includes:

* R v4.5.0
* required R and Bioconductor packages
* system libraries required for package compilation and execution
* SMR

Scripts should be run through Apptainer using the required bind path. For example:

```bash
apptainer exec --cleanenv \
  --bind /data/menzies_projects/ms-epi:/data/menzies_projects/ms-epi \
  EBV_snRNAseq_v4.sif \
  Rscript Scripts/DEG_analysis/5_Stimulated_MS_vs_control.R
```

On systems where the project data are stored in a different location, users should modify the `--bind` path and the input paths within the scripts accordingly.

The container recipe and build script used to generate the image are provided in:

```text
container/scRNAseq_v4.def
container/EBV_scRNAseq_docker.sh
```

The `renv.lock` file is also included to document the R package versions used in the project. However, `renv.lock` alone is not sufficient to fully reproduce the environment because the analysis also depends on system libraries and command-line tools included in the container.

# Container Test

The Docker Hub-pulled image was tested using:

```bash
apptainer pull EBV_snRNAseq_v4_from_dockerhub.sif oras://registry-1.docker.io/pantehak/ebv_snrnaseq_v4:v4
```

The container was then launched with:

```bash
apptainer exec --cleanenv \
  --bind /data/menzies_projects/ms-epi:/data/menzies_projects/ms-epi \
  EBV_snRNAseq_v4_from_dockerhub.sif \
  Rscript -e "sessionInfo()"
```

This confirmed that the Docker Hub-pulled image launches successfully with **R v4.5.0**.

Example script execution:

```bash
apptainer exec --cleanenv \
  --bind /data/menzies_projects/ms-epi:/data/menzies_projects/ms-epi \
  EBV_snRNAseq_v4.sif \
  Rscript Scripts/MS_signature/3_Pheatmap_hierarchial_clustering.R
```

```bash
apptainer exec --cleanenv \
  --bind /data/menzies_projects/ms-epi:/data/menzies_projects/ms-epi \
  EBV_snRNAseq_v4.sif \
  Rscript Scripts/Clustering/4_feature_plots.R
```

```bash
apptainer exec --cleanenv \
  --bind /data/menzies_projects/ms-epi:/data/menzies_projects/ms-epi \
  EBV_snRNAseq_v4.sif \
  Rscript Scripts/DEG_analysis/5_Stimulated_MS_vs_control.R
```

# Files Required for Reproducibility

The `.sif` image itself is not stored directly in this GitHub repository due to file size constraints. Instead, it is available through Docker Hub using the pull command above.

This GitHub repository includes:

```text
renv.lock
container/scRNAseq_v4.def
container/EBV_scRNAseq_docker.sh
Scripts subfolders/
README.Rmd
README.md
sessionInfo.txt
```

The following file is hosted externally on Docker Hub:

```text
EBV_snRNAseq_v4.sif
```

Docker Hub / ORAS reference:

```text
oras://registry-1.docker.io/pantehak/ebv_snrnaseq_v4:v4
```

# Analysis Workflow

The analysis is divided into seven main stages:

1. Quality control and metadata merging
2. Integration, clustering and plotting
3. Multimer repertoire analysis with correlations
4. VDJ analysis
5. Differential gene expression analysis
6. MS signature phenotype analysis
7. SMR analysis using relevant GWAS and eQTL datasets

Scripts are organised into folders corresponding to each stage.

# 1. Quality Control

Directory: `Scripts/Quality_control/`

These scripts process Cell Ranger outputs, perform QC filtering, detect doublets, and integrate VDJ and multimer metadata.

Scripts should be executed in the following order:

## `1_Read_data.R`

* Reads 10x filtered feature-barcode matrices
* Creates initial Seurat objects

## `2_QC_filtering.R`

* Calculates QC metrics
* Filters low-quality cells

## `3_Doublet_removal.R`

* Detects and removes doublets using DoubletFinder

## `4_Cell_cycle_scoring.R`

* Calculates cell cycle scores

## `5_Tetramer_merging.R`

* Integrates tetramer or multimer staining data

## `6_VDJ_merge.R`

* Integrates TCR VDJ sequencing data

## `7_Merging_all_data.R`

* Combines donor metadata with cell-level metadata

## `8_Generate_TCR_combined_rds.R`

* Creates integrated TCR and multimer metadata

## `9_TCR_multimer_QC.R`

* Performs final QC of TCR and multimer annotations

# 2. Clustering

Directory: `Scripts/Clustering/`

These scripts perform integration, clustering, and module scoring.

## `1_Harmony_integration_clustering.R`

* Performs SCTransform normalisation
* Runs Harmony integration
* Performs dimensionality reduction and clustering
* Performs module scoring for cytotoxicity-related genes

## `2_cluster_composition.R`

* Determines cluster composition
* Performs relevant statistical analyses

## `3_cluster_umaps.R`

* Plots UMAPs for clusters and cell types

## `4_feature_plots.R`

* Plots feature plots of canonical cell type markers

# 3. Multimer Repertoire Analysis

Directory: `Scripts/Multimer_repertoire/`

These scripts analyse EBV-specific T cell frequency differences and correlate these with markers of EBV immunity.

## `1_Estimate_tetramer_frequency.R`

* Estimates frequencies of multimer-positive cells

## `2_Spearman_correlation.R`

* Performs Spearman correlation analyses between estimated frequencies and markers of EBV burden

## `3_HLA_stratified_response.R`

* Plots estimated T cell frequencies based on protective/risk alleles in MS

## `4_Cytotoxicity_score.R`

* Generates median cytotoxicity scores per multimer-specific T cell population per sample

# 4. T Cell Receptor Repertoire Analysis

Directory: `Scripts/VDJ_analysis/`

These scripts analyse TCR diversity, expansion, clonality, and relationships between antigen-specific T cells.

## `1_Expansion.R`

* Estimates clonal expansion in antigen-specific T cells with statistical testing

## `2_D50_diversity.R`

* Calculates D50 as a measure of diversity and performs statistical testing

## `3_Clone_size.R`

* Plots log-based clone size on UMAPs

## `4_Clone_density.R`

* Quantifies and compares clonal enrichment across clusters

## `5_Clonal_sharing.R`

* Quantifies clonal sharing of multimer-specific T cells between samples

# 5. Differential Gene Expression Analysis

Directory: `Scripts/DEG_analysis/`

These scripts define transcriptional signatures and perform differential expression analyses.

## `1_Baseline_EBV_vs_other.R`

* Defines baseline EBV-specific transcriptional signatures

## `2_Baseline_lytic_vs_latent.R`

* Defines baseline lytic and latent transcriptional signatures

## `3_Baseline_MS_vs_control.R`

* Performs pseudobulk differential expression analysis comparing MS and control samples in lytic, latent, and non-EBV-specific T cells within each cluster

## `4_Stimulated_vs_baseline.R`

* Defines transcriptional signatures for stimulated latent EBV-specific T cells

## `5_Stimulated_MS_vs_control.R`

* Performs pseudobulk differential expression analysis for stimulated latent-specific T cells in each cluster

## `6_hdWGCNA_pathway_enrichment.R`

* Performs weighted gene co-expression network analysis and pathway enrichment analysis on stimulated latent-specific T cells

# 6. MS Signature Phenotype

Directory: `Scripts/MS_signature/`

These scripts define and analyse an MS-specific transcriptional signature based on differential expression analysis of stimulated latent-specific T cells.

## `1_MS_signature_module_score.R`

* Defines an MS signature and exports the mean signature score of stimulated latent-specific T cells per sample with relevant statistical testing

## `2_MS_GWAS_gene_expression.R`

* Exports expression levels of overlapping genes between the MS signature and MS GWAS-enriched genes

## `3_Pheatmap_hierarchial_clustering.R`

* Generates an unbiased hierarchical clustering heatmap visualising the MS-associated gene signature across samples

## `4_GLM_CNS_injury.R`

* Models the impact of MS signature and EBV-associated biomarkers on sGFAP and sNfL using GLMs with covariates
* Plots forest plots with statistical reporting

# 7. SMR Analysis

Directory: `Scripts/SMR_analysis/`

These scripts analyse and visualise overlap between SMR-eQTL results and differentially expressed genes identified after stimulation of latent-specific T cells.

## `1_SMR_analysis.R`

* Performs SMR analysis on relevant GWAS and eQTL datasets

## `2_Chord_diagram.R`

* Performs FDR and HEIDI filtering on SMR results
* Plots a chord diagram

# Data Availability

Raw sequencing data are not included in this repository due to file size and data-access restrictions.

Required input data include:

* 10x Genomics filtered feature-barcode matrices
* TCR VDJ contig annotations
* sample metadata
* processed Seurat/RDS objects
* supplementary inputs required by the analysis scripts
* GWAS and eQTL summary statistics used for SMR analysis

Raw sequencing data and processed single-cell matrices will be made available through EGA. Supplementary script inputs will be made available through Zenodo.

Input paths are defined within the individual scripts and should be updated by users according to their local directory structure.

# Software

All analyses were performed using **R v4.5.0** inside the `EBV_snRNAseq_v4.sif` Apptainer/Singularity container.

The container was built using:

```text
container/scRNAseq_v4.def
container/EBV_scRNAseq_docker.sh
```

The R package versions are recorded in:

```text
renv.lock
```

The container can be pulled using:

```bash
apptainer pull EBV_snRNAseq_v4.sif oras://registry-1.docker.io/pantehak/ebv_snrnaseq_v4:v4
```

# Notes

This repository is intended to support reproducibility of the EBV-MS single-cell RNA-seq analysis. Users should update file paths in scripts to match their local data locations before running analyses.

For HPC systems, analyses should be run using Apptainer with the appropriate `--bind` path so that input and output directories are visible inside the container.
