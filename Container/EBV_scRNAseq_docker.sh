#!/bin/bash

### INFORMATION
### DATE: 2026-July-07
### AUTHOR: YUAN ZHOU

#############################################################################
#### Create a Docker image (for full reproducibility)
#############################################################################

##########
## STEP 2: create the docker def file: scRNAseq_v4.def
## this will build R + renv (restored R packages from renv.lock) +SMR

## on HPC interactive node , load apptanier
module load rosalind apptainer/ApptainerSigularityAlias

## then build the container
cd /data/menzies_projects/ms-epi/tmp/testYuan/docker

apptainer build EBV_snRNAseq_v4.sif scRNAseq_v4.def

## run the container
apptainer shell EBV_snRNAseq_v4.sif

#Inside:
# R with renv-restored packages
#R

# PLINK
#plink --help

# Python + LDSC
# ldsc --help
