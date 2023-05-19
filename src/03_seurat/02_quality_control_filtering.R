# ------------------------------------------------------------------------------
# Quality control and filtering of low-quality cells
# March 27, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Description -----
# Basic quality control metrics and filtering out non-cell barcodes.

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)

# Load 10X filtered data
dat_10x <- readRDS(
  here::here("data/processed/seurat_object/01_dat_10x_cells.rds")
)

# Calculate Nucleosome Signal and TSS Enrichment for ATAC QC metrics
DefaultAssay(dat_10x) <- "ATAC"
dat_10x <- NucleosomeSignal(dat_10x)
dat_10x <- TSSEnrichment(dat_10x)

# Calculate the percentage of RNA reads that map to mitochondrial genes 
DefaultAssay(dat_10x) <- "RNA"
dat_10x[["percent.mt"]] <- PercentageFeatureSet(dat_10x, 
                                                pattern = "^mt:")

# Save object with Nucleosome TSS, and mito metadata added 
saveRDS(dat_10x, "data/processed/seurat_object/02_dat_10x_cells.rds")

# FILTER OUT DOUBLETS ----------------------------------------------------------
# # Filter out 
# sce <- scds::cxds_bcds_hybrid(as.SingleCellExperiment(dat_10x))
# 
# singleCellTK::runScrublet(sce)
