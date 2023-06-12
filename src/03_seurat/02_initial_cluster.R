# ------------------------------------------------------------------------------
# Dimension reduction and clustering
# May 11, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)
library(clustree)

# Load data
dat <- readRDS(
  here::here("data/processed/seurat_object/01_dat_10x_filtered.rds")
)

# Run the standard workflow for visualization and clustering -------------------
# This initial clustering is just for the multiplets plots

# Set RNA to the default assay set for the following set of operations 
DefaultAssay(dat) <- "RNA"

# Log-Normalize data
dat <- NormalizeData(dat)

# Find variable features
dat <- FindVariableFeatures(dat)

# Scale and center data
dat <- ScaleData(dat)

# Run principle component analysis 
dat <- RunPCA(dat)

# Run UMAP
dat <- RunUMAP(dat, dims = 1:30)

# Construct nearest-neighbor graph
dat <- FindNeighbors(dat, dims = 1:30)

# Finding the correct resolution of cluster -----
# Run on a bunch of different resolutions
dat <- FindClusters(dat, resolution = seq(0.01, 0.1, 0.01))
# Look at the cluster tree
clustree::clustree(dat)

# Remove the RNA_snn columns added
dat@meta.data <- dat@meta.data |>
  select(!contains("RNA_snn"))

# Set the final clusters
dat <- FindClusters(dat, resolution = 0.03)

# Save data
saveRDS(
  dat, 
  here::here("data/processed/seurat_object/02_dat_initial_cluster.rds")
)
