# ------------------------------------------------------------------------------
# Dimension reduction and clustering
# May 11, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(tidyverse)

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/02_dat_cells.rds"))

# Analyze data ----

# Run the standard workflow for visualization and clustering
dat <- ScaleData(dat, verbose = FALSE)
dat <- RunPCA(dat, npcs = 30, verbose = FALSE)

# t-SNE and Clustering
dat <- RunUMAP(dat, reduction = "pca", dims = 1:20)
dat <- FindNeighbors(dat, reduction = "pca", dims = 1:20)
dat <- FindClusters(dat, resolution = 0.5)

# Plot 
DimPlot(dat,
        group.by = "acc_temp")
