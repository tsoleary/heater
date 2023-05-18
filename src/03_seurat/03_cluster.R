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
dat <- readRDS(here::here("data/processed/seurat_object/01_dat_10x_cells.rds"))

# Run the standard workflow for visualization and clustering -------------------

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

# Finding the correct resolution of clusters
clustree::clustree(dat)

# Determine the clusters within the data
dat <- FindClusters(dat, resolution = 0.1)

# Plot UMAP
DimPlot(dat,
        split.by = "acc_temp")

# Count per cluster
dat@meta.data %>%
  group_by(seurat_clusters) %>%
  tally()

# Count per cluster per acclimation state
dat@meta.data %>%
  group_by(seurat_clusters, acc_temp) %>%
  tally() %>%
  pivot_wider(values_from = n, names_from = acc_temp)

# Quick bar plot counting the number of cells for each cluster
dat@meta.data %>%
  ggplot(aes(x = seurat_clusters,
             fill = acc_temp)) +
  geom_bar(position = "dodge",
           color = "grey50") +
  labs(y = "Number of cells",
       x = "Cluster") +
  scale_fill_manual(name = element_blank(),
                     values = c("#43aa8b", "#f3722c")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  cowplot::theme_minimal_hgrid()

# Log-Normalize ATAC Reads ------
dat <- NormalizeData(dat, assay = "ATAC")

# Save data
saveRDS(dat, here::here("data/processed/seurat_object/03_dat_clustered.rds"))

