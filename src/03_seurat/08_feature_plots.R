# ------------------------------------------------------------------------------
# Create feature plots
# May 16, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/03_dat_clustered.rds"))

# Gene
gene <- "sesB"

# Plot expression level on top of the UMAP projection
FeaturePlot(dat, 
            gene, 
            cols = c("grey95", "darkgreen"),
            keep.scale = "all",
            split.by = "acc_temp")

# Violin Plot
Idents(dat) <- "seurat_clusters"
VlnPlot(dat,
        gene,
        cols = c("#43aa8b", "#f3722c"),
        split.by = "acc_temp", 
        pt.size = 0)


Idents(dat) <- "clust_acc"
RidgePlot(dat, 
          gene)



