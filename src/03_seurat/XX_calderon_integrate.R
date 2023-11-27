# ------------------------------------------------------------------------------
# Integrate data with Calderon -- Super quick for visualization
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(SeuratDisk)

# Load data and set seq style
dat <- readRDS(here::here("data/processed/seurat_object/10_dat_linked.rds"))
DefaultAssay(dat) <- "SCT"
cal <- readRDS(here::here("scratch/calderon/data/dat_5.rds"))

dat_list <- list(dat, cal)

features <- SelectIntegrationFeatures(object.list = dat_list)

dat_anchors <- FindIntegrationAnchors(object.list = dat_list, anchor.features = features)
# this command creates an 'integrated' data assay
dat_combined <- IntegrateData(anchorset = dat_anchors)


# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)