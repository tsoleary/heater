# ------------------------------------------------------------------------------
# Dimension reduction and clustering
# May 30, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)
library(clustree)

# Load data
dat <- readRDS(
  here::here("data/processed/seurat_object/06_dat_qc.rds")
)

# Joint Clustering based on Signac tutorial url below --------------------------
# https://stuartlab.org/signac/articles/pbmc_multiomic.html

# RNA data processing
DefaultAssay(dat) <- "RNA"
dat <- SCTransform(dat)
dat <- RunPCA(dat)

# ATAC data processing
DefaultAssay(dat) <- "peaks"
dat <- FindTopFeatures(dat, min.cutoff = 5)
dat <- RunTFIDF(dat)
dat <- RunSVD(dat)

# The first lsi component often captures sequencing depth rather than biological 
# variation so we should not use this first component in further analysis.
# See plot below for ~-1 corr with depth and the first component
DepthCor(dat)

# Build a joint neighbor graph using both assays
dat <- FindMultiModalNeighbors(
  object = dat,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  verbose = TRUE
)

# Build a joint UMAP
dat <- RunUMAP(
  object = dat,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

# Build a joint t-SNE
dat <- RunTSNE(
  object = dat,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

# Find Clusters based on the LSI
dat <- FindClusters(dat, 
                    resolution = 0.2,
                    graph.name = "wknn", 
                    algorithm = 3)

# Save data
saveRDS(dat, here::here("data/processed/seurat_object/07_dat_cluster.rds"))


# RNA-only clustering ----------------------------------------------------------
# Load data
dat <- readRDS(
  here::here("data/processed/seurat_object/07_dat_cluster.rds")
)

# RNA data processing
DefaultAssay(dat) <- "RNA"
dat <- SCTransform(dat)
dat <- RunUMAP(RunPCA(dat), dims = c(1:50))

# Save dat with RNA only UMAP
saveRDS(dat,
  here::here("data/processed/seurat_object/07_dat_rna_umap.rds")
)

# ATAC-only clustering ---------------------------------------------------------
# Load data
dat <- readRDS(
  here::here("data/processed/seurat_object/07_dat_cluster.rds")
)
 
# ATAC UMAP quick
DefaultAssay(dat) <- "peaks"
dat <- dat |> 
  ScaleData() |>
  RunPCA() |>
  RunUMAP(dims = c(1:50))

# Save dat with ATAC only UMAP
saveRDS(dat,
  here::here("data/processed/seurat_object/07_dat_atac_umap.rds")
)

dat |> 
  DimPlot(
    split.by = "acc_temp")
