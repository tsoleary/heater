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

# Joint Clustering based on Signac tutorial ------------------------------------
# https://stuartlab.org/signac/articles/pbmc_multiomic.html

# RNA data processing
DefaultAssay(dat) <- "RNA"
dat <- SCTransform(dat)
dat <- RunPCA(dat)

<<<<<<< HEAD
# Quick elbow plot of the first 50 PCs
ElbowPlot(dat, ndims = 50) +
  cowplot::theme_minimal_grid()

ggsave(here::here("output/figs/qc/elbow_plot.png"),
       height = 20,
       width = 25,
       units = "cm")
ggsave(here::here("output/figs/qc/elbow_plot.pdf"),
       height = 20,
       width = 25,
       units = "cm")

# Keeping all 50 PCs for a few reasons -----
# The elbow plot does not seem to entirely tail off before 50 PCs and the 
# SCTransform is a more robust normalization method that is less likely to 
# carry artifacts do to technical variation. Please see the below two links
# https://satijalab.org/seurat/archive/v3.0/sctransform_vignette.html
# https://www.biostars.org/p/423306/

=======
>>>>>>> a8de26aef8a500cdbd82fa277045d3e778c75e21
# ATAC data processing
DefaultAssay(dat) <- "peaks"
dat <- FindTopFeatures(dat, min.cutoff = 5)
dat <- RunTFIDF(dat)
dat <- RunSVD(dat)

# The first lsi component often captures sequencing depth rather than biological 
# variation so we should not use this first component in further analysis.
<<<<<<< HEAD
# See plot below for ~1 corr with depth and the first component
DepthCor(dat)

ggsave(here::here("output/figs/qc/lsi_depthcor.png"),
       height = 20,
       width = 25,
       units = "cm")
ggsave(here::here("output/figs/qc/lsi_depthcor.pdf"),
       height = 20,
       width = 25,
       units = "cm")

=======
# See plot below for ~-1 corr with depth and the first component
DepthCor(dat)

>>>>>>> a8de26aef8a500cdbd82fa277045d3e778c75e21
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

<<<<<<< HEAD
# Find Clusters on the wknn
dat <- FindClusters(
  dat, 
  resolution = seq(from = 0.1, to = 1, by = 0.1),
  graph.name = "wknn", 
  algorithm = 3
)
=======
# Find Clusters based on the LSI
dat <- FindClusters(dat, 
                    resolution = 0.2,
                    graph.name = "wknn", 
                    algorithm = 3)
>>>>>>> a8de26aef8a500cdbd82fa277045d3e778c75e21

# Save data
saveRDS(dat, here::here("data/processed/seurat_object/07_dat_cluster.rds"))

################################################################################
# To see how and if the dimensional reduction of only ATAC or RNA libraries 
# looks different -- save these data below

# RNA-only dim-reduction -------------------------------------------------------
# Load data
dat <- readRDS(
  here::here("data/processed/seurat_object/07_dat_cluster.rds")
)

<<<<<<< HEAD
# Write over UMAP reduction with RNA-only UMAP
dat <- RunUMAP(
  dat, 
  assay = "SCT",
  dims = c(1:50)
)
# Write over t-SNE reduction with RNA-only t-SNE
dat <- RunTSNE(
  dat, 
  assay = "SCT",
  dims = c(1:50)
)

# Save dat with RNA only UMAP
saveRDS(dat,
  here::here("data/processed/seurat_object/07_dat_rna_only_dim.rds")
=======
# RNA data processing
DefaultAssay(dat) <- "RNA"
dat <- SCTransform(dat)
dat <- RunUMAP(RunPCA(dat), dims = c(1:50))

# Save dat with RNA only UMAP
saveRDS(dat,
  here::here("data/processed/seurat_object/07_dat_rna_umap.rds")
>>>>>>> a8de26aef8a500cdbd82fa277045d3e778c75e21
)

# ATAC-only dim-reduction ------------------------------------------------------
# Load data
dat <- readRDS(
  here::here("data/processed/seurat_object/07_dat_cluster.rds")
)
<<<<<<< HEAD

DefaultAssay(dat) <- "peaks"
dat <- ScaleData(dat)
dat <- RunPCA(dat)
dat <- RunUMAP(dat, dims = c(1:50))
dat <- RunTSNE(dat, dims = c(1:50))

# Save dat with ATAC-only UMAP
saveRDS(dat,
  here::here("data/processed/seurat_object/07_dat_atac_only_dim.rds")
=======
 
# ATAC UMAP quick
DefaultAssay(dat) <- "peaks"
dat <- dat |> 
  ScaleData() |>
  RunPCA() |>
  RunUMAP(dims = c(1:50))

# Save dat with ATAC only UMAP
saveRDS(dat,
  here::here("data/processed/seurat_object/07_dat_atac_umap.rds")
>>>>>>> a8de26aef8a500cdbd82fa277045d3e778c75e21
)
