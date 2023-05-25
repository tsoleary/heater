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
  here::here("data/processed/seurat_object/05_dat_filtered.rds")
)

# RNA data processing
DefaultAssay(dat) <- "RNA"
dat <- SCTransform(dat)
dat <- RunPCA(dat)

# ATAC data processing
DefaultAssay(dat) <- "ATAC"
dat <- FindTopFeatures(dat, min.cutoff = 5)
dat <- RunTFIDF(dat)
dat <- RunSVD(dat)

# The first lsi compenent often captures sequencing depth rather than biological 
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

# Annotate from reference
reference <- readRDS(here::here("scratch/calderon/data/dat_5.rds")) %>%
  SCTransform()

reference <- LoadH5Seurat(
  here::here("scratch/calderon/data/pbmc_multimodal.h5seurat")
)

# transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(
  reference = reference,
  query = dat,
  normalization.method = "SCT",
  reference.reduction = "spca",
  recompute.residuals = FALSE,
  dims = 1:50
)

predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = reference$celltype.l2,
  weight.reduction = pbmc[['pca']],
  dims = 1:50
)

pbmc <- AddMetaData(
  object = pbmc,
  metadata = predictions
)


# Save data
saveRDS(dat, here::here("data/processed/seurat_object/06_dat_cluster.rds"))
