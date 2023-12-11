# ------------------------------------------------------------------------------
# Convert files to format for SCENIC+ analysis
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/14_dat_tf_motif.rds"))

# # Make SCT data for AnnData GEX
DefaultAssay(dat) <- "RNA"

# Save as an h5Seurat file
SeuratDisk::SaveH5Seurat(
  DietSeurat(dat, 
             assays = "RNA", 
             dimreducs = "umap"), 
  filename = here::here("data/processed/scenic/00_dat_gex.h5Seurat"),
  overwrite = TRUE
)

# Convert file to h5ad file
SeuratDisk::Convert(
  here::here("data/processed/scenic/00_dat_gex.h5Seurat"), 
  dest = "h5ad",
  overwrite = TRUE
)
