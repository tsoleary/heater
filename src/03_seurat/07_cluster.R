# ------------------------------------------------------------------------------
# Dimension reduction and clustering
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)

# Load data --------------------------------------------------------------------
dat <- readRDS(
  here::here("data/processed/seurat_object/06_dat_qc.rds")
)
# Save S phase genes in a vector
s_genes <- read_csv(
  here::here("data/raw/annot/dmel_cell-cycle_genes.csv")) |> 
  filter(phase == "S") |> 
  select(gene) |> 
  deframe()
# Save G2 or M phase genes in a vector
g2m_genes <- read_csv(
  here::here("data/raw/annot/dmel_cell-cycle_genes.csv")) |> 
  filter(phase == "G2/M") |> 
  select(gene) |> 
  deframe()

# Output dir
out_dir <- "data/processed/seurat_object"

# Joint Clustering based on Signac tutorial ------------------------------------
# https://stuartlab.org/signac/articles/pbmc_multiomic.html
# https://satijalab.org/seurat/archive/v3.0/cell_cycle_vignette.html

# RNA data processing -----

# RNA data pre-processing before cell-cycle scoring
DefaultAssay(dat) <- "RNA"
dat <- SCTransform(dat)

# Score cell cycle
dat <- CellCycleScoring(
  dat,
  s.features = s_genes,
  g2m.features = g2m_genes,
  set.ident = TRUE
)

# Set assay back to RNA and re-run with S.Score and G2M.Score regressed out 
DefaultAssay(dat) <- "RNA"
dat <- SCTransform(
  dat, 
  vars.to.regress = c("S.Score", "G2M.Score"),
  variable.features.n = 3000 # Deafults to 3k, should we change to less?
)
dat <- RunPCA(dat)

# ATAC data processing -----
DefaultAssay(dat) <- "peaks"
dat <- FindTopFeatures(dat, min.cutoff = 5)
dat <- RunTFIDF(dat)
dat <- RunSVD(dat)

# Build a joint neighbor graph using both assays -----
dat <- FindMultiModalNeighbors(
  object = dat,
  k.nn = 100,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  verbose = TRUE
)

# Build a joint UMAP
dat <- RunUMAP(
  object = dat,
  n.neighbors = 100L,
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

# Find clusters on the wknn at resolutions ranging from 0.1 to 1
dat <- FindClusters(
  dat, 
  resolution = seq(from = 0.1, to = 1, by = 0.1),
  graph.name = "wknn", 
  algorithm = 3
)

# Save data
saveRDS(dat, here::here(out_dir, "07_dat_cluster.rds"))


################################################################################
# Creating RNA-only and ATAC-only dimension reductions -------------------------
# To see how and if the dimensional reduction look different 

# RNA-only dim-reduction -------------------------------------------------------
# Load data
dat <- readRDS(
  here::here("data/processed/seurat_object/07_dat_cluster.rds")
)

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
  here::here(out_dir, "07_dat_rna_only_dim.rds")
)

# ATAC-only dim-reduction ------------------------------------------------------
# Load data
dat <- readRDS(
  here::here("data/processed/seurat_object/07_dat_cluster.rds")
)

DefaultAssay(dat) <- "peaks"
dat <- ScaleData(dat)
dat <- RunPCA(dat)
dat <- RunUMAP(dat, dims = c(1:50))
dat <- RunTSNE(dat, dims = c(1:50))

# Save dat with ATAC-only UMAP
saveRDS(dat,
  here::here(out_dir, "07_dat_atac_only_dim.rds")
)
