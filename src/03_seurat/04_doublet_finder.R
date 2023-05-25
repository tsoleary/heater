# ------------------------------------------------------------------------------
# Detecting multiplets in the multiome data sets
# DoubleFinder and amulet
# May 22, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(DoubletFinder)

# DoubletFinder for RNA doublets -----------------------------------------------

# Define expected doublet rate from 10X estimates
doublet_rate <- 0.008/1000
exp_double_rate <- doublet_rate*4000

# Each sample must be done separately per the instructions of DoubletFinder

# 18C_Rep1 ---------------------------------------------------------------------

# Load data
dat <- readRDS(
  here::here("data/processed/seurat_object/02_dat_10x_filtered.rds")) %>%
  subset(sample_name == "18C_Rep1")

# Pre-process Seurat object ----------------------------------------------------
DefaultAssay(dat) <- "RNA"
dat <- NormalizeData(dat)
dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
dat <- ScaleData(dat)
dat <- RunPCA(dat)
dat <- RunUMAP(dat, dims = 1:10)
dat <- FindNeighbors(dat, dims = 1:10)
dat <- FindClusters(dat, resolution = 0.1)

# pK Identification ------------------------------------------------------------
sweep.res.list_dat <- paramSweep_v3(dat, PCs = 1:10, sct = FALSE)
sweep.stats_dat <- summarizeSweep(sweep.res.list_dat, GT = FALSE)
bcmvn_dat <- find.pK(sweep.stats_dat)

# Homotypic Doublet Proportion Estimate ----------------------------------------
homotypic.prop <- modelHomotypic(dat@meta.data$seurat_clusters)
nExp_poi <- round(exp_double_rate * nrow(dat@meta.data))
nExp_poi.adj <- round(nExp_poi*(1 - homotypic.prop))

# Run DoubletFinder with varying classification stringencies -------------------
dat <- doubletFinder_v3(dat, 
                        PCs = 1:10,
                        pK = 0.23,
                        nExp = nExp_poi.adj,
                        reuse.pANN = FALSE, 
                        sct = FALSE)

# Save the classified doublets
doublets_18C_Rep1 <- dat@meta.data %>%
  subset(DF.classifications_0.25_0.23_61 == "Doublet") %>%
  rownames()

# 18C_Rep2 ---------------------------------------------------------------------

# Load data
dat <- readRDS(
  here::here("data/processed/seurat_object/02_dat_10x_filtered.rds")) %>%
  subset(sample_name == "18C_Rep2")

# Pre-process Seurat object ----------------------------------------------------
DefaultAssay(dat) <- "RNA"
dat <- NormalizeData(dat)
dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
dat <- ScaleData(dat)
dat <- RunPCA(dat)
dat <- RunUMAP(dat, dims = 1:10)
dat <- FindNeighbors(dat, dims = 1:10)
dat <- FindClusters(dat, resolution = 0.1)

# pK Identification ------------------------------------------------------------
sweep.res.list_dat <- paramSweep_v3(dat, PCs = 1:10, sct = FALSE)
sweep.stats_dat <- summarizeSweep(sweep.res.list_dat, GT = FALSE)
bcmvn_dat <- find.pK(sweep.stats_dat)

# Homotypic Doublet Proportion Estimate ----------------------------------------
homotypic.prop <- modelHomotypic(dat@meta.data$seurat_clusters)
nExp_poi <- round(exp_double_rate * nrow(dat@meta.data))
nExp_poi.adj <- round(nExp_poi*(1 - homotypic.prop))

# Run DoubletFinder with varying classification stringencies -------------------
dat <- doubletFinder_v3(dat, 
                        PCs = 1:10,
                        pK = 0.005,
                        nExp = nExp_poi.adj,
                        reuse.pANN = FALSE, 
                        sct = FALSE)

# Save the classified doublets
doublets_18C_Rep2 <- dat@meta.data %>%
  subset(DF.classifications_0.25_0.005_103 == "Doublet") %>%
  rownames()


# 25C_Rep1 ---------------------------------------------------------------------

# Load data
dat <- readRDS(
  here::here("data/processed/seurat_object/02_dat_10x_filtered.rds")) %>%
  subset(sample_name == "25C_Rep1")

# Pre-process Seurat object ----------------------------------------------------
DefaultAssay(dat) <- "RNA"
dat <- NormalizeData(dat)
dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
dat <- ScaleData(dat)
dat <- RunPCA(dat)
dat <- RunUMAP(dat, dims = 1:10)
dat <- FindNeighbors(dat, dims = 1:10)
dat <- FindClusters(dat, resolution = 0.1)

# pK Identification ------------------------------------------------------------
sweep.res.list_dat <- paramSweep_v3(dat, PCs = 1:10, sct = FALSE)
sweep.stats_dat <- summarizeSweep(sweep.res.list_dat, GT = FALSE)
bcmvn_dat <- find.pK(sweep.stats_dat)

# Homotypic Doublet Proportion Estimate ----------------------------------------
homotypic.prop <- modelHomotypic(dat@meta.data$seurat_clusters)
nExp_poi <- round(exp_double_rate * nrow(dat@meta.data))
nExp_poi.adj <- round(nExp_poi*(1 - homotypic.prop))

# Run DoubletFinder with varying classification stringencies -------------------
dat <- doubletFinder_v3(dat, 
                        PCs = 1:10,
                        pK = 0.3,
                        nExp = nExp_poi.adj,
                        reuse.pANN = FALSE, 
                        sct = FALSE)

# Save the classified doublets
doublets_25C_Rep1 <- dat@meta.data %>%
  subset(DF.classifications_0.25_0.3_63 == "Doublet") %>%
  rownames()

# 25C_Rep2 ---------------------------------------------------------------------

# Load data
dat <- readRDS(
  here::here("data/processed/seurat_object/02_dat_10x_filtered.rds")) %>%
  subset(sample_name == "25C_Rep2")

# Pre-process Seurat object ----------------------------------------------------
DefaultAssay(dat) <- "RNA"
dat <- NormalizeData(dat)
dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
dat <- ScaleData(dat)
dat <- RunPCA(dat)
dat <- RunUMAP(dat, dims = 1:10)
dat <- FindNeighbors(dat, dims = 1:10)
dat <- FindClusters(dat, resolution = 0.1)

# pK Identification ------------------------------------------------------------
sweep.res.list_dat <- paramSweep_v3(dat, PCs = 1:10, sct = FALSE)
sweep.stats_dat <- summarizeSweep(sweep.res.list_dat, GT = FALSE)
bcmvn_dat <- find.pK(sweep.stats_dat)

# Homotypic Doublet Proportion Estimate ----------------------------------------
homotypic.prop <- modelHomotypic(dat@meta.data$seurat_clusters)
nExp_poi <- round(exp_double_rate * nrow(dat@meta.data))
nExp_poi.adj <- round(nExp_poi*(1 - homotypic.prop))

# Run DoubletFinder with varying classification stringencies -------------------
dat <- doubletFinder_v3(dat, 
                        PCs = 1:10,
                        pK = 0.01,
                        nExp = nExp_poi.adj,
                        reuse.pANN = FALSE, 
                        sct = FALSE)

# Save the classified doublets
doublets_25C_Rep2 <- dat@meta.data %>%
  subset(DF.classifications_0.25_0.01_72 == "Doublet") %>%
  rownames()

doublets <- c(doublets_18C_Rep1,
              doublets_18C_Rep2,
              doublets_25C_Rep1,
              doublets_25C_Rep2)

# Save doublets
saveRDS(doublets, here::here("data/processed/qc/doublet_finder.rds"))
