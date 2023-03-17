# ------------------------------------------------------------------------------
# Early vs Early & Late vs Late Differential expression
# November 29, 2022
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(tidyverse)
set.seed(1)

# Load data
dat <- DietSeurat(readRDS(here::here("calderon/data/main.Rds")))
dat$runif <- runif(length(Cells(dat)))

# Create subsets
dat_1  <- subset(dat, NNv1_age >= 1 & NNv1_age < 2 & runif < 0.114)
dat_2  <- subset(dat, NNv1_age >= 2 & NNv1_age < 3 & runif < 0.32)
dat_17 <- subset(dat, NNv1_age >= 17 & NNv1_age < 18 & runif < 0.216)
dat_18 <- subset(dat, NNv1_age >= 18 & NNv1_age < 19 & runif < 0.5)

# Number of Cells in each downsampled dataset
c(length(Cells(dat_1)), length(Cells(dat_2)),
  length(Cells(dat_17)), length(Cells(dat_18)))

# Create list of subsetted data
dats_early <- list(dat_1, dat_2)
dats_late <- list(dat_17, dat_18)

# Perform integration ----------------------------------------------------------
dat_anchor_early <- FindIntegrationAnchors(object.list = dats_early, 
                                      dims = 1:20)
dat_early <- IntegrateData(anchorset = dat_anchor_early,
                          dims = 1:20)

# Dimensional reduction and UMAP projection for subsetted data

# Run the standard workflow for visualization and clustering
dat_early <- ScaleData(dat_early, verbose = FALSE)
dat_early <- RunPCA(dat_early, npcs = 30, verbose = FALSE)

# t-SNE and Clustering
dat_early <- RunUMAP(dat_early, reduction = "pca", dims = 1:20)
dat_early <- FindNeighbors(dat_early, reduction = "pca", dims = 1:20)
dat_early <- FindClusters(dat_early, resolution = 0.5)

# Perform integration ----------------------------------------------------------
dat_anchor_late <- FindIntegrationAnchors(object.list = dats_late, 
                                           dims = 1:20)
dat_late <- IntegrateData(anchorset = dat_anchor_late,
                           dims = 1:20)

# Dimensional reduction and UMAP projection for subsetted data

# Run the standard workflow for visualization and clustering
dat_late <- ScaleData(dat_late, verbose = FALSE)
dat_late <- RunPCA(dat_late, npcs = 30, verbose = FALSE)

# t-SNE and Clustering
dat_late <- RunUMAP(dat_late, reduction = "pca", dims = 1:20)
dat_late <- FindNeighbors(dat_late, reduction = "pca", dims = 1:20)
dat_late <- FindClusters(dat_late, resolution = 0.5)

# Find differentially expressed genes between hours ----------------------------
dat_early$time_window <- as.character(floor(dat_early$NNv1_age))

# Set identity to time window
Idents(dat_early) <- dat_early$time_window

# Identify DEGs
degs_early <- FindMarkers(dat_early, 
                          ident.1 = "1", 
                          ident.2 = "2",
                          logfc.threshold = 0.1,
                          min.pct = 0.02)

# Find differentially expressed genes between hours ----------------------------
dat_late$time_window <- as.character(floor(dat_late$NNv1_age))

# Set identity to time window
Idents(dat_late) <- dat_late$time_window

# Identify DEGs
degs_late <- FindMarkers(dat_late, 
                         ident.1 = "17", 
                         ident.2 = "18",
                         logfc.threshold = 0.1,
                         min.pct = 0.02)

degs <- bind_rows(list("early" = degs_early,
                       "late" = degs_late),
                       .id = "comparison")
degs
