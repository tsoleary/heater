# ------------------------------------------------------------------------------
# Find differentially expressed genes between conditions
# May 17, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/10_dat_linked.rds"))

# Set default assay
DefaultAssay(dat) <- "RNA"

# Differential expression bulk between two conditions regardless of cluster ----
Idents(dat) <- "acc_temp"
degs <- FindMarkers(dat, 
                    ident.1 = "18°C",
                    ident.2 = "25°C")

# Filter only padj < 0.05
degs %>%
  filter(p_val_adj < 0.05)

# Save degs regardless of cluster
saveRDS(degs, here::here("output/degs/degs.rds"))

# Differential expression between conditions within clusters -------------------

# Create combined groups of clusters and acclimation temperatures
dat$clust_acc <- paste(dat$seurat_clusters, dat$acc_temp, sep = "_")

# Set that new variable as the ident
Idents(dat) <- "clust_acc"

# Find markers genes between different clusters 
degs <- list()

# Loop through each cluster
for (i in levels(dat$seurat_clusters)) {
  degs[[i]] <- FindMarkers(dat, 
                      ident.1 = paste0(i, "_18°C"),
                      ident.2 = paste0(i, "_25°C")) %>%
    rownames_to_column("gene")
}

# Combine all clusters together
degs <- bind_rows(degs, .id = "cluster")

# Save the degs_clusters
saveRDS(degs, here::here("output/degs/degs_clusters.rds"))

# Total number of DEGs
degs %>%
  filter(p_val_adj < 0.05) %>%
  tally()

# Number of DEGs per cluster -----
degs %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  tally()
