# ------------------------------------------------------------------------------
# Differential expression and accessibility
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/07_dat_cluster.rds"))

# Output dir
out_dir <- "output/markers"

# Set resolution for the seurat_clusters
dat@meta.data$seurat_clusters <- dat@meta.data$wknn_res.0.7
Idents(dat) <- "seurat_clusters"

# Find markers of each cluster that are conserved between the two conditions ---

# Create empty list to populate in the for loop
markers <- list()

# Loop through each cluster individually
for (i in levels(dat)) {
  markers[[i]] <- FindConservedMarkers(dat, 
                                       ident.1 = i,
                                       grouping.var = "acc_temp",
                                       only.pos = TRUE) |> 
    rownames_to_column("gene")
}

# Combine markers for all clusters
markers <- bind_rows(markers, 
                     .id = "cluster")

# Number of markers per cluster
markers |> 
  mutate(cluster = as.numeric(cluster)) |> 
  filter(max_pval < 0.05) |> 
  group_by(cluster) |> 
  tally()

# Save markers object
saveRDS(markers, here::here(out_dir, "cluster_markers.rds"))
