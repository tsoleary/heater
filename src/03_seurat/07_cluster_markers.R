# ------------------------------------------------------------------------------
# Differential expression and accessibility
# May 30, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/06_dat_cluster.rds"))

# Find markers of each cluster that are conserved between the two conditions ---

# Create empty list to populate in the for loop
markers <- list()

# Loop through each cluster individually
for (i in levels(dat)) {
  markers[[i]] <- FindConservedMarkers(dat, 
                                       ident.1 = i,
                                       grouping.var = "acc_temp",
                                       only.pos = TRUE) %>%
    rownames_to_column("gene")
}

# Combine markers for all clusters
markers <- bind_rows(markers, 
                     .id = "cluster")

# Number of markers per cluster
markers %>%
  filter(max_pval < 0.05) %>%
  group_by(cluster) %>%
  tally()

# Save markers object
saveRDS(markers, here::here("data/processed/genes/markers.rds"))

# # Trying to annotate the clusters based on their cell type ---------------------
# 
# # Load annotations from BDGP
# annot <- read_csv(here::here("data/raw/annot/insitu_annot.csv"))
# 
# # Top 10 markers per cluster
# top_markers <- markers %>%
#   filter(max_pval < 10^-4) %>%
#   group_by(cluster) %>%
#   slice_min(max_pval, n = 10, with_ties = FALSE)
# 
# # This doesn't work -- lol
# x <- top_markers %>%
#   left_join(annot) %>%
#   filter(!is.na(annot)) %>%
#   select(cluster, gene, annot) %>%
#   group_by(cluster) %>%
#   distinct() %>%
#   reframe(annots = paste(annot, sep = "_"))
