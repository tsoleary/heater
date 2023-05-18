# ------------------------------------------------------------------------------
# Find differentially expressed genes between conditions
# May 17, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(tidyverse)

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/03_dat_clustered.rds"))

# Set default assay
DefaultAssay(dat) <- "RNA"

# Reset primary idents to the acclimation temperatures
Idents(dat) <- "acc_temp"

# Average expression of RNA across all cells
avg_exp <- log10(AverageExpression(dat, verbose = FALSE)$RNA)
avg_exp <- avg_exp %>%
  as_tibble(rownames = "gene")

# Scatter plot of the average expression for all cells -------------------------
p <- avg_exp %>%
  ggplot() +
  geom_point(aes(x = `25C`,
                 y = `18C`,
                 label = gene)) +
  labs(x = "18°C",
       y = "25°C") +
  theme_minimal()

# Quick plotly with gene labels
plotly::ggplotly(p)

# Differential expression bulk between two conditions regardless of cluster ----
degs <- FindMarkers(dat, 
                    ident.1 = "18C",
                    ident.2 = "25C")

degs %>%
  filter(p_val_adj < 0.05)

# I---
# AggregateExpression(dat, assay = "RNA")

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
                      ident.1 = paste0(i, "_18C"),
                      ident.2 = paste0(i, "_25C")) %>%
    rownames_to_column("gene")
}

# Combine all clusters together
degs <- bind_rows(degs, .id = "cluster")

# Total number of DEGs
degs %>%
  filter(p_val_adj < 0.05) %>%
  tally()

# Number of DEGs per cluster -----
degs %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  tally()

# Filter to just the genes with an adjusted p-value less than 0.05 -------------
degs <- degs %>%
  filter(p_val_adj < 0.05)
