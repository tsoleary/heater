# ------------------------------------------------------------------------------
# Find differentially expressed genes between conditions
# May 17, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/09_dat_annot.rds"))

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
dat$celltype_acc <- paste(dat$cell_type, dat$acc_temp, sep = "_")

# Set that new variable as the ident
Idents(dat) <- "celltype_acc"

# Find markers genes between different clusters 
degs <- list()

# Loop through each cluster
for (i in unique(dat$cell_type)) {
  degs[[i]] <- FindMarkers(dat, 
                      ident.1 = paste0(i, "_18°C"),
                      ident.2 = paste0(i, "_25°C")) %>%
    rownames_to_column("gene")
}

# Combine all clusters together
degs <- bind_rows(degs, .id = "cluster")

# Rename cols for later and adjust p-vals again to account for all clusters 
# being tested at the same time groups
degs <- degs |> 
  rename(avg_log2FC_18_25 = avg_log2FC,
         pct.18 = pct.1,
         pct.25 = pct.2) |> 
  mutate(padj = p.adjust(p_val_adj, method = "BH"))

# Save the degs_clusters
saveRDS(degs, here::here("output/degs/degs_clusters.rds"))

# Total number of DEGs
degs %>%
  filter(padj < 0.05) %>%
  tally()

# Number of DEGs per cluster -----
degs %>%
  filter(padj < 0.05) %>%
  group_by(cluster) %>%
  tally() |> 
  arrange(desc(n))

