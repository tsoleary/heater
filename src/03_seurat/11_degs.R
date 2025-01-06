# ------------------------------------------------------------------------------
# Find differentially expressed genes between conditions
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/09_dat_annot.rds"))

# Output dir
out_dir <- "output/degs"

# Set default assay
DefaultAssay(dat) <- "SCT"

# Differential expression bulk between two conditions regardless of cluster ----
Idents(dat) <- "acc_temp"

degs <- FindMarkers(
  dat, 
  ident.1 = "18°C",
  ident.2 = "25°C",
  test.use = "MAST",
  min.pct = 0,
  logfc.threshold = 0,
  min.cells.feature = 0
)

degs |> 
  filter(abs(avg_log2FC) >= 0.25 &
           p_val_adj < 0.05) |> 
  filter(pct.1 >= 0.1 | pct.2 >= 0.1) |> 
  tally()

bind_rows(degs, degs_1, .id = "test") |> 
  as.data.frame(row.names = "gene") |> 
  pivot_wider(names_from = "gene")
  ggplot() +
  geom_histogram(aes(x = avg_log2FC)) +
  scale_y_continuous(trans = "log10") + facet_wrap(~test)

# Save degs regardless of cluster
saveRDS(degs, here::here(out_dir, "degs_bulk.rds"))

# Differential expression between conditions within clusters -------------------

# Create combined groups of clusters and acclimation temperatures
dat$celltype_acc <- paste(dat$cell_type, dat$acc_temp, sep = "_")

# Set that new variable as the ident
Idents(dat) <- "celltype_acc"

# Find markers genes between different clusters 
degs <- list()

# Loop through each cluster
for (i in unique(dat$cell_type)) {
  degs[[i]] <- FindMarkers(
    dat, 
    ident.1 = paste0(i, "_18°C"),
    ident.2 = paste0(i, "_25°C"),
    test.use = "MAST",
    min.pct = 0,
    logfc.threshold = 0,
    min.cells.feature = 0
    ) |>
    rownames_to_column("gene")
}

# Combine all clusters together
degs <- bind_rows(degs, .id = "cell_type")

# Rename cols for later and adjust p-vals again to account for all clusters 
# being tested at the same time groups
degs <- degs |> 
  dplyr::rename("avg_log2FC_18_25" = "avg_log2FC",
         "pct.18" = "pct.1",
         "pct.25" = "pct.2")

# Save the degs cell types
saveRDS(degs, here::here(out_dir, "degs_cell-type_MAST.rds"))

# Total number of DEGs
degs |>
  filter(abs(avg_log2FC_18_25) > 0.25 &
         p_val_adj < 0.05) |> 
  filter(pct.18 > 0.1 | pct.25 > 0.1) |> 
  tally()

# Number of DEGs per cell-type -----
degs |>
  group_by(cell_type) |> 
  filter(abs(avg_log2FC_18_25) > 0.25 &
           p_val_adj < 0.05) |> 
  filter(pct.18 > 0.1 | pct.25 > 0.1) |> 
  group_by(cell_type) |>
  tally()