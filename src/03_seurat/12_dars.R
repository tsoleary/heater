# ------------------------------------------------------------------------------
# Find differentially accessible regions within cell-types between acclimation
# temperatures
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/09_dat_annot.rds"))

# Output dir
out_dir <- "output/dars"

# Set default assay
DefaultAssay(dat) <- "peaks"

# Differential accessibility bulk between two conditions regardless of cluster -
Idents(dat) <- "acc_temp"
dars <- FindMarkers(
  dat, 
  ident.1 = "18°C",
  ident.2 = "25°C",
  test.use = "MAST",
  min.pct = 0,
  logfc.threshold = 0,
  min.cells.feature = 0
)

# Filter only padj < 0.05
dars |> 
  filter(abs(avg_log2FC) >= 0.25 &
           p_val_adj < 0.05) |> 
  filter(pct.1 >= 0.1 | pct.2 >= 0.1) |> 
  tally()

# Save dars regardless of cluster
saveRDS(dars, here::here(out_dir, "dars_bulk_MAST.rds"))

# Differential expression between conditions within clusters -------------------

# Create combined groups of clusters and acclimation temperatures
dat$celltype_acc <- paste(dat$cell_type, dat$acc_temp, sep = "_")

# Set that new variable as the ident
Idents(dat) <- "celltype_acc"

# Find markers genes between different clusters 
dars <- list()

# Loop through each cluster
for (i in unique(dat$cell_type)) {
  dars[[i]] <- FindMarkers(
    dat, 
    ident.1 = paste0(i, "_18°C"),
    ident.2 = paste0(i, "_25°C"),
    test.use = "MAST",
    min.pct = 0,
    logfc.threshold = 0,
    min.cells.feature = 0) |> 
    rownames_to_column("region")
}

# Combine all cell-types together
dars <- bind_rows(dars, .id = "cell_type")

# Rename cols for later and adjust p-vals again to account for all clusters 
# being tested at the same time groups
dars <- dars |> 
  dplyr::rename(
    avg_log2FC_18_25 = avg_log2FC,
    pct.18 = pct.1,
    pct.25 = pct.2)

# Save the dars cell-type
saveRDS(dars, here::here(out_dir, "dars_cell-type_MAST.rds"))

# Total number of DARs
dars |> 
  filter(abs(avg_log2FC_18_25) >= 0.25 &
           p_val_adj < 0.05) |> 
  filter(pct.18 >= 0.1 | pct.25 >= 0.1) |> 
  tally()

# Number of DARs per cell-type
dars |>
  filter(abs(avg_log2FC_18_25) >= 0.25 &
           p_val_adj < 0.05) |> 
  filter(pct.18 >= 0.1 | pct.25 >= 0.1) |> 
  group_by(cell_type) |>
  tally() |> 
  arrange(desc(n))

