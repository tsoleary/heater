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
dars <- FindMarkers(dat, 
                    ident.1 = "18°C",
                    ident.2 = "25°C")

# Filter only padj < 0.05
dars |> 
  filter(p_val_adj < 0.05)

# Save dars regardless of cluster
saveRDS(dars, here::here(out_dir, "dars.rds"))

# Differential expression between conditions within clusters -------------------

# Create combined groups of clusters and acclimation temperatures
dat$celltype_acc <- paste(dat$cell_type, dat$acc_temp, sep = "_")

# Set that new variable as the ident
Idents(dat) <- "celltype_acc"

# Find markers genes between different clusters 
dars <- list()

# Loop through each cluster
for (i in unique(dat$cell_type)) {
  dars[[i]] <- FindMarkers(dat, 
                           ident.1 = paste0(i, "_18°C"),
                           ident.2 = paste0(i, "_25°C")) |> 
    rownames_to_column("region")
}

# Combine all cell-types together
dars <- bind_rows(dars, .id = "cell_type")

# Rename cols for later and adjust p-vals again to account for all clusters 
# being tested at the same time groups
dars <- dars |> 
  dplyr::rename(avg_log2FC_18_25 = avg_log2FC,
         pct.18 = pct.1,
         pct.25 = pct.2) |> 
  mutate(padj = p.adjust(p_val_adj, method = "BH"))

# Save the dars cell-type
saveRDS(dars, here::here(out_dir, "dars_cell-type.rds"))


# Total number of DARs
dars |> 
  filter(padj < 0.05) |>
  tally()

# Number of DARs per cell-type
dars |>
  filter(padj < 0.05) |>
  group_by(cell_type) |>
  tally() |> 
  arrange(desc(n))

# If not used double adjusted p-values -------

# Total number of DARs
dars |> 
  filter(p_val_adj < 0.05) |>
  tally()

# Number of DARs per cell-type
dars |>
  filter(p_val_adj < 0.05) |>
  group_by(cell_type) |>
  tally() |> 
  arrange(desc(n))

