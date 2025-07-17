# ------------------------------------------------------------------------------
# Sample exclusion analysis of DARs
# Find differentially accessible regions between conditions
#   Subset out a 18°C sample to compare DAR lists
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/09_dat_annot.rds"))

# Set default assay
DefaultAssay(dat) <- "peaks"

# Output dir
out_dir <- "output/dars/sample_exclusion"

# Subset out samples to compare across replicates
# Remove 18°C Rep2
dat_18_Rep1 <- subset(dat, subset = sample_name != "18C_Rep2")
Idents(dat_18_Rep1) <- "acc_temp"
DefaultAssay(dat_18_Rep1) <- "peaks"

# Remove 18°C Rep1
dat_18_Rep2 <- subset(dat, subset = sample_name != "18C_Rep1")
Idents(dat_18_Rep2) <- "acc_temp"
DefaultAssay(dat_18_Rep2) <- "peaks"

# Test pseudobulk expression ---------------------------------------------------

## 18C Rep1 ------
dars_18_Rep1 <- FindMarkers(
  dat_18_Rep1, 
  ident.1 = "18°C",
  ident.2 = "25°C",
  test.use = "MAST",
  min.pct = 0,
  logfc.threshold = 0,
  min.cells.feature = 0
)

## 18C Rep2 ------
dars_18_Rep2 <- FindMarkers(
  dat_18_Rep2, 
  ident.1 = "18°C",
  ident.2 = "25°C",
  test.use = "MAST",
  min.pct = 0,
  logfc.threshold = 0,
  min.cells.feature = 0
)

# Save dars regardless of cluster
saveRDS(dars_18_Rep1, here::here(out_dir, "dars_bulk_18_Rep1.rds"))
saveRDS(dars_18_Rep2, here::here(out_dir, "dars_bulk_18_Rep2.rds"))


# Save all degs together for plotting later ----
# Combine all for comparison
dars_bulk <- full_join(
  dars, dars_18_Rep1,
  by = c("cell_type", "region"),
  suffix = c("", ".18.1")) |> 
  full_join(
    dars_18_Rep2,
    by = c("cell_type", "region"),
    suffix = c(".all", ".18.2")) |> 
  mutate(shared = case_when((
    avg_log2FC_18_25.18.1 > 0 & avg_log2FC_18_25.18.2 > 0) | 
      (avg_log2FC_18_25.18.1 < 0 & avg_log2FC_18_25.18.2 < 0) ~ TRUE,
    .default = FALSE))

saveRDS(dars_bulk, here::here(out_dir, "dars_bulk_MAST_exclusion.rds"))


# Cell-type specific differential expression between conditions ----------------

# Create combined groups of clusters and acclimation temperatures
dat_18_Rep1$celltype_acc <- paste(
  dat_18_Rep1$cell_type, 
  dat_18_Rep1$acc_temp, 
  sep = "_"
)

dat_18_Rep2$celltype_acc <- paste(
  dat_18_Rep2$cell_type, 
  dat_18_Rep2$acc_temp, 
  sep = "_"
)

# Set that new variable as the ident
Idents(dat_18_Rep1) <- "celltype_acc"
Idents(dat_18_Rep2) <- "celltype_acc"

## 18C Rep1 ------

# Find markers genes between different clusters 
dars <- list()

# Loop through each cluster for 18C Rep 1
for (i in unique(dat_18_Rep1$cell_type)) {
  dars[[i]] <- FindMarkers(
    dat_18_Rep1, 
    ident.1 = paste0(i, "_18°C"),
    ident.2 = paste0(i, "_25°C"),
    test.use = "MAST",
    min.pct = 0,
    logfc.threshold = 0,
    min.cells.feature = 0
  ) |>
    rownames_to_column("region")
}

# Combine all clusters together
dars <- bind_rows(dars, .id = "cell_type")

# Rename cols for later and adjust p-vals again to account for all clusters 
# being tested at the same time groups
dars <- dars |> 
  dplyr::rename("avg_log2FC_18_25" = "avg_log2FC",
                "pct.18" = "pct.1",
                "pct.25" = "pct.2")

# Save the dars cell types
saveRDS(dars, here::here(out_dir, "dars_cell-type_MAST_18_Rep1.rds"))


## 18C Rep2 ------

# Find markers genes between different clusters 
dars <- list()

# Loop through each cluster for 18C Rep 1
for (i in unique(dat_18_Rep2$cell_type)) {
  dars[[i]] <- FindMarkers(
    dat_18_Rep2, 
    ident.1 = paste0(i, "_18°C"),
    ident.2 = paste0(i, "_25°C"),
    test.use = "MAST",
    min.pct = 0,
    logfc.threshold = 0,
    min.cells.feature = 0
  ) |>
    rownames_to_column("region")
}

# Combine all clusters together
dars <- bind_rows(dars, .id = "cell_type")

# Rename cols for later and adjust p-vals again to account for all clusters 
# being tested at the same time groups
dars <- dars |> 
  dplyr::rename("avg_log2FC_18_25" = "avg_log2FC",
                "pct.18" = "pct.1",
                "pct.25" = "pct.2")

# Save the dars cell types
saveRDS(dars, here::here(out_dir, "dars_cell-type_MAST_18_Rep2.rds"))


# Save all dars together for plotting later ----
# Combine all for comparison
dars_all <- full_join(
  dars, dars_18_Rep1,
  by = c("cell_type", "region"),
  suffix = c("", ".18.1")) |> 
  full_join(
    dars_18_Rep2,
    by = c("cell_type", "region"),
    suffix = c(".all", ".18.2")) |> 
  mutate(shared = case_when((
    avg_log2FC_18_25.18.1 > 0 & avg_log2FC_18_25.18.2 > 0) | 
      (avg_log2FC_18_25.18.1 < 0 & avg_log2FC_18_25.18.2 < 0) ~ TRUE,
    .default = FALSE))

saveRDS(dars_all, here::here(out_dir, "dars_cell-type_MAST_exclusion.rds"))





