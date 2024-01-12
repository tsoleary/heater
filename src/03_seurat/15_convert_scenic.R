# ------------------------------------------------------------------------------
# Convert files to format for SCENIC+ analysis
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/14_dat_tf_motif.rds"))

# # Make SCT data for AnnData GEX
DefaultAssay(dat) <- "SCT"

# Save as an h5Seurat file
SeuratDisk::SaveH5Seurat(
  DietSeurat(dat, 
             assays = "SCT", 
             dimreducs = "umap"), 
  filename = here::here("data/processed/scenic/gex.h5Seurat"),
  overwrite = TRUE
)

# Convert file to h5ad file
SeuratDisk::Convert(
  here::here("data/processed/scenic/gex.h5Seurat"), 
  dest = "h5ad",
  overwrite = TRUE
)

# Save meta data for cells
dat@meta.data |> 
  write.csv(here::here("data/processed/scenic/meta_data.csv"))


# Save the count matrix

  # Create peaks_matrix object
peaks_matrix <- dat@assays$peaks@counts

# Change region naming to chr:111-112 format
rownames(peaks_matrix) <- str_replace(
  rownames(peaks_matrix), 
  "mitochondrion-genome", 
  "mt_g") |> 
  str_replace("-", ":") |> 
  str_replace("mt_g", "mitochondrion-genome")

peaks_matrix |> 
  as.data.frame() |> 
  write.csv(here::here("data/processed/scenic/peaks_matrix.csv"))

# Save DARs in format that can be moved to python for SCENIC+ ------------------

dars <- readRDS(here::here("output/dars/dars.rds"))

dars_celltype <- readRDS(here::here("output/dars/dars_cell-type.rds"))


dars |> 
  dplyr::filter(p_val_adj < 0.05) |> 
  rownames_to_column("region") |> 
  mutate(region = str_replace(region, "-", ":")) |> 
  column_to_rownames("region") |> 
  write.csv(here::here("output/dars/dars.csv"))

dars_celltype |> 
  select(-padj) |> 
  dplyr::filter(p_val_adj < 0.05) |> 
  mutate(region = str_replace(region, "-", ":")) |> 
  write_csv(here::here("output/dars/dars_cell-type.csv"))


# Get the TSS sites
TSS_sites <- GetTSSPositions(
  dat@assays$ATAC@annotation
)

library(tidyverse)
ggplot(dat@meta.data,
  aes(
    y = nFeature_RNA,
    x = nCount_RNA
  )
) +
  geom_point()



