# ------------------------------------------------------------------------------
# Basic quality control metrics and filtering out non-cell barcodes from the 10X 
# Cell Ranger ARC called cells.
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)

# Load 10X filtered data
dat <- readRDS(
  here::here("data/processed/seurat_object/00_dat_10x_cells.rds")
)

# Output dir
out_dir <- "data/processed/seurat_object"

# Calculate Nucleosome Signal and TSS Enrichment for ATAC QC metrics
DefaultAssay(dat) <- "ATAC"
dat <- NucleosomeSignal(dat)
dat <- TSSEnrichment(dat)

# Calculate the percentage of RNA reads that map to mitochondrial genes 
DefaultAssay(dat) <- "RNA"
dat[["percent.mt"]] <- PercentageFeatureSet(
  dat,
  pattern = "^mt:"
)
# Calculate the percentage of RNA reads that map to ribosomal genes
dat[["percent.ribo"]] <- PercentageFeatureSet(
  dat, 
  pattern = "^Rp[S|L]"
)

# Save object with Nucleosome, TSS, and mitochondria meta.data added 
saveRDS(
  dat, 
  here::here(out_dir, "01_dat_10x_cells.rds")
)

# Filter out barcodes with low counts or high mitochondrial content ------------

# Define filtering thresholds
low_ATAC <- 800
low_RNA <- 200
max_mt <- 5
max_rb <- 25

# Subset 10x data based on thresholds
dat <- dat |>
  subset(
    nCount_RNA >= low_RNA &
      nCount_ATAC >= low_ATAC & 
      percent.mt < max_mt &
      percent.ribo < max_rb
    )

# Save object
saveRDS(
  dat, 
  here::here(out_dir, "01_dat_10x_filtered.rds")
)

# Create a filtered fragment file to be used in the Amulet doublet finder ------

# Path to fragment file with all barcodes included
data_dir <- "data/processed/seq/all/outs"
frag_path <- here::here(data_dir, "atac_fragments.tsv.gz")

# Write the file with fragments from only the filtered barcodes
FilterCells(
  frag_path,
  Cells(dat),
  outfile = here::here(data_dir, "atac_fragments_cells.tsv.gz")
)
