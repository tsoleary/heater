# ------------------------------------------------------------------------------
# Quality control and filtering of low-quality cells
# March 27, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Description -----
# Basic quality control metrics and filtering out non-cell barcodes from the 10X 
# Cell Ranger ARC called cells.

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)

# Load 10X filtered data
dat <- readRDS(
  here::here("data/processed/seurat_object/00_dat_10x_cells.rds")
)

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
<<<<<<< HEAD
  dat, 
=======
  data, 
>>>>>>> a8de26aef8a500cdbd82fa277045d3e778c75e21
  pattern = "^Rp[S|L]"
)

# Save object with Nucleosome, TSS, and mitochondria meta.data added 
<<<<<<< HEAD
saveRDS(
  dat, 
  here::here("data/processed/seurat_object/01_dat_10x_cells.rds")
)
=======
saveRDS(dat, "data/processed/seurat_object/01_dat_10x_cells.rds")
>>>>>>> a8de26aef8a500cdbd82fa277045d3e778c75e21

# Filter out barcodes with low counts or high mitochondrial content ------------

# Define filtering thresholds
low_ATAC <- 800
low_RNA <- 200
max_mt <- 5 #### Calderon did 25 % for both mt and rb -- ask Seth? Maybe
<<<<<<< HEAD
max_rb <- 5
=======
#max_rb <- 5
>>>>>>> a8de26aef8a500cdbd82fa277045d3e778c75e21

# Subset 10x data based on thresholds
dat <- dat |>
  subset(
    nCount_RNA >= low_RNA &
      nCount_ATAC >= low_ATAC & 
<<<<<<< HEAD
      percent.mt < max_mt &
      percent.ribo < max_rb
    )

# Save object
saveRDS(
  dat, 
  here::here("data/processed/seurat_object/01_dat_10x_filtered.rds")
)
=======
      percent.mt < max_mt # &
      # percent.ribo < max_rb
    )

# Save object
saveRDS(dat, "data/processed/seurat_object/01_dat_10x_filtered.rds")
>>>>>>> a8de26aef8a500cdbd82fa277045d3e778c75e21

# Create a filtered fragment file to be used in the amulet doublet finder ------

# Path to fragment file with all barcodes included
data_dir <- here::here("data/processed/seq/all/outs/")
frag_path <- paste0(data_dir, "atac_fragments.tsv.gz")

# Write the file with fragments from only the filtered barcodes
FilterCells(frag_path,
            Cells(dat),
            outfile = paste0(data_dir, "atac_fragments_cells.tsv.gz"))
