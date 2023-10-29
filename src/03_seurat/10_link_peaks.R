# ------------------------------------------------------------------------------
# Link peaks to genes -- find peaks that are correlated with nearby genes
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(SeuratDisk)

# Load data and set seq style
dat <- readRDS(here::here("data/processed/seurat_object/09_dat_annot.rds"))
seqlevelsStyle(BSgenome.Dmelanogaster.UCSC.dm6) <- "NCBI"
seqlevelsStyle(dat@assays$peaks@annotation) <- "NCBI"

# Set assay
DefaultAssay(dat) <- "peaks"

# Compute the GC content for each peak
dat <- RegionStats(
  dat, 
  genome = BSgenome.Dmelanogaster.UCSC.dm6
)

# Link peaks to genes
dat <- LinkPeaks(
  object = dat,
  peak.assay = "peaks",
  expression.assay = "SCT"
)

# Save data
saveRDS(dat, here::here("data/processed/seurat_object/10_dat_linked.rds"))
SeuratDisk::as.loom(
  dat, 
  filename = here::here("data/processed/seurat_object/10_dat_linked.loom")
)
