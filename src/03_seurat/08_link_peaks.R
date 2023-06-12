# ------------------------------------------------------------------------------
# Link peaks to genes
# June 06, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)
library(BSgenome.Dmelanogaster.UCSC.dm6)

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/07_dat_cluster.rds"))
seqlevelsStyle(BSgenome.Dmelanogaster.UCSC.dm6) <- "NCBI"

# Set assay
DefaultAssay(dat) <- "ATAC"

# Warning to me --------
# There is a problem with the peaks assay where a called peak from MACS2 CallPeaks
# give a peak that is overlapping beyond the end of 3L
# possible to remove my hand https://github.com/stuart-lab/signac/issues/1050


# Compute the GC content for each peak
dat <- RegionStats(
  dat, 
  genome = BSgenome.Dmelanogaster.UCSC.dm6
)

# Link peaks to genes
dat <- LinkPeaks(
  object = dat,
  peak.assay = "ATAC",
  expression.assay = "SCT"
)

# Save data
saveRDS(dat, here::here("data/processed/seurat_object/08_dat_linked.rds"))
