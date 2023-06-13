# ------------------------------------------------------------------------------
# Quality control and filtering of low-quality cells
# May 25, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Description -----
# Quality control metrics of final cells

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)
library(AnnotationHub)

# Load 10X filtered data
dat <- readRDS(
  here::here("data/processed/seurat_object/05_dat_filtered.rds")
)

# Load path to fragment file for later counting
data_dir <- here::here("data/processed/seq/all/outs/")
frag_path <- paste0(data_dir, "atac_fragments.tsv.gz")

# Get gene annotations for dm6 -------------------------------------------------
BDGP6.32 <- query(AnnotationHub(), 
                  c("EnsDb", "Drosophila melanogaster", "109"))[[1]]
annotation <- GetGRangesFromEnsDb(BDGP6.32)
seqlevelsStyle(annotation) <- 'UCSC'

# Calculate additional QC metrics ----------------------------------------------

# Calculate TSS Enrichment matrix for the cells that are left
DefaultAssay(dat) <- "ATAC"
dat <- TSSEnrichment(dat, fast = FALSE)

# Counting fraction of reads in peaks -----

# Count fragments
total_fragments <- CountFragments(
  fragments = frag_path, 
  cells = Cells(dat)
)

# Add to metadata
dat$fragments <- total_fragments[total_fragments$CB == colnames(dat), 
                                 "frequency_count"]

# Peak calling -----------------------------------------------------------------
# Installed macs2 using the R package Herper

# Path to miniconda
path_to_miniconda <-
  "/slipstream_old/home/thomasoleary/.local/share/r-miniconda" 
# Path to macs2 for CallPeaks function
macs2_path <- paste0(path_to_miniconda, 
                     "/envs/PeakCalling_analysis/bin/macs2")

# Call peaks
peaks <- CallPeaks(
  dat, 
  macs2.path = macs2_path,
  effective.genome.size = 1.25e8
)

# Create the peak matrix
peak_matrix <- FeatureMatrix(
  fragments = Fragments(dat),
  features = peaks
)

# Create a new assay using the MACS2 peak set and add it to the Seurat object
dat[["peaks"]] <- CreateChromatinAssay(
  counts = peak_matrix,
  fragments = frag_path,
  annotation = annotation
)

# Calculate fraction of reads in peaks
dat <- FRiP(
  dat,
  assay = "peaks", 
  total.fragments = "fragments"
)

# Fraction of reads in TSS region ----------------------------------------------

# Get the TSS sites
TSS_sites <- GetTSSPositions(
  dat@assays$ATAC@annotation
)

# # Get the 2000 bp upstream and 200 bp downstream region with promoter
# # This is how TSS was defined in Calderon et al. 2022
promoter_region <- promoters(TSS_sites)

dat$nCount_TSS <- CountsInRegion(
  dat, 
  assay = "ATAC", 
  regions = TSS_sites
)

dat$nCount_promoter <- CountsInRegion(
  dat, 
  assay = "ATAC", 
  regions = promoter_region
)

# Calculate
dat$FRiT <- dat$nCount_TSS/dat$fragments

# Save seurat object with added info
saveRDS(dat, here::here("data/processed/seurat_object/06_dat_qc.rds"))
