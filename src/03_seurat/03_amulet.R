# ------------------------------------------------------------------------------
# AMULET to detect doublets in scATAC-Seq > 2 overlapping reads
# May 22, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(scDblFinder)

# needs higher read depth maybe?? (10-15K) max performance at 25k/cell
# median(dat@meta.data$nCount_ATAC) is ~6,032

# Path to filtered fragment file
data_dir <- here::here("data/processed/seq/all/outs/")
frag_path <- paste0(data_dir, "atac_fragments_cells.tsv.gz")

# Regions to exclude -- sex chromosomes, mito genome, non chr fragments
# Add highly repetitive regions?
exclude_regions <- GenomicRanges::GRanges(c("X", "Y", "mitochondrion_genome",
                                            "Unmapped_Scaffold_8_D1580_D1567",
                                            "211000022280328",
                                            "211000022278760"),
                                          IRanges(1L, width = 10^9))

# Run amulet on fragments excluding regions
doublet_amulet <- amulet(frag_path, regionsToExclude = exclude_regions)

# Save the list of doublets from amulet
saveRDS(doublet_amulet, here::here("data/processed/qc/doublet_amulet.rds"))
