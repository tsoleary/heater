# ------------------------------------------------------------------------------
# AMULET to detect doublets in scATAC-Seq > 2 overlapping reads
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(scDblFinder)

# Path to filtered fragment file
data_dir <- here::here("data/processed/seq/all/outs/")
frag_path <- paste0(data_dir, "atac_fragments_cells.tsv.gz")

# Output dir
out_dir <- "data/processed/qc"

# Regions to exclude:
# Sex chromosomes, mitochonrial genome, non-chromosomal fragments
exclude_regions <- GenomicRanges::GRanges(
  c("X", "Y", "mitochondrion_genome", "Unmapped_Scaffold_8_D1580_D1567",
    "211000022280328", "211000022278760"),
  IRanges::IRanges(1L, width = 10^9)
)

# Run amulet on fragments excluding regions
doublet_amulet <- amulet(frag_path, regionsToExclude = exclude_regions)

# Save the list of doublets from amulet
saveRDS(doublet_amulet, here::here(out_dir, "doublet_amulet.rds"))
