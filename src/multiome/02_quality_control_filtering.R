# ------------------------------------------------------------------------------
# Quality control and filtering of low-quality cells
# March 27, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Description -----
# Create Quality Control Figures for supplement and filter out 

# Filter criteria
### Figure out flexible way to do this
min_ATAC_count <- 100000
min_RNA_count <- 25000
max_nucleosome_signal <- 2
min_TSS.enrichment <- 1


# Load libraries
library(Seurat)

# Load  data
dat <- readRDS("data/processed/seurat_object_01_initial.rds")

# Analyze data
DefaultAssay(pbmc) <- "ATAC"

pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

# Save QC Violin Plot figure
ggsave2()

# Filter out low quality cells
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
pbmc

# Save data
saveRDS(dat, "data/processed/seurat_object_02_filtered.rds")