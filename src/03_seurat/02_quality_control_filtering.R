# ------------------------------------------------------------------------------
# Quality control and filtering of low-quality cells
# March 27, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Description -----
# Basic quality control metrics and filtering out non-cell barcodes.

# Load libraries
library(Seurat)
library(Signac)
library(tidyverse)

# Load  data
dat <- readRDS(here::here("data/processed/seurat_object/01_dat_raw.rds"))


# Per sample raw sequencing read data information
dat@meta.data %>%
  group_by(sample_name) %>%
  summarize(total_RNA_Reads = sum(nCount_RNA),
            total_ATAC_Reads = sum(nCount_ATAC))

# Knee-plot for ATAC -----------------------------------------------------------
# Create ranked data.frame
knee_atac_df <- tibble(
  nCount_ATAC = sort(dat$nCount_ATAC, 
                     decreasing = TRUE)
  ) %>%
  rownames_to_column("rank") %>%
  mutate(rank = as.numeric(rank)) %>%
  head(40000)

# Plot the knee_plot
knee_atac_df %>%
  ggplot() +
  geom_point(aes(x = rank,
                 y = nCount_ATAC),
             color = "grey50",
             shape = 21) +
  geom_hline(aes(yintercept = 10000),
             linetype = 2,
             color = "orange",
             alpha = 0.8,
             linewidth = 1.1) +
  geom_hline(aes(yintercept = 500),
             linetype = 2,
             color = "orange",
             alpha = 0.8,
             linewidth = 1.1) +
  scale_y_continuous(trans = "log10") +
  cowplot::theme_minimal_grid()

# Knee-plot for RNA ------------------------------------------------------------
# Create ranked data.frame
knee_rna_df <- tibble(
  nCount_RNA = sort(dat$nCount_RNA, 
                    decreasing = TRUE)
  ) %>%
  rownames_to_column("rank") %>%
  mutate(rank = as.numeric(rank)) %>%
  head(40000)

# Plot the knee_plot
knee_rna_df %>%
  ggplot() +
  geom_point(aes(x = rank,
                 y = nCount_RNA),
             color = "grey50",
             shape = 21) +
  geom_hline(aes(yintercept = 10000),
             linetype = 2,
             color = "orange",
             alpha = 0.8,
             linewidth = 1.1) +
  geom_hline(aes(yintercept = 500),
             linetype = 2,
             color = "orange",
             alpha = 0.8,
             linewidth = 1.1) +
  scale_y_continuous(trans = "log10") +
  cowplot::theme_minimal_grid()

# Save knee plot for the 
ggsave(here::here("output/figs/rna_knee_plot.pdf"))

# Filter out low quality cells and non-cell barcodes ---------------------------

# Filter criteria just based on ATAC?
dat <- subset(
  x = dat,
  subset = nCount_ATAC < 10000 &
    nCount_ATAC > 500 &
    nCount_RNA > 100 &
    nCount_RNA < 5000
)

# Total number of nuclei
dat@meta.data %>%
  tally()

# Number of nuclei at each acclimation temperature
dat@meta.data %>%
  group_by(acc_temp) %>%
  tally()

# Number of nuclei in each sample
dat@meta.data %>%
  group_by(sample_name) %>%
  tally()

# Does calderon have these knee plots -- if not I should just make them if it is possible
# Albright has good knee plots but that is pre cellularization so it is possible 
# is why they were able to reduce ambient RNA in their sample that 


# DecontX to account for ambient RNAs ------------------------------------------ 

# Convert RNA counts to a single cell experiment object 
counts <- GetAssayData(object = dat, slot = "counts")
sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))

# Run decontX on the RNA counts
sce <- celda::decontX(sce)

# Create a new assay object with dexcontX fixed counts
# Optionally write over the old counts ---??
dat[["decontX_counts"]] <- CreateAssayObject(
  counts = celda::decontXcounts(sce)
)

# Save data
saveRDS(dat, "data/processed/seurat_object/02_dat_cells.rds")

### SCRATCH
# Calculate Nucleosome Signal and TSS Enrichment for ATAC QC metrics
# DefaultAssay(dat) <- "ATAC"
# dat <- NucleosomeSignal(dat)
# dat <- TSSEnrichment(dat) # Doesn't work?

# Create simple Violin Plot with cell meta.data information before filtering
# VlnPlot(
#   object = dat,
#   features = c("nCount_RNA", 
#                "nCount_ATAC", 
#                "TSS.enrichment", 
#                "nucleosome_signal"),
#   ncol = 4,
#   pt.size = 0
# )

# VlnPlot(
#   object = dat,
#   features = c("nCount_RNA",
#                "nCount_ATAC"),
#   ncol = 4,
#   pt.size = 0
# )