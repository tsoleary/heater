# ------------------------------------------------------------------------------
# Remove multiplets
# May 23, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(tidyverse)
require(Seurat)

# Load multiplet barcodes
doublet_amulet <- readRDS(
  here::here("data/processed/qc/doublet_amulet.rds")
  ) %>%
  filter(q.value < 0.05) %>%
  rownames()
doublet_finder <- readRDS(here::here("data/processed/qc/doublet_finder.rds"))

# Load seurat object
dat <- readRDS(
  here::here("data/processed/seurat_object/02_dat_initial_cluster.rds")
  ) %>%
  DietSeurat()

# Mark amulet described multiplet barcodes
dat$doublet_amulet <- ifelse(colnames(dat) %in% doublet_amulet, 
                             "Multiplet",
                             "Singlet")

# Mark amulet described multiplet barcodes
dat$doublet_finder <- ifelse(colnames(dat) %in% doublet_finder,
                             "Multiplet",
                             "Singlet")

# Save dat with multiplets marked
saveRDS(dat, here::here("data/processed/seurat_object/05_dat_multiplets.rds"))

# Remove all multiplets from the data set
dat <- dat %>%
  subset(doublet_amulet == "Singlet" & doublet_finder == "Singlet")

# Remove cluster info
dat@meta.data <- dat@meta.data %>%
  select(-seurat_clusters) %>%
  select(!contains("snn"))

# Save dat with multiplets removed
saveRDS(dat, here::here("data/processed/seurat_object/05_dat_filtered.rds"))
