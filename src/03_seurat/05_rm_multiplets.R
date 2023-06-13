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
  ) |>
  filter(q.value < 0.05) |>
  rownames()
doublet_finder <- readRDS(here::here("data/processed/qc/doublet_finder.rds"))

# Load seurat object
dat <- readRDS(
  here::here("data/processed/seurat_object/02_dat_initial_cluster.rds")
) 

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
dat <- dat |>
  subset(doublet_amulet == "Singlet" & doublet_finder == "Singlet")

# There are still two outlier barcodes with exceptionally high nCount_RNA when
# compared to the remaining high quality nuclei – 12-13k compared to 6-7k for 
# the next highest cells -- so let's remove those two outliers now
dat <- dat |> 
  subset(nCount_RNA < 10000)

# Remove cluster info from that initial clustering
dat@meta.data <- dat@meta.data |>
  select(-seurat_clusters) |>
  select(!contains("snn"))

# Save dat with multiplets removed
saveRDS(dat, here::here("data/processed/seurat_object/05_dat_filtered.rds"))
