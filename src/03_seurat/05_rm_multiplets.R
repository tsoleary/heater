# ------------------------------------------------------------------------------
# Remove multiplets
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)

# Load multiplet barcodes
doublet_amulet <- readRDS(
  here::here("data/processed/qc/doublet_amulet.rds")
  ) |>
  dplyr::filter(q.value < 0.05) |>
  base::rownames()
doublet_finder <- readRDS(here::here("data/processed/qc/doublet_finder.rds"))

# Load seurat object
dat <- readRDS(
  here::here("data/processed/seurat_object/02_dat_initial_cluster.rds")
) 

# Output dir
out_dir <- "data/processed/seurat_object"

# Mark amulet described multiplet barcodes
dat$doublet_amulet <- ifelse(colnames(dat) %in% doublet_amulet, 
                             "Multiplet",
                             "Singlet")

# Mark amulet described multiplet barcodes
dat$doublet_finder <- ifelse(colnames(dat) %in% doublet_finder,
                             "Multiplet",
                             "Singlet")

dat@meta.data <- dat@meta.data |> 
  mutate(multiplet_call = case_when(doublet_amulet == "Singlet" & 
                                    doublet_finder == "Singlet" ~ "Singlet",
                                  .default = "Multiplet")) 



# Remove putative multiplets from the data set -----

# Make clusters at a more fine resolution
dat <- FindClusters(dat, resolution = 1.0)

# Plot of subclusters
DimPlot(dat, label = TRUE) + 
  theme_void() +
  theme(title = element_blank(),
        legend.position = "none")

# Plot of individual multiplets marked
DimPlot(dat, group = "multiplet_call") +
  scale_color_manual(name = element_blank(),
                     values = c("firebrick", "grey80")) +
  theme_void() +
  theme(title = element_blank(),
        legend.position = "bottom")

# Identify proportion of multiplets within subclusters
multiplet_prop <- dat@meta.data |> 
  group_by(seurat_clusters, multiplet_call) |> 
  tally() |> 
  pivot_wider(names_from = multiplet_call,
              values_from = n) |> 
  mutate(percent_multiplet = Multiplet/(Singlet + Multiplet)) |> 
  arrange(desc(percent_multiplet))

# Identify all clusters with greater than 15% multiplets
multiplet_clusters <- multiplet_prop |> 
  filter(percent_multiplet > 0.15) |> 
  pull(seurat_clusters)

# Identify all cells marked individually as multiplets or cells that are in 
#   clusters with at least 15% multiplets
dat@meta.data <- dat@meta.data |> 
  mutate(multiplet_rm = 
           case_when(multiplet_call == "Multiplet" |
                       seurat_clusters %in% multiplet_clusters ~ "Multiplet",
                     .default = "Singlet"))

# Save dat with multiplets marked
saveRDS(dat, here::here(out_dir, "05_dat_multiplets.rds"))

# Plot all multiplets that will be removed
DimPlot(dat, group = "multiplet_rm") +
  scale_color_manual(name = element_blank(),
                     values = c("firebrick", "grey80")) +
  theme_void() +
  theme(title = element_blank(),
        legend.position = "bottom")

# Remove all those multiplets
dat <- dat |>
  base::subset(multiplet_rm != "Multiplet")

# Remove cluster info from that initial clustering
dat@meta.data <- dat@meta.data |>
  dplyr::select(-seurat_clusters) |>
  dplyr::select(!contains("snn")) |>
  dplyr::select(!contains("doublet")) |> 
  dplyr::select(!contains("multiplet"))

# Save dat with multiplets removed
saveRDS(dat, here::here(out_dir, "05_dat_filtered.rds"))
