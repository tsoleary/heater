# ------------------------------------------------------------------------------
# Explore differential expression across time
# November 14, 2022
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(tidyverse)
require(Seurat)

# Print message
print("Starting script")

# Load data
dat <- DietSeurat(readRDS(here::here("calderon/data/main.Rds")))

# Create subset list of cells and save
set.seed(1)
subset_cells_100k <- sample(Cells(dat), 100000)
saveRDS(subset_cells_100k,here::here("calderon/data/subset_cells_100k.rds"))

# Subset
dat <- subset(dat, cell %in% subset_cells_100k)

# Subset data into hour long chunks -----

# Create list to store subsets
dats <- vector(mode = "list", length = 18)

# Loop through and create subsets
for (i in 1:length(dats)){
  # Subset and save in dats
  print(paste("subsetting hour", i))
  dats[[i]] <- DietSeurat(subset(dat, NNv1_age >= i & NNv1_age < i + 1))
}

# Print message
print("Done subsetting")

# Perform integration
dat_anchors <- FindIntegrationAnchors(object.list = dats, 
                                      dims = 1:20)
dat_subs <- IntegrateData(anchorset = dat_anchors,
                          dims = 1:20)

# Dimensional reduction and UMAP projection for subsetted data -----------------

# Run the standard workflow for visualization and clustering
dat_subs <- ScaleData(dat_subs, verbose = FALSE)
dat_subs <- RunPCA(dat_subs, npcs = 30, verbose = FALSE)

# t-SNE and Clustering
dat_subs <- RunUMAP(dat_subs, reduction = "pca", dims = 1:20)
dat_subs <- FindNeighbors(dat_subs, reduction = "pca", dims = 1:20)
dat_subs <- FindClusters(dat_subs, resolution = 0.5)

# Save dat_subs object
saveRDS(dat_subs, here::here("calderon/data/dat_subs.rds"))

# Print message
print("Done script")

# Add time windows to meta.data
# Split the data by the two time windows for plotting
# Not sure if this info exists somewhere else?
# dat_comb@meta.data <- dat_comb@meta.data %>%
#   mutate(time_window = NNv1_age)


# Visualization of UMAP projection


