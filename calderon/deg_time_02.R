# ------------------------------------------------------------------------------
# Data subsets -- Visualization and differential expression
# November 18, 2022
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(tidyverse)
require(Seurat)

# Load data
dat_subs <- readRDS(here::here("calderon/data/dat_subs.rds"))

# Add time windows to meta.data
dat_subs$time_window <- as.character(floor(dat_subs$NNv1_age))

# Visualization of UMAP projection
dat_subs@meta.data %>%
  group_by(time_window) %>%
  tally() 

# Differential expression between adjacent time points -------------------------

# Set identity to time window
Idents(dat_subs) <- dat_subs$time_window

degs_all <- tibble()

for (i in 1:17) {
  
  # Print loop
  print(i)
  
  # Identify DEGs
  degs <- FindMarkers(dat_subs, 
                      ident.1 = as.character(i), 
                      ident.2 = as.character(i + 1),
                      logfc.threshold = 0.1,
                      min.pct = 0.02)

  
  if (nrow(degs) > 0) {
    # Change rownames to column
    degs <- degs %>%
      rownames_to_column("gene")
    
    # Save comparison
    degs$comp <- paste(i, "|", i+1)
    
    # Save comparisons
    degs_all <- bind_rows(degs_all, degs)
  }

}

# DEGs

degs_all %>%
  filter(p_val_adj < 0.05) %>%
  group_by(comp) %>%
  tally()

# Plot UMAP over time ----------------------------------------------------------

dat_umap <- dat_subs@reductions$umap@cell.embeddings %>%
  as_tibble() %>%
  rownames_to_column("cell")

dat_umap$time_window <- dat_subs$time_window  

dat_umap %>%
  mutate(time_window = as.numeric(time_window)) %>%
  ggplot() +
  geom_point(aes(x = UMAP_1,
                 y = UMAP_2, 
                 color = time_window)) +
  gganimate::transition_states(time_window,
                               transition_length = 2,
                               state_length = 1) +
  labs(title = "Hour: {frame}")







