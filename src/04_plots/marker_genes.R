# ------------------------------------------------------------------------------
# Create plots for marker genes
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)

# Load plot themes 
source(here::here("src/04_plots/00_plot_themes.R"))

# # Output fig dir
# fig_dir <- "output/figs/markers"

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/07_dat_cluster.rds"))
markers <- readRDS(here::here("output/markers/cluster_markers.rds"))

# Set resolution that was used for marker gene detection
Idents(dat) <- "wknn_res.0.7"

# Top three markers for each cluster
top_markers <- markers |> 
  group_by(cluster) |> 
  slice_min(max_pval, n = 3, with_ties = FALSE)

# Output fig dir
fig_dir <- "output/figs/marker_genes"

# Gene example -----------------------------------------------------------------

# Set assay for expression data and clusters
DefaultAssay(dat) <- "SCT"

# Loop through all top genes
for (i in 1:nrow(top_markers)) {
  gene <- top_markers$gene[i]
  
  # Violin Plot
  p1 <- VlnPlot(dat,
                top_markers$gene[i],
                slot = "counts",
                cols = acc_colors,
                split.by = "acc_temp", 
                pt.size = 0) +
    labs(
      title = 
        paste0(top_markers$gene[i], 
               ": Cluster ", 
               top_markers$cluster[i],
               " marker")) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(name = "Cluster") +
    cowplot::theme_minimal_hgrid()
  
  # Expression level on top of the UMAP projection
  p2 <- FeaturePlot(dat, 
                    top_markers$gene[i],
                    cols = c("grey85", "darkgreen"),
                    reduction = "umap",
                    keep.scale = "all",
                    split.by = "acc_temp") &
    theme_void() &
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold")) 
  
  # Combine the two plots
  cowplot::plot_grid(
    rel_heights = c(1, 1.25),
    p1, p2,
    nrow = 2
  )
  
  #Save each plot
  ggsave(
    here::here(
      fig_dir,
      paste0(top_markers$cluster[i],
             "_cluster_", 
             top_markers$gene[i],
             ".png")),
         height = 30,
         width = 25,
         units = "cm")
}
