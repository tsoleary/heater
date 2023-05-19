# ------------------------------------------------------------------------------
# Cluster and expression plots
# May 16, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/03_dat_clustered.rds"))

# Quick bar plot counting the number of cells for each cluster -----------------
dat@meta.data %>%
  ggplot(aes(x = seurat_clusters,
             fill = acc_temp)) +
  geom_bar(position = "dodge",
           color = "grey50") +
  labs(y = "Number of cells",
       x = "Cluster") +
  scale_fill_manual(name = element_blank(),
                    values = c("#43aa8b", "#f3722c")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  cowplot::theme_minimal_hgrid()

# Save plot
ggsave(here::here("output/figs/cluster/cells_per_cluster.pdf"),
       height = 15,
       width = 30,
       units = "cm")
ggsave(here::here("output/figs/cluster/cells_per_cluster.png"),
       height = 15,
       width = 30,
       units = "cm")

# Acclimation temperature conditions plotted on top of each other --------------
DimPlot(dat, 
        group.by = "acc_temp") +
  labs(title = element_blank()) +
  scale_color_manual(name = element_blank(),
                     values = c("#43aa8b", "#f3722c")) +
  theme_void() +
  theme(legend.position = "bottom")

# Save plot
ggsave(here::here("output/figs/cluster/umap_18_25_overlay.pdf"),
       height = 20,
       width = 20,
       units = "cm")
ggsave(here::here("output/figs/cluster/umap_18_25_overlay.png"),
       height = 20,
       width = 20,
       units = "cm")

# Acclimation temperature split apart ------------------------------------------
DimPlot(dat, 
        split.by = "acc_temp",
        label = TRUE,
        label.color = "grey20") +
  theme_void() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(byrow = TRUE,
                              nrow = 1,
                              override.aes = list(size = 2)))

# Save plot
ggsave(here::here("output/figs/cluster/umap_18_25_split_labeled.pdf"),
       height = 15,
       width = 30,
       units = "cm")
ggsave(here::here("output/figs/cluster/umap_18_25_split_labeled.png"),
       height = 15,
       width = 30,
       units = "cm")

# Without labels on the plot
DimPlot(dat, 
        split.by = "acc_temp") +
  theme_void() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(byrow = TRUE,
                              nrow = 1,
                              override.aes = list(size = 2)))

# Save plot
ggsave(here::here("output/figs/cluster/umap_18_25_split.pdf"),
       height = 15,
       width = 30,
       units = "cm")
ggsave(here::here("output/figs/cluster/umap_18_25_split.png"),
       height = 15,
       width = 30,
       units = "cm")
