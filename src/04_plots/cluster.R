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
dat <- readRDS(here::here("data/processed/seurat_object/07_dat_cluster.rds"))
dat <- readRDS(here::here("data/processed/seurat_object/08_dat_linked.rds"))

# Quick bar plot counting the number of cells for each cluster -----------------
dat@meta.data |>
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

# # Joint UMAP projection ------------------------------------------------------
# Acclimation temperature conditions plotted on top of each other 
DimPlot(dat, 
        group.by = "acc_temp") +
  scale_color_manual(name = element_blank(),
                    values = c("#43aa8b", "#f3722c")) +
  labs(title = element_blank()) +
  theme_void() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(byrow = TRUE,
                              nrow = 1,
                              override.aes = list(size = 2)))

# Save plot
ggsave(here::here("output/figs/cluster/umap_18_25_overlay.pdf"),
       height = 20,
       width = 20,
       units = "cm")
ggsave(here::here("output/figs/cluster/umap_18_25_overlay.png"),
       height = 20,
       width = 20,
       units = "cm")

# Acclimation temperature split apart 
DimPlot(dat, 
        split.by = "acc_temp",
        label = TRUE,
        repel = TRUE,
        label.color = "grey20") +
  theme_void() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 16, face = "bold"),
        strip.background = element_rect(fill = "grey95", color = "grey95")
        ) +
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


# RNA only UMAP projection -----------------------------------------------------
# Load RNA UMAP
dat <- readRDS(
  here::here("data/processed/seurat_object/07_dat_rna_umap.rds")
)

# Acclimation temperature conditions plotted on top of each other 
DimPlot(dat, 
        group.by = "acc_temp") +
  scale_color_manual(name = element_blank(),
                     values = c("#43aa8b", "#f3722c")) +
  labs(title = element_blank()) +
  theme_void() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(byrow = TRUE,
                              nrow = 1,
                              override.aes = list(size = 2)))

# Save plot
ggsave(here::here("output/figs/cluster/umap_18_25_overlay_rna.pdf"),
       height = 20,
       width = 20,
       units = "cm")
ggsave(here::here("output/figs/cluster/umap_18_25_overlay_rna.png"),
       height = 20,
       width = 20,
       units = "cm")

# Acclimation temperature split apart 
DimPlot(dat, 
        split.by = "acc_temp",
        label = TRUE,
        repel = TRUE,
        label.color = "grey20") +
  theme_void() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 16, face = "bold"),
        strip.background = element_rect(fill = "grey95", color = "grey95")
  ) +
  guides(color = guide_legend(byrow = TRUE,
                              nrow = 1,
                              override.aes = list(size = 2)))

# Save plot
ggsave(here::here("output/figs/cluster/umap_18_25_split_rna.pdf"),
       height = 15,
       width = 30,
       units = "cm")
ggsave(here::here("output/figs/cluster/umap_18_25_split_rna.png"),
       height = 15,
       width = 30,
       units = "cm")


# ATAC only UMAP projection ----------------------------------------------------
# Load ATAC UMAP
dat <- readRDS(
  here::here("data/processed/seurat_object/07_dat_atac_umap.rds")
)

# Acclimation temperature conditions plotted on top of each other 
DimPlot(dat, 
        group.by = "acc_temp") +
  scale_color_manual(name = element_blank(),
                     values = c("#43aa8b", "#f3722c")) +
  labs(title = element_blank()) +
  theme_void() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(byrow = TRUE,
                              nrow = 1,
                              override.aes = list(size = 2)))

# Save plot
ggsave(here::here("output/figs/cluster/umap_18_25_overlay_atac.pdf"),
       height = 20,
       width = 20,
       units = "cm")
ggsave(here::here("output/figs/cluster/umap_18_25_overlay_atac.png"),
       height = 20,
       width = 20,
       units = "cm")

# Acclimation temperature split apart 
DimPlot(dat, 
        split.by = "acc_temp",
        label = TRUE,
        repel = TRUE,
        label.color = "grey20") +
  theme_void() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 16, face = "bold"),
        strip.background = element_rect(fill = "grey95", color = "grey95")
  ) +
  guides(color = guide_legend(byrow = TRUE,
                              nrow = 1,
                              override.aes = list(size = 2)))

# Save plot
ggsave(here::here("output/figs/cluster/umap_18_25_split_atac.pdf"),
       height = 15,
       width = 30,
       units = "cm")
ggsave(here::here("output/figs/cluster/umap_18_25_split_atac.png"),
       height = 15,
       width = 30,
       units = "cm")
