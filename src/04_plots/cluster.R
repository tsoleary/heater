# ------------------------------------------------------------------------------
# Cluster and expression plots
# May 16, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)
<<<<<<< HEAD
library(clustree)

# Load data
# Joint dim reduction
dat <- readRDS(
  here::here("data/processed/seurat_object/07_dat_cluster.rds")
)
# ATAC-only dim reduction 
dat_atac <- readRDS(
  here::here("data/processed/seurat_object/07_dat_atac_only_dim.rds")
)
# RNA-only dim reduction
dat_rna <- readRDS(
  here::here("data/processed/seurat_object/07_dat_rna_only_dim.rds")
)

# Look at the different clustering
clustree(
  dat,
  prefix = "wknn_res."
)

ggsave(here::here("output/figs/cluster/clustree.pdf"),
       height = 24,
       width = 24,
       units = "cm")
ggsave(here::here("output/figs/cluster/clustree.png"),
       height = 24,
       width = 24,
       units = "cm")

# DimPlot projections ----------------------------------------------------------

# Joint Dim Reduction
# Joint UMAP
DimPlot(dat,
        group.by = "acc_temp") +
  scale_color_manual(name = element_blank(),
                     values = c("#43aa8b", "#f3722c")) +
=======

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
>>>>>>> a8de26aef8a500cdbd82fa277045d3e778c75e21
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

<<<<<<< HEAD
# Joint t-SNE
DimPlot(dat,
        reduction = "tsne",
        group.by = "acc_temp",
        pt.size = 0.5) +
  scale_color_manual(name = element_blank(),
                     values = c("#43aa8b", "#f3722c")) +
  labs(title = element_blank()) +
  theme_void() +
  theme(legend.position = "bottom") +
=======
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
>>>>>>> a8de26aef8a500cdbd82fa277045d3e778c75e21
  guides(color = guide_legend(byrow = TRUE,
                              nrow = 1,
                              override.aes = list(size = 2)))

# Save plot
<<<<<<< HEAD
ggsave(here::here("output/figs/cluster/tsne_18_25_overlay.pdf"),
       height = 20,
       width = 20,
       units = "cm")
ggsave(here::here("output/figs/cluster/tsne_18_25_overlay.png"),
       height = 20,
       width = 20,
       units = "cm")

# RNA only -----

# RNA-only UMAP
DimPlot(dat_rna,
=======
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
>>>>>>> a8de26aef8a500cdbd82fa277045d3e778c75e21
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

<<<<<<< HEAD

# RNA-only t-SNE
DimPlot(dat_rna,
        reduction = "tsne",
        group.by = "acc_temp",
        pt.size = 0.5) +
  scale_color_manual(name = element_blank(),
                     values = c("#43aa8b", "#f3722c")) +
  labs(title = element_blank()) +
  theme_void() +
  theme(legend.position = "bottom") +
=======
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
>>>>>>> a8de26aef8a500cdbd82fa277045d3e778c75e21
  guides(color = guide_legend(byrow = TRUE,
                              nrow = 1,
                              override.aes = list(size = 2)))

# Save plot
<<<<<<< HEAD
ggsave(here::here("output/figs/cluster/tsne_18_25_overlay_rna.pdf"),
       height = 20,
       width = 20,
       units = "cm")
ggsave(here::here("output/figs/cluster/tsne_18_25_overlay_rna.png"),
       height = 20,
       width = 20,
       units = "cm")

# ATAC only -----

# Acclimation temperature conditions plotted on top of each other
DimPlot(dat_atac,
=======
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
>>>>>>> a8de26aef8a500cdbd82fa277045d3e778c75e21
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

<<<<<<< HEAD
# ATAC-only t-SNE
DimPlot(dat_atac,
        reduction = "tsne",
        group.by = "acc_temp",
        pt.size = 0.5) +
  scale_color_manual(name = element_blank(),
                     values = c("#43aa8b", "#f3722c")) +
  labs(title = element_blank()) +
  theme_void() +
  theme(legend.position = "bottom") +
=======
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
>>>>>>> a8de26aef8a500cdbd82fa277045d3e778c75e21
  guides(color = guide_legend(byrow = TRUE,
                              nrow = 1,
                              override.aes = list(size = 2)))

# Save plot
<<<<<<< HEAD
ggsave(here::here("output/figs/cluster/tsne_18_25_overlay_atac.pdf"),
       height = 20,
       width = 20,
       units = "cm")
ggsave(here::here("output/figs/cluster/tsne_18_25_overlay_atac.png"),
       height = 20,
       width = 20,
       units = "cm")


# Cluster resolution specific plotting -----------------------------------------

# Set resolutions
res_list <- as.character(seq(from = 0.1, to = 1, by = 0.1))

# Loop through and plot all resolutions
for (res in res_list) {
  
  print(paste("starting plots for", res))
  
  # Set seurat clusters to the specific resolution
  dat@meta.data$seurat_clusters <- dat@meta.data[, paste0("wknn_res.", res)]
  dat_rna@meta.data$seurat_clusters <- 
    dat_rna@meta.data[, paste0("wknn_res.", res)]
  dat_atac@meta.data$seurat_clusters <- 
    dat_atac@meta.data[, paste0("wknn_res.", res)]
  Idents(dat) <- "seurat_clusters"
  Idents(dat_rna) <- "seurat_clusters"
  Idents(dat_atac) <- "seurat_clusters"
  
  # Set output fig_dir
  fig_dir <- paste0("output/figs/cluster/res_", res)
  
  # UMAP acclimation temperature split apart ---------
  # Joint
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
  ggsave(here::here(fig_dir, "umap_18_25_split.pdf"),
         height = 15,
         width = 30,
         units = "cm")
  ggsave(here::here(fig_dir, "umap_18_25_split.png"),
         height = 15,
         width = 30,
         units = "cm")
  
  # RNA-only UMAP
  DimPlot(dat_rna, 
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
  ggsave(here::here(fig_dir, "rna_umap_18_25_split.pdf"),
         height = 15,
         width = 30,
         units = "cm")
  ggsave(here::here(fig_dir, "rna_umap_18_25_split.png"),
         height = 15,
         width = 30,
         units = "cm")
  
  # RNA-only UMAP
  DimPlot(dat_atac, 
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
  ggsave(here::here(fig_dir, "atac_umap_18_25_split.pdf"),
         height = 15,
         width = 30,
         units = "cm")
  ggsave(here::here(fig_dir, "atac_umap_18_25_split.png"),
         height = 15,
         width = 30,
         units = "cm")
  
  # Quick bar plot counting the number of cells for each cluster --------
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
  ggsave(here::here(fig_dir, "cells_per_cluster.pdf"),
         height = 15,
         width = 30,
         units = "cm")
  ggsave(here::here(fig_dir, "cells_per_cluster.png"),
         height = 15,
         width = 30,
         units = "cm")
  
  # Cluster QC -----
  # RNA Count 
  dat@meta.data |>
    group_by(seurat_clusters) |> 
    ggplot(aes(y = nCount_RNA,
               x = seurat_clusters)) +
    geom_violin(color = "grey20",
                fill = "grey80",
                width = 0.8)+
    scale_y_continuous(expand = c(0, 0),
                       name = "RNA Counts") +
    cowplot::theme_cowplot() 
  
  # Save plot
  ggsave(here::here(fig_dir, "rna_clust_violin.pdf"),
         height = 10,
         width = 20,
         units = "cm")
  ggsave(here::here(fig_dir, "rna_clust_violin.png"),
         height = 10,
         width = 20,
         units = "cm")
  
  # ATAC Count 
  dat@meta.data |>
    group_by(seurat_clusters) |> 
    ggplot(aes(y = nCount_ATAC,
               x = seurat_clusters)) +
    geom_violin(color = "grey20",
                fill = "grey80",
                width = 0.8)+
    scale_y_continuous(expand = c(0, 0),
                       name = "ATAC Counts") +
    cowplot::theme_cowplot() 
  
  # Save plot
  ggsave(here::here(fig_dir, "atac_clust_violin.pdf"),
         height = 10,
         width = 20,
         units = "cm")
  ggsave(here::here(fig_dir, "atac_clust_violin.png"),
         height = 10,
         width = 20,
         units = "cm")
  
  
  # FRiP Count 
  dat@meta.data |>
    group_by(seurat_clusters) |> 
    ggplot(aes(y = FRiP,
               x = seurat_clusters)) +
    geom_violin(color = "grey20",
                fill = "grey80",
                width = 0.8)+
    scale_y_continuous(expand = c(0, 0),
                       name = "FRiP") +
    cowplot::theme_cowplot() 
  
  # Save plot
  ggsave(here::here(fig_dir, "frip_clust_violin.pdf"),
         height = 10,
         width = 20,
         units = "cm")
  ggsave(here::here(fig_dir, "frip_clust_violin.png"),
         height = 10,
         width = 20,
         units = "cm")
  
  # FRiT Count 
  dat@meta.data |>
    group_by(seurat_clusters) |> 
    ggplot(aes(y = FRiT,
               x = seurat_clusters)) +
    geom_violin(color = "grey20",
                fill = "grey80",
                width = 0.8)+
    scale_y_continuous(expand = c(0, 0),
                       name = "FRiT") +
    cowplot::theme_cowplot() 
  
  # Save plot
  ggsave(here::here(fig_dir, "frit_clust_violin.pdf"),
         height = 10,
         width = 20,
         units = "cm")
  ggsave(here::here(fig_dir, "frit_clust_violin.png"),
         height = 10,
         width = 20,
         units = "cm")
  
  # percent.mt Count 
  dat@meta.data |>
    group_by(seurat_clusters) |> 
    ggplot(aes(y = percent.mt,
               x = seurat_clusters)) +
    geom_violin(color = "grey20",
                fill = "grey80",
                width = 0.8)+
    scale_y_continuous(expand = c(0, 0),
                       name = "percent.mt") +
    cowplot::theme_cowplot() 
  
  # Save plot
  ggsave(here::here(fig_dir, "percent.mt_clust_violin.pdf"),
         height = 10,
         width = 20,
         units = "cm")
  ggsave(here::here(fig_dir, "percent.mt_clust_violin.png"),
         height = 10,
         width = 20,
         units = "cm")
  
  print(paste("done plots for", res))
}






=======
ggsave(here::here("output/figs/cluster/umap_18_25_split_atac.pdf"),
       height = 15,
       width = 30,
       units = "cm")
ggsave(here::here("output/figs/cluster/umap_18_25_split_atac.png"),
       height = 15,
       width = 30,
       units = "cm")
>>>>>>> a8de26aef8a500cdbd82fa277045d3e778c75e21
