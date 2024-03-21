# ------------------------------------------------------------------------------
# Cluster and expression plots
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)
library(clustree)

# Output fig dir
fig_dir <- "output/figs/cluster"

# Load plot themes 
source(here::here("src/04_plots/00_plot_themes.R"))

# Load data --------------------------------------------------------------------
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

# Look at the different clustering resolutions 
clustree(
  dat,
  prefix = "wknn_res."
)

ggsave(here::here(fig_dir, "qc", "clustree.pdf"),
       height = 24,
       width = 24,
       units = "cm")
ggsave(here::here(fig_dir, "qc", "clustree.png"),
       height = 24,
       width = 24,
       units = "cm")

# DepthCorr for plot for which components to include downstream  ---------------
# The first lsi component often captures sequencing depth rather than biological 
# variation so we should not use this first component in further analysis.
# See plot below for ~1 corr with depth and the first component

DepthCor(dat) +
  cowplot::theme_minimal_grid()

ggsave(here::here(fig_dir, "qc", "lsi_depthcor.pdf"),
       height = 20,
       width = 25,
       units = "cm")
ggsave(here::here(fig_dir, "qc", "lsi_depthcor.png"),
       height = 20,
       width = 25,
       units = "cm")

# Quick elbow plot of the first 50 PCs -----------------------------------------
# Keeping all 50 PCs
# The elbow plot does not seem to entirely tail off before 50 PCs and the 
# SCTransform is a more robust normalization method that is less likely to 
# carry artifacts do to technical variation. Please see the below two links
# https://satijalab.org/seurat/archive/v3.0/sctransform_vignette.html
# https://www.biostars.org/p/423306/

ElbowPlot(dat, ndims = 50) +
  cowplot::theme_minimal_grid()

ggsave(here::here(fig_dir, "qc", "elbow_plot.png"),
       height = 20,
       width = 25,
       units = "cm")
ggsave(here::here(fig_dir, "qc", "elbow_plot.pdf"),
       height = 20,
       width = 25,
       units = "cm")

# DimPlot projections ----------------------------------------------------------

# Joint Dim Reduction
# Joint UMAP
DimPlot(dat,
        group.by = "acc_temp") +
  scale_color_manual(name = element_blank(),
                     values = acc_colors) +
  labs(title = element_blank()) +
  theme_void() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(byrow = TRUE,
                              nrow = 1,
                              override.aes = list(size = 3)))

# Save plot
ggsave(here::here(fig_dir, "overlay", "joint_umap.pdf"),
       height = 20,
       width = 20,
       units = "cm")
ggsave(here::here(fig_dir, "overlay", "joint_umap.png"),
       height = 20,
       width = 20,
       units = "cm")

# Joint t-SNE
DimPlot(dat,
        reduction = "tsne",
        group.by = "acc_temp",
        pt.size = 0.2) +
  scale_color_manual(name = element_blank(),
                     values = acc_colors) +
  labs(title = element_blank()) +
  theme_void() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(byrow = TRUE,
                              nrow = 1,
                              override.aes = list(size = 2)))

# Save plot
ggsave(here::here(fig_dir, "overlay", "joint_tsne.pdf"),
       height = 20,
       width = 20,
       units = "cm")
ggsave(here::here(fig_dir, "overlay", "joint_tsne.png"),
       height = 20,
       width = 20,
       units = "cm")

# RNA only -----

# RNA-only UMAP
DimPlot(dat_rna,
        group.by = "acc_temp") +
  scale_color_manual(name = element_blank(),
                     values = acc_colors) +
  labs(title = element_blank()) +
  theme_void() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(byrow = TRUE,
                              nrow = 1,
                              override.aes = list(size = 2)))

# Save plot
ggsave(here::here(fig_dir, "overlay", "rna_umap.pdf"),
       height = 20,
       width = 20,
       units = "cm")
ggsave(here::here(fig_dir, "overlay", "rna_umap.png"),
       height = 20,
       width = 20,
       units = "cm")


# RNA-only t-SNE
DimPlot(dat_rna,
        reduction = "tsne",
        group.by = "acc_temp",
        pt.size = 0.2) +
  scale_color_manual(name = element_blank(),
                     values = acc_colors) +
  labs(title = element_blank()) +
  theme_void() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(byrow = TRUE,
                              nrow = 1,
                              override.aes = list(size = 2)))

# Save plot
ggsave(here::here(fig_dir, "overlay", "rna_tsne.pdf"),
       height = 20,
       width = 20,
       units = "cm")
ggsave(here::here(fig_dir, "overlay", "rna_tsne.png"),
       height = 20,
       width = 20,
       units = "cm")

# ATAC only -----

# Acclimation temperature conditions plotted on top of each other
DimPlot(dat_atac,
        group.by = "acc_temp") +
  scale_color_manual(name = element_blank(),
                     values = acc_colors) +
  labs(title = element_blank()) +
  theme_void() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(byrow = TRUE,
                              nrow = 1,
                              override.aes = list(size = 2)))

# Save plot
ggsave(here::here(fig_dir, "overlay", "atac_umap.pdf"),
       height = 20,
       width = 20,
       units = "cm")
ggsave(here::here(fig_dir, "overlay", "atac_umap.png"),
       height = 20,
       width = 20,
       units = "cm")

# ATAC-only t-SNE
DimPlot(dat_atac,
        reduction = "tsne",
        group.by = "acc_temp",
        pt.size = 0.5) +
  scale_color_manual(name = element_blank(),
                     values = acc_colors) +
  labs(title = element_blank()) +
  theme_void() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(byrow = TRUE,
                              nrow = 1,
                              override.aes = list(size = 2)))

# Save plot
ggsave(here::here(fig_dir, "overlay", "atac_tsne.pdf"),
       height = 20,
       width = 20,
       units = "cm")
ggsave(here::here(fig_dir, "overlay", "atac_tsne.png"),
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
  fig_dir_res <- paste0(fig_dir, "/res/", "res_", res)
  
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
  ggsave(here::here(fig_dir_res, "joint_umap.pdf"),
         height = 15,
         width = 30,
         units = "cm")
  ggsave(here::here(fig_dir_res, "joint_umap.png"),
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
  ggsave(here::here(fig_dir_res, "rna_umap.pdf"),
         height = 15,
         width = 30,
         units = "cm")
  ggsave(here::here(fig_dir_res, "rna_umap.png"),
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
  ggsave(here::here(fig_dir_res, "atac_umap.pdf"),
         height = 15,
         width = 30,
         units = "cm")
  ggsave(here::here(fig_dir_res, "atac_umap.png"),
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
                      values = acc_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    cowplot::theme_minimal_hgrid()
  
  # Save plot
  ggsave(here::here(fig_dir_res, "cells_per_cluster.pdf"),
         height = 15,
         width = 30,
         units = "cm")
  ggsave(here::here(fig_dir_res, "cells_per_cluster.png"),
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
  ggsave(here::here(fig_dir_res, "rna_vln.pdf"),
         height = 10,
         width = 20,
         units = "cm")
  ggsave(here::here(fig_dir_res, "rna_vln.png"),
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
  ggsave(here::here(fig_dir_res, "atac_vln.pdf"),
         height = 10,
         width = 20,
         units = "cm")
  ggsave(here::here(fig_dir_res, "atac_vln.png"),
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
  ggsave(here::here(fig_dir_res, "frip_vln.pdf"),
         height = 10,
         width = 20,
         units = "cm")
  ggsave(here::here(fig_dir_res, "frip_vln.png"),
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
  ggsave(here::here(fig_dir_res, "frit_vln.pdf"),
         height = 10,
         width = 20,
         units = "cm")
  ggsave(here::here(fig_dir_res, "frit_vln.png"),
         height = 10,
         width = 20,
         units = "cm")
  
  # percent.mt 
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
  ggsave(here::here(fig_dir_res, "percent.mt_vln.pdf"),
         height = 10,
         width = 20,
         units = "cm")
  ggsave(here::here(fig_dir_res, "percent.mt_vln.png"),
         height = 10,
         width = 20,
         units = "cm")
  
  print(paste("done plots for", res))
}
