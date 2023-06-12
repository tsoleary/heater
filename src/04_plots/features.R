# ------------------------------------------------------------------------------
# Create feature plots
# May 16, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/07_dat_cluster.rds"))

# nCount_RNA -----
FeaturePlot(dat, 
            "nCount_RNA", 
            cols = c("grey95", "darkgreen"),
            keep.scale = "feature") +
  cowplot::theme_nothing() +
  theme(title = element_blank(),
        legend.position = "right")

# Save plot
ggsave(here::here("output/figs/features/qc/nCount_RNA.png"),
       height = 15,
       width = 15,
       units = "cm")

# nCount_ATAC -----
FeaturePlot(dat, 
            "nCount_ATAC", 
            cols = c("grey95", "darkgreen"),
            keep.scale = "feature") +
  cowplot::theme_nothing() +
  theme(title = element_blank(),
        legend.position = "right")

# Save plot
ggsave(here::here("output/figs/features/qc/nCount_ATAC.png"),
       height = 15,
       width = 15,
       units = "cm")

# TSS.enrichment -----
FeaturePlot(dat, 
            "TSS.enrichment", 
            cols = c("grey95", "darkgreen"),
            keep.scale = "feature") +
  cowplot::theme_nothing() +
  theme(title = element_blank(),
        legend.position = "right")

# Save plot
ggsave(here::here("output/figs/features/qc/TSS.enrichment.png"),
       height = 15,
       width = 15,
       units = "cm")

# percent.mt -----
FeaturePlot(dat, 
            "percent.mt", 
            cols = c("grey95", "darkgreen"),
            keep.scale = "feature") +
  cowplot::theme_nothing() +
  theme(title = element_blank(),
        legend.position = "right")

# Save plot
ggsave(here::here("output/figs/features/qc/percent.mt.png"),
       height = 15,
       width = 15,
       units = "cm")

# FRiP -----
FeaturePlot(dat, 
            "FRiP", 
            cols = c("grey95", "darkgreen"),
            keep.scale = "feature") +
  cowplot::theme_nothing() +
  theme(title = element_blank(),
        legend.position = "right")

# Save plot
ggsave(here::here("output/figs/features/qc/FRiP.png"),
       height = 15,
       width = 15,
       units = "cm")

# FRiT -----
FeaturePlot(dat, 
            "FRiT", 
            cols = c("grey95", "darkgreen"),
            keep.scale = "feature") +
  cowplot::theme_nothing() +
  theme(title = element_blank(),
        legend.position = "right")

# Save plot
ggsave(here::here("output/figs/features/qc/FRiT.png"),
       height = 15,
       width = 15,
       units = "cm")

# Gene example -----------------------------------------------------------------

DefaultAssay(dat) <- "SCT"
gene <- "lmd"

# Plot expression level on top of the UMAP projection
FeaturePlot(dat, 
            gene,
            cols = c("grey95", "darkgreen"),
            reduction = "umap",
            keep.scale = "all",
            split.by = "acc_temp")

# Violin Plot
Idents(dat) <- "seurat_clusters"
VlnPlot(dat,
        gene,
        cols = c("#43aa8b", "#f3722c"),
        split.by = "acc_temp", 
        pt.size = 1) +
  scale_y_continuous(expand = c(0, 0))

