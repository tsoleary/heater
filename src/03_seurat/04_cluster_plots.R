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

# Acclimation temperature conditions plotted on top of each other --------------
DimPlot(dat, 
        group.by = "acc_temp") +
  labs(title = element_blank()) +
  scale_color_manual(name = element_blank(),
                     values = c("#43aa8b", "#f3722c")) +
  theme_void() +
  theme(legend.position = "bottom")

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

# Without labels on the plot
DimPlot(dat, 
        split.by = "acc_temp") +
  theme_void() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(byrow = TRUE,
                              nrow = 1,
                              override.aes = list(size = 2)))


