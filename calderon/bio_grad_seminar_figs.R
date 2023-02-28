# ------------------------------------------------------------------------------
# HEATER - Calderon figs for BioLunch
# February 28, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(tidyverse)
require(Seurat)

# Load data
dat <- readRDS(here::here("calderon/data/main.Rds"))

# Subset the data to between 5 and 6 hours according to the NNv1 model
dat_5 <- subset(dat, NNv1_age > 5 & NNv1_age < 6)

# UMAP with clusters labeled with number
DimPlot(dat_5, label = TRUE)
