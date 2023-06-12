# ------------------------------------------------------------------------------
# Coverage plot of linked peaks
# June 07, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/08_dat_linked.rds"))

# Coverage plot example -----
Idents(dat) <- "seurat_clusters"
gene <- "Ldh"
p1 <- CoveragePlot(
  object = dat |> subset(acc_temp == "18°C"),
  region = gene,
  features = gene,
  expression.assay = "SCT",
  extend.upstream = 500,
  extend.downstream = 1000
)

p2 <- CoveragePlot(
  object = dat |> subset(acc_temp == "25°C"),
  region = gene,
  features = gene,
  expression.assay = "SCT",
  extend.upstream = 500,
  extend.downstream = 1000
) + cowplot::theme_cowplot(font_size = 24)

# Plot the coverage
p <- cowplot::plot_grid(p1, p2,
                        labels = c("A. 18°C", "B. 25°C"))

# Create title for the plot
title <- cowplot::ggdraw() +
  cowplot::draw_label(
    paste0("Chromatin accessibility and expression of ", gene),
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 7)
  )

# Put plot together with the title
cowplot::plot_grid(
  title, p,
  ncol = 1,
  rel_heights = c(0.1, 1)
)