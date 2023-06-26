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
dat <- readRDS(here::here("data/processed/seurat_object/10_dat_linked.rds"))

# Coverage plot example -----
Idents(dat) <- "cell_type"
gene <- "CG13427"

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

ggsave(here::here("output/figs/degs", gene, "coverage.png"),
       height = 25,
       width = 40,
       units = "cm")

# dat |> 
#   subset(cell_type == "mesoderm prim.") |> 
#   CoveragePlot(
#     region = gene,
#     features = gene,
#     split.by = "acc_temp",
#     expression.assay = "SCT",
#     extend.upstream = 500,
#     extend.downstream = 1000
#   )



# Violin Plot
VlnPlot(dat,
        features = gene,
        split.by = "acc_temp",
        pt.size = 0,
        assay = "RNA") +
  scale_y_continuous(expand = c(0, 0.05),
                     position = "right") +
  scale_x_discrete(name = element_blank()) +
  scale_fill_manual(values = c("#43aa8b", "#f3722c")) +
  coord_flip() +
  cowplot::theme_cowplot()


ggsave(here::here("output/figs/degs", gene, "vln.png"),
       height = 25,
       width = 20,
       units = "cm")

