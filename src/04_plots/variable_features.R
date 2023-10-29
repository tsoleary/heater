# ------------------------------------------------------------------------------
# Plots: Variable features: heatmap, histogram
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)
library(cowplot)
source(here::here("src/04_plots/00_plot_themes.R"))

# Load data
dat <- readRDS(
  here::here("data/processed/seurat_object/10_dat_linked.rds")
)

# Output fig dir
fig_dir <- "output/figs"

# Make default assay the normalized gene expression
DefaultAssay(dat) <- "SCT"

# Make heatmap of RNA expression
DoHeatmap(
  dat,
  group.by = "seurat_clusters",
  features = VariableFeatures(dat)[1:300],
  cells = sample(Cells(dat), 1000)
) &
  colorspace::scale_fill_continuous_diverging()

# # Create a mini data set
# mini <- dat |> 
#   subset(cells = sample(Cells(dat), 10))
# 
# # Make a heatmap of that data
# n_cells <- 1240
# n_genes <- 300
# 
# top_genes <- HVFInfo(dat) |> 
#   mutate(rank = dense_rank(desc(residual_variance))) |> 
#   slice_min(rank, n = n_genes) |> 
#   rownames()
# 
# dat@assays$SCT@scale.data[top_genes, ] |> 
#   pheatmap::pheatmap(
#     scale = "row",
#     color = colorspace::diverge_hcl(palette = "Blue-Red", n = 100),
#     show_rownames = FALSE,
#     show_colnames = FALSE, 
#     border_color = NA)

# Save
ggsave(here::here(fig_dir, "heatmap.png"),
       height = 20,
       width = 20,
       units = "cm")
ggsave(here::here(fig_dir, "heatmap.pdf"),
       height = 20,
       width = 20,
       units = "cm")

DefaultAssay(dat) <- "SCT"
HVFInfo(dat) |> arrange(desc(gmean))

HVFInfo(dat) |> 
  ggplot(aes(x=gmean, y = residual_variance)) + 
  geom_point() +
  scale_y_log10() +
  scale_x_log10()

# Residual variance in rank order
HVFInfo(dat) |> 
  mutate(rank = dense_rank(desc(residual_variance))) |> 
  ggplot() +
  geom_point(aes(x = rank,
                 y= residual_variance)) +
  geom_vline(xintercept = 200,
             color = "orange",
             linetype = 2,
             linewidth = 1.1) +
  scale_x_log10(name = "Genes in rank order of variability") +
  scale_y_continuous(name = "Residual variance") +
  theme_minimal_grid()

# Residual variance on a log scale in rank order
HVFInfo(dat) |> 
  mutate(rank = dense_rank(desc(residual_variance))) |> 
  ggplot() +
  geom_point(aes(x = rank,
                 y= residual_variance)) +
  geom_vline(xintercept = 2000,
             color = "orange",
             linetype = 2,
             linewidth = 1.1) +
  scale_x_log10(name = "Genes in rank order of variability") +
  scale_y_log10(name = "Residual variance") +
  theme_minimal_grid()


# Variable Feature Plot
VariableFeaturePlot(dat)



