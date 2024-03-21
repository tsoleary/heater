# ------------------------------------------------------------------------------
# Plot WGCNA results
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)
library(cowplot)
library(hdWGCNA)

# Load plot themes 
source(here::here("src/04_plots/00_plot_themes.R"))

# Output fig directory
fig_dir <- here::here("output/figs/wgcna")

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/13_dat_wgcna.rds"))
DMEs <- readRDS(here::here("output/wgcna/DMEs.rds"))
red <- read_csv(here::here("output/wgcna/red.csv")) |> pull()
brown <- read_csv(here::here("output/wgcna/brown.csv")) |> pull()

# Plot the soft power analysis results
PlotSoftPowers(dat) |>
  wrap_plots(ncol = 2)

# Save
ggsave(here::here(fig_dir, "softpower_threshold.png"),
       height = 20,
       width = 20,
       units = "cm")

# Plot the dendrogram
png(here::here(fig_dir, "dendrogram.png"),
    res = 300,
    height = 15,
    width = 20,
    units = "cm")
PlotDendrogram(
  dat,
  main = "hdWGCNA Dendrogram"
)
dev.off()

# Lollipop DME plot ----
mod_tally <- GetModules(dat) |> 
  group_by(module) |> 
  tally()

x <- DMEs |>
  rownames_to_column("mod") |> 
  mutate(mod = fct_relevel(mod, DMEs |> arrange(avg_log2FC) |> pull(module)))

mod_colors <- c(
  "red", 
  "brown", 
  "blue", 
  "turquoise", 
  "black", 
  "yellow", 
  "green"
) |> colorspace::desaturate(amount = 0.5)

left_join(x, mod_tally) |> 
  ggplot(aes(y = mod,
             x = avg_log2FC)) +
  geom_segment(aes(xend = 0, 
                   yend = mod,
                   color = mod),
               linewidth = 1) +
  annotate(geom = "text", 
           y = "red", 
           x = 0.18, 
           color = "grey20",
           label = "*", 
           vjust = 0.8, 
           size = 10) +
  annotate(geom = "text", 
           y = "brown", 
           x = 0.18, 
           color = "grey20",
           label = "*", 
           vjust = 0.8, 
           size = 10) +
  geom_point(aes(fill = mod, 
                 size = n), 
             color = "grey80",
             shape = 21) +
  geom_vline(xintercept = 0, color = "grey50") +
  scale_x_continuous(name = "log2(fold-change)",
                     limits = c(-0.25, 0.25)) +
  scale_fill_manual(values = rev(mod_colors)) +
  scale_color_manual(values = rev(mod_colors)) +
  scale_size_continuous(range = c(5, 10)) +
  labs(y = element_blank()) +
  cowplot::theme_minimal_vgrid() +
  theme(legend.position = "none",
        panel.grid.major.y = element_line(color = "grey95"),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

# Save
ggsave(here::here(fig_dir, "DME_lolipops.png"),
       height = 12,
       width = 15,
       units = "cm")


# Feature plots for brown and red modules

# Make module umap 
MEs_umap <- GetMEs(dat) |> 
  rownames_to_column("cell") |> 
  full_join(dat@reductions[["umap"]]@cell.embeddings |> 
              as_tibble(rownames = "cell")) |> 
  full_join(dat@meta.data|> 
              as_tibble(rownames = "cell")) |> 
  select(cell, UMAP_1, UMAP_2, brown, red, acc_temp) |> 
  tibble()


# Brown module feature plot
MEs_umap_rank <- MEs_umap |> 
  mutate(rank = dense_rank(desc(brown))) 


ggplot() +
  geom_point(
    data = MEs_umap_rank |> filter(rank > 100),
    aes(x = UMAP_1,
        y = UMAP_2,
        color = brown),
        size = 0.5) +
  geom_point(
    data = MEs_umap_rank |> filter(rank < 100),
    aes(x = UMAP_1,
        y = UMAP_2,
        color = brown),
    size = 0.5) +
  scale_color_gradient2(
    low = "grey95", 
    mid = "grey20", 
    high = "grey10",
    midpoint = 10) +
  facet_wrap(~acc_temp) +
  theme_void() +
  theme(strip.background = element_rect(fill = "grey90"))

# Save
ggsave(here::here(fig_dir, "brown_umap.png"),
       height = 10,
       width = 20,
       units = "cm")

# Red module feature map
MEs_umap |> 
  tibble() |> 
  ggplot() +
  geom_point(aes(x = UMAP_1, 
                 y = UMAP_2, 
                 color = red),
             size = 0.5) +
  scale_color_gradient2(
    low = "grey95", 
    mid = "grey80",
    high = "grey20",
    midpoint = 0) +
  facet_wrap(~acc_temp) +
  theme_void() +
  theme(strip.background = element_rect(fill = "grey90"))

# Save
ggsave(here::here(fig_dir, "red_umap.png"),
       height = 10,
       width = 20,
       units = "cm")

# # get  reduction from seurat obj
# umap <- dat@reductions[["umap"]]@cell.embeddings
# x_name <- colnames(umap)[1]
# y_name <- colnames(umap)[2]
# FeaturePlot(
#   dat,
#   GetMEs(dat)
# )

# Make featureplot of MEs for each module
ModuleFeaturePlot(
  dat,
  features = "MEs",
  order = TRUE) |>
  wrap_plots(ncol = 3)

# Save
ggsave(here::here(fig_dir, "module_feature_plot.png"),
       height = 20,
       width = 20,
       units = "cm")

GetModules(dat)

# Module correlogram
png(here::here(fig_dir, "module_correlogram.png"),
    res = 300,
    height = 20,
    width = 20,
    units = "cm")
ModuleCorrelogram(dat)
dev.off()

MEs <- GetMEs(dat) |> 
  rownames_to_column("cell") |> 
  full_join(dat@meta.data |> 
              rownames_to_column("cell") |> 
              select(cell, acc_temp, cell_type))

MEs |>
  pivot_longer(rownames(DMEs), 
               names_to = "module",
               values_to = "eigengene") |> 
  mutate(module = factor(module, levels = rownames(DMEs))) |> 
  ggplot() +
  geom_density(aes(x = eigengene,
                   fill = acc_temp),
               alpha = 0.5) +
  scale_x_continuous(name = "Eigengene",
                     expand = c(0, 0)) +
  scale_y_continuous(name = "Density",
                     expand = c(0, 0)) +
  scale_fill_manual(values = acc_colors) +
  facet_wrap(~ module) +
  theme_minimal_hgrid() +
  theme(legend.position = "none")

# Fetch data
FetchData(object = dat, 
          vars = red[1:20],
          slot = "data") |>
  rownames_to_column("cell") |> 
  pivot_longer(-cell, names_to = "feature", values_to = "data") |> 
  full_join(meta) |> 
  ggplot() +
  geom_density(aes(x = data,
                   fill = acc_temp),
               alpha = 0.5) +
  scale_x_continuous(name = "Feature count",
                     expand = c(0, 0)) +
  scale_y_continuous(name = "Density",
                     expand = c(0, 0)) +
  scale_fill_manual(values = acc_colors) +
  facet_wrap(~feature) +
  theme_minimal_hgrid() +
  theme(legend.position = "none")


# Reorder brown by highest expression
brown <- HVFInfo(dat) |> 
  rownames_to_column("gene") |> 
  filter(gene %in% brown) |> 
  arrange(desc(gmean)) |>
  pull("gene")

# Top 20 highest expressed
FetchData(object = dat, 
          vars = "ATPsynE",
          slot = "data") |>
  rownames_to_column("cell") |> 
  pivot_longer(-cell, names_to = "feature", values_to = "data") |> 
  mutate(feature = factor(feature, levels = brown)) |> 
  full_join(dat@meta.data |> rownames_to_column("cell")) |> 
  ggplot() +
  geom_density(aes(x = data,
                   fill = acc_temp,
                   color = acc_temp),
               alpha = 0.5) +
  scale_x_continuous(name = "Log-normalized expression",
                     expand = c(0, 0)) +
  scale_y_continuous(name = "Density",
                     expand = c(0, 0)) +
  scale_fill_manual(values = acc_colors) +
  scale_color_manual(values = acc_colors) +
  # facet_wrap(~feature,
  #            scales = "free_y") +
  theme_minimal_hgrid() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey95"))

# Save
ggsave(here::here(fig_dir, "brown", "brown_top_301-346.png"),
       height = 60,
       width = 60,
       bg = "white",
       limitsize = FALSE,
       units = "cm")

# Reorder red by highest expression
red <- HVFInfo(dat) |> 
  rownames_to_column("gene") |> 
  filter(gene %in% red) |> 
  arrange(desc(gmean)) |>
  pull("gene")

# Top 20 highest expressed
p <- FetchData(object = dat, 
          vars = red,
          slot = "data") |>
  rownames_to_column("cell") |> 
  pivot_longer(-cell, names_to = "feature", values_to = "data") |> 
  mutate(feature = factor(feature, levels = red)) |> 
  full_join(meta) |> 
  ggplot() +
  geom_density(aes(x = data,
                   fill = acc_temp,
                   color = acc_temp),
               alpha = 0.5) +
  scale_x_continuous(name = "Log-normalized expression",
                     expand = c(0, 0)) +
  scale_y_continuous(name = "Density",
                     expand = c(0, 0)) +
  scale_fill_manual(values = acc_colors) +
  scale_color_manual(values = acc_colors) +
  facet_wrap(~feature,
             scales = "free_y") +
  theme_minimal_hgrid() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey95"))

# Save
ggsave(here::here(fig_dir, "red", "red_all.png"),
       height = 100,
       width = 100,
       bg = "white",
       limitsize = FALSE,
       units = "cm")

meta <- dat@meta.data |> rownames_to_column("cell")
# Top hub genes for brown
FetchData(object = dat, 
          vars = c("Non2", "Fkbp39", "Df31", "Set"),
          slot = "data") |>
  rownames_to_column("cell") |> 
  pivot_longer(-cell, names_to = "feature", values_to = "data") |> 
  mutate(feature = factor(feature, levels = brown)) |> 
  full_join(meta) |> 
  ggplot() +
  geom_density(aes(x = data,
                   fill = acc_temp,
                   color = acc_temp),
               alpha = 0.5) +
  scale_x_continuous(name = "Log-normalized expression",
                     expand = c(0, 0)) +
  scale_y_continuous(name = "Density",
                     expand = c(0, 0)) +
  scale_fill_manual(values = acc_colors) +
  scale_color_manual(values = acc_colors) +
  facet_wrap(~feature,
             scales = "free_y") +
  theme_minimal_hgrid() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey95"))

# Save
ggsave(here::here(fig_dir, "brown", "brown_top_kMEs.png"),
       height = 20,
       width = 25,
       bg = "white",
       limitsize = FALSE,
       units = "cm")

# Top hub genes for red
FetchData(object = dat, 
          vars = c("RpS20", "RpL18A", "RpL13A", "RpL27A"),
          slot = "data") |>
  rownames_to_column("cell") |> 
  pivot_longer(-cell, names_to = "feature", values_to = "data") |> 
  mutate(feature = factor(feature, levels = red)) |> 
  full_join(meta) |> 
  ggplot() +
  geom_density(aes(x = data,
                   fill = acc_temp,
                   color = acc_temp),
               alpha = 0.5) +
  scale_x_continuous(name = "Log-normalized expression",
                     expand = c(0, 0)) +
  scale_y_continuous(name = "Density",
                     expand = c(0, 0)) +
  scale_fill_manual(values = acc_colors) +
  scale_color_manual(values = acc_colors) +
  facet_wrap(~feature,
             scales = "free_y") +
  theme_minimal_hgrid() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey95"))

# Save
ggsave(here::here(fig_dir, "red", "red_top_kMEs.png"),
       height = 20,
       width = 25,
       bg = "white",
       limitsize = FALSE,
       units = "cm")

# ATAC -------------------------------------------------------------------------
# Run a WGCNA on ATAC data per Brent's suggestion...

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/13_dat_wgcna_atac.rds"))
DMEs <- readRDS(here::here("output/wgcna/atac/DMEs.rds"))
blue <- read_csv(here::here("output/wgcna/atac/blue.csv")) |> pull()
brown <- read_csv(here::here("output/wgcna/atac/brown.csv")) |> pull()

# Output fig directory
fig_dir <- here::here("output/figs/wgcna/atac")

# Plot the soft power analysis results
PlotSoftPowers(dat) |> 
  wrap_plots(ncol = 2)

# Save
ggsave(here::here(fig_dir, "softpower_threshold.png"),
       height = 20,
       width = 20,
       units = "cm")

# Plot the dendrogram for each
png(here::here(fig_dir, "dendrogram.png"),
    res = 300,
    height = 15,
    width = 20,
    units = "cm")
PlotDendrogram(
  dat, 
  main = "hdWGCNA Dendrogram"
)
dev.off()

mod_tally <- GetModules(dat) |> 
  group_by(module) |> 
  tally()

x <- DMEs |>
  rownames_to_column("mod") |> 
  mutate(mod = fct_relevel(mod, DMEs |> arrange(avg_log2FC) |> pull(module)))

mod_colors <- c(
  "blue", 
  "magenta",
  "green",
  "red",
  "purple",
  "yellow", 
  "pink",
  "black", 
  "turquoise", 
  "brown"
) |> colorspace::desaturate(amount = 0.5)

left_join(x, mod_tally) |> 
  ggplot(aes(y = mod,
             x = avg_log2FC)) +
  geom_segment(aes(xend = 0, 
                   yend = mod,
                   color = mod),
               linewidth = 1) +
  annotate(geom = "text", 
           y = "blue", 
           x = 0.18, 
           color = "grey20",
           label = "*", 
           vjust = 0.8, 
           size = 10) +
  annotate(geom = "text", 
           y = "brown", 
           x = -0.18, 
           color = "grey20",
           label = "*", 
           vjust = 0.8, 
           size = 10) +
  geom_point(aes(fill = mod, 
                 size = n), 
             color = "grey80",
             shape = 21) +
  geom_vline(xintercept = 0, color = "grey50") +
  scale_x_continuous(name = "log2(fold-change)",
                     limits = c(-0.2, 0.2)) +
  scale_fill_manual(values = rev(mod_colors)) +
  scale_color_manual(values = rev(mod_colors)) +
  scale_size_continuous(range = c(5, 10)) +
  labs(y = element_blank()) +
  cowplot::theme_minimal_vgrid() +
  theme(legend.position = "none",
        panel.grid.major.y = element_line(color = "grey95"),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

# Save
ggsave(here::here(fig_dir, "DME_lolipops.png"),
       height = 20,
       width = 20,
       units = "cm")


# Make featureplot of MEs for each module
ModuleFeaturePlot(
  dat,
  features = "MEs",
  order = TRUE) |> 
  patchwork::wrap_plots(ncol = 3)

# Save
ggsave(here::here(fig_dir, "module_feature_plot.png"),
       height = 20,
       width = 20,
       units = "cm")

# Module correlogram
png(here::here(fig_dir, "module_correlogram.png"),
    res = 300,
    height = 20,
    width = 20,
    units = "cm")
ModuleCorrelogram(dat)
dev.off()

## Messing with plot types for individual genes etc ----------------------------


# Reorder red by highest expression
red <- HVFInfo(dat) |> 
  rownames_to_column("gene") |> 
  filter(gene %in% red) |> 
  arrange(desc(gmean)) |>
  pull("gene")

# Top 20 highest expressed
p <- FetchData(object = dat, 
               vars = red,
               slot = "data") |>
  rownames_to_column("cell") |> 
  pivot_longer(-cell, names_to = "feature", values_to = "data") |> 
  mutate(feature = factor(feature, levels = red)) |> 
  full_join(meta) |> 
  ggplot() +
  geom_density(aes(x = data,
                   fill = acc_temp,
                   color = acc_temp),
               alpha = 0.5) +
  scale_x_continuous(name = "Log-normalized expression",
                     expand = c(0, 0)) +
  scale_y_continuous(name = "Density",
                     expand = c(0, 0)) +
  scale_fill_manual(values = acc_colors) +
  scale_color_manual(values = acc_colors) +
  facet_wrap(~feature,
             scales = "free_y") +
  theme_minimal_hgrid() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey95"))

# Save
ggsave(here::here(fig_dir, "red_all.png"),
       height = 100,
       width = 100,
       bg = "white",
       limitsize = FALSE,
       units = "cm")

# Module eigengene values ----

# Per-cell MEs for each group
MEs <- GetMEs(dat) |>
  rownames_to_column("cell") |>
  full_join(dat@meta.data |>
              rownames_to_column("cell") |>
              select(cell, acc_temp, cell_type))

MEs |>
  pivot_longer(rownames(DMEs),
               names_to = "module",
               values_to = "eigengene") |>
  mutate(module = factor(module, levels = rownames(DMEs))) |>
  ggplot() +
  geom_density(aes(x = eigengene,
                   fill = acc_temp),
               alpha = 0.5) +
  scale_x_continuous(name = "Eigengene",
                     expand = c(0, 0)) +
  scale_y_continuous(name = "Density",
                     expand = c(0, 0)) +
  scale_fill_manual(values = acc_colors) +
  facet_wrap(~ module) +
  theme_minimal_hgrid() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey95"))

# Save
ggsave(here::here(fig_dir, "module_eigengenes.png"),
       height = 20,
       width = 30,
       bg = "white",
       limitsize = FALSE,
       units = "cm")

# Random 20 turquoise
FetchData(object = dat, 
          vars = turquoise[1:20],
          slot = "data") |>
  rownames_to_column("cell") |> 
  pivot_longer(-cell, names_to = "feature", values_to = "data") |> 
  full_join(meta) |> 
  ggplot() +
  geom_density(aes(x = data,
                   fill = acc_temp,
                   color = acc_temp),
               alpha = 0.5) +
  scale_x_continuous(name = "Feature count",
                     expand = c(0, 0)) +
  scale_y_continuous(name = "Density",
                     expand = c(0, 0)) +
  scale_fill_manual(values = acc_colors) +
  scale_color_manual(values = acc_colors) +
  facet_wrap(~feature) +
  theme_minimal_hgrid() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey95"))

# Save
ggsave(here::here(fig_dir, "turquoise_20.png"),
       height = 20,
       width = 30,
       bg = "white",
       limitsize = FALSE,
       units = "cm")


# Random 20 brown
FetchData(object = dat, 
          vars = brown[1:20],
          slot = "data") |>
  rownames_to_column("cell") |> 
  pivot_longer(-cell, names_to = "feature", values_to = "data") |> 
  full_join(meta) |> 
  ggplot() +
  geom_density(aes(x = data,
                   fill = acc_temp,
                   color = acc_temp),
               alpha = 0.5) +
  scale_x_continuous(name = "Feature count",
                     expand = c(0, 0)) +
  scale_y_continuous(name = "Density",
                     expand = c(0, 0)) +
  scale_fill_manual(values = acc_colors) +
  scale_color_manual(values = acc_colors) +
  facet_wrap(~feature) +
  theme_minimal_hgrid() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey95"))

# Save
ggsave(here::here(fig_dir, "brown_20.png"),
       height = 20,
       width = 30,
       bg = "white",
       limitsize = FALSE,
       units = "cm")
