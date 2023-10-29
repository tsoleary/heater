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

# Plot DME lolipop
PlotDMEsLollipop(
  dat,
  DMEs,
  wgcna_name = "wgcna",
  pvalue = "p_val_adj"
)

mod_tally <- GetModules(dat) |> 
  group_by(module) |> 
  tally()

x <- DMEs |>
  rownames_to_column("mod") |> 
  mutate(mod = fct_relevel(mod, DMEs |> arrange(avg_log2FC) |> pull(module)))

left_join(x, mod_tally) |> 
  ggplot(aes(y = mod,
             x = avg_log2FC, 
             alpha = p_val_adj < 0.05)) +
  geom_segment(aes(xend = 0, 
                   yend = mod,
                   color = mod),
               size = 0.5) +
  geom_point(aes(fill = mod, 
                 size = n), 
             color = "grey80",
             shape = 21) +
  geom_vline(xintercept = 0, color = "grey50") +
  scale_x_continuous(name = "log2(fold-change)") +
  scale_fill_manual(values = rev(x$module)) +
  scale_alpha_discrete(range = c(0.3, 1)) +
  scale_color_manual(values = rev(x$module)) +
  labs(y = element_blank()) +
  cowplot::theme_minimal_grid() +
  theme(legend.position = "none",
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
  wrap_plots(ncol = 3)

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
          vars = brown[301:346],
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

# ATAC -------------------------------------------------------------------------
# Run a WGCNA on ATAC data per Brent's suggestion...

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/13_dat_wgcna_atac.rds"))
DMEs <- readRDS(here::here("output/wgcna/atac/DMEs.rds"))
turquoise <- read_csv(here::here("output/wgcna/atac/turquoise.csv")) |> pull()
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


# Plot DME lolipop
PlotDMEsLollipop(
  dat, 
  DMEs, 
  wgcna_name = "wgcna", 
  pvalue = "p_val_adj"
)

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
  wrap_plots(ncol = 3)

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
