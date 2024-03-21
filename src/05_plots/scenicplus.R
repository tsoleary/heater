# ------------------------------------------------------------------------------
# SCENIC+ plots generated in R
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(ggh4x)
library(Seurat)
library(Signac)

# Output dir
fig_dir <- here::here("output/figs/scenicplus")

# Load data
dotplot_dat <- read_csv(
  here::here("data/processed/scenic/dotplot.csv")) |> 
  mutate(eRegulon_type = factor(str_extract(eRegulon_name, "(?<=_)\\+|-(?=_)"), levels = c("+", "-"))) |> 
  mutate(eRegulon_name = str_replace_all(eRegulon_name, "_[\\+|-]_", " "))

# Create eRegulon order based on TF expression per cell-type
eRegulon_order <- dotplot_dat |> 
  group_by(eRegulon_name, eRegulon_type) |> 
  distinct() |>
  slice_max(color_val, n = 1) |> 
  arrange(eRegulon_type, index) |> 
  pull(eRegulon_name)

# Make the plot for heatmap
dotplot_dat |> 
  distinct() |> 
  filter(eRegulon_type == "+") |> 
  separate(index, into = c("cell_type", "acclimation"), 
           sep = " (?=[[:digit:]])") |>
  mutate(cell_type = ordered(cell_type)) |> 
  mutate(eRegulon_name = factor(eRegulon_name, 
                                levels = eRegulon_order)) |> 
  ggplot(aes(y = interaction(acclimation, cell_type, sep = " "), 
             x = eRegulon_name,
             fill = color_val)) +
  geom_point(aes(size = size_val),
             shape = 21,
             stroke = 0.1) +
  scale_fill_gradient(low = "grey95", 
                      # mid = "grey85", 
                       high = "forestgreen",
                       #midpoint = 0.5,
                      name = "TF exp") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.98),
        axis.title = element_blank(),
        strip.background = element_rect(fill = "grey90", color = "grey90")) +
  scale_y_discrete(labels = rep(c("18°C", "25°C"), 9)) +
  facet_grid(rows = vars(cell_type),
             cols = vars(eRegulon_type),
             switch = "y",
             scales = "free") +
  force_panelsizes(cols = c(1, 0.4)) +
  scale_size_continuous(range = c(0.25, 4),
                        name = "Target region\nenrichment") +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, 
                                         hjust = 1, 
                                         margin = margin(r = 5)),
        strip.background.y = element_rect(fill ="grey98", color = "grey98"))

ggsave(here::here(fig_dir, "eRegulon_heatmap_dots_activator.png"),
       bg = "white",
       width = 10,
       height = 6,
       units = "in") 

dotplot_dat |> 
  filter(TF == "Mef2")


scplus_region <- read_csv(here::here("data/processed/scenic/scplus_region.csv")) |> 
  mutate(Cell = str_remove_all(Cell, "___cisTopic"))

dat <- readRDS(here::here("data/processed/seurat_object/10_dat_linked.rds"))

 

dat@meta.data <- dat@meta.data |> 
  rownames_to_column("Cell") |> 
  full_join(scplus_region) |> 
  column_to_rownames("Cell")


FeaturePlot(
  dat,
  "Six4_extended_+_(37r)",
  cols = c("grey85", "darkgreen"),
  split.by = "acc_temp") &
  theme_void() &
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold")) 

ggsave(here::here(fig_dir, "Six4_feature.png"),
       bg = "white",
       width = 10,
       height = 5,
       units = "in")


FeaturePlot(
  dat,
  "Mef2_+_(56r)",
  cols = c("grey85", "darkgreen"),
  split.by = "acc_temp") &
  theme_void() &
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold")) 

ggsave(here::here(fig_dir, "Mef2_feature.png"),
       bg = "white",
       width = 10,
       height = 5,
       units = "in")


