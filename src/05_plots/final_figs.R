# ------------------------------------------------------------------------------
# Final figures
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)
library(cowplot)

# Output fig dir
fig_dir <- "output/figs/final"

# Load plot themes & modified functions
source(here::here("src/05_plots/00_plot_themes.R"))
source(here::here("src/05_plots/plot_functions.R"))

# Phenotype data
pheno <- read_csv(
  here::here("data/raw/pheno/acc_hs_survival.csv")) |> 
  mutate(acc_temp_factor = factor(paste0(acc_temp, "°C")))

# Load raw data
dat_raw <- readRDS(
  here::here("data/processed/seurat_object/00_dat_raw.rds")
)

# Load final seurat object
dat <- readRDS(
  here::here("data/processed/seurat_object/14_dat_tf_motif.rds")
)

# Recode sample name for prettier plotting
dat@meta.data <- dat@meta.data |> 
  mutate(sample_name_fig = case_when(sample_name == "18C_Rep1" ~ "18°C Rep 1",
                                     sample_name == "18C_Rep2" ~ "18°C Rep 2",
                                     sample_name == "25C_Rep1" ~ "25°C Rep 1",
                                     sample_name == "25C_Rep2" ~ "25°C Rep 2")) |> 
  mutate(sample_name_fig = factor(sample_name_fig, levels = c(c("25°C Rep 2",
                                                        "25°C Rep 1",
                                                        "18°C Rep 2",
                                                        "18°C Rep 1"))))
# Create dataframe with cell names as column for plotting
meta <- dat@meta.data |> 
  rownames_to_column("cell")

# Load all cell-type specific acclimation comparisons
degs <- readRDS(here::here("output/degs/degs_cell-type_MAST.rds")) |> 
  mutate(sig = ifelse(abs(avg_log2FC_18_25) >= 0.25 & p_val_adj < 0.05 &
                        (pct.18 >= 0.10 | pct.25 >= 0.10), 
         TRUE, FALSE)) 
dars <- readRDS(here::here("output/dars/dars_cell-type_MAST.rds")) |> 
  mutate(sig = ifelse(abs(avg_log2FC_18_25) >= 0.25 & p_val_adj < 0.05 &
                        (pct.18 >= 0.10 | pct.25 >= 0.10), 
                      TRUE, FALSE)) 

# eRegulon data
eReg <- read_csv(here::here("data/processed/scenic/eRegulon_data.csv"))

# ------------------------------------------------------------------------------
# Manuscript figures -----------------------------------------------------------
# ------------------------------------------------------------------------------

# Figure 1 ---------------------------------------------------------------------

# Acclimation design schematic
p_exp <- ggdraw() + 
  draw_image(here::here("output/figs/exp_design/acc_design.png"))

# Acclimation egg hatching phenotype
p_pheno <- pheno |>
  ggplot(aes(x = acc_temp,
             y = n_hatched/n_eggs*100,
             group = acc_temp_factor)) +
  geom_boxplot(width = 1.5,
               color = "grey20",
               outlier.shape = NA) +
  ggbeeswarm::geom_beeswarm(size = 3,
                            cex = 2,
                            shape = 21,
                            fill = "grey50",
                            color = "grey20",
                            alpha = 0.9) +
  xlab("Acclimation temperature") +
  scale_y_continuous(name = "Hatching success\nafter 45 min @ 39°C",
                     limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.05)),
                     labels = function(x) paste0(x, "%")) + 
  scale_x_continuous(breaks = c(18, 25, 30),
                     labels = levels(pheno$acc_temp_factor)) +
  cowplot::theme_minimal_hgrid()

# Knee plots together
cowplot::plot_grid(
  p_exp, 
  p_pheno,
  labels = c("A", "B"),
  rel_heights = c(1, 2),
  nrow = 2
)

# Save
ggsave(here::here(fig_dir, "fig_1.png"),
       height = 15,
       width = 15,
       units = "cm")
ggsave(here::here(fig_dir, "fig_1.pdf"),
       height = 15,
       width = 15,
       units = "cm")

# Figure 2 ---------------------------------------------------------------------

# UMAP acclimation treatment overlay
p_overlay <- DimPlot(
  dat,
  group.by = "acc_temp",
  reduction = "umap") +
  scale_color_manual(name = element_blank(),
                     values = acc_colors) +
  labs(title = element_blank()) +
  theme_void() +
  theme(legend.position = c(0.8, 0.2)) +
  guides(color = guide_legend(byrow = TRUE,
                              nrow = 2,
                              override.aes = list(size = 3)))

# UMAP clusters
p_clust <- DimPlot(
  dat,
  group.by = "seurat_clusters",
  reduction = "umap") +
  labs(title = element_blank()) +
  scale_color_manual(
    values = scales::hue_pal()(12) |> 
      colorspace::desaturate(amount = 0.5)
  ) +
  theme_void() +
  theme(legend.position = "none") +
  guides(color = guide_legend(byrow = TRUE,
                              nrow = 1,
                              override.aes = list(size = 3)))
  
p_clust_labeled <- LabelClusters(
  p_clust, 
  id = "seurat_clusters",
  fontface = "bold",
  color = "grey20",
  alpha = 0.8,
  box = TRUE,
  repel = TRUE,
  family = "Arial"
)


# UMAP cell_type
# Set cell-type as the primary identity
Idents(dat) <- "cell_type"

# UMAP colored by cell-types
p_cell_type <- DimPlot(dat, 
                reduction = "umap",
                repel = TRUE,
                label.color = "grey20",
                pt.size = 0.2) +
  theme_void() +
  theme(legend.position = "none",
        strip.text = element_text(size = 16, face = "bold"),
        strip.background = element_rect(fill = "grey95", color = "grey95")) +
  scale_color_manual(
    values = c("grey90",
               "#ADD9F4",
               "#57A4B2",
               "#D39C9E",
               "#FEF29A",
               "#F9DCEE",
               "#819FC5",
               "#A7BF9B",
               "#bfa3a4")
  )

p_cell_type_label <- LabelClusters(plot = p_cell_type, 
              id = "ident",
              box = TRUE,
              size = 3,
              position = "nearest",
              label.padding = unit(0.1, "lines"),
              alpha = 0.75,
              color = c("grey10"))


# Combining the top row
f2_row1 <- cowplot::plot_grid(
  p_overlay,
  p_clust_labeled,
  p_cell_type_label,
  nrow = 1,
  labels = c("A", "B", "C"),
  label_colour = "grey20"
)

# Marker gene dot plot
markers <- read_csv(here::here("output/markers/cluster_markers_all.csv"))

# Top three markers for each cluster
top_markers <- markers |> 
  group_by(cluster) |> 
  slice_min(max_pval, n = 3, with_ties = FALSE) |> 
  mutate(cluster = factor(cluster, levels = as.character(0:11))) |> 
  arrange(cluster)

# Set default assay for fetching data
DefaultAssay(dat) <- "SCT"
Idents(dat) <- "seurat_clusters"

# Percent nuclei with expressed gene
pct_exp <- FetchData(dat, 
                     vars = c(top_markers$gene, "seurat_clusters"), 
                     layer = "counts") |> 
  pivot_longer(cols = -c("seurat_clusters"),
               names_to = "gene",
               values_to = "count") |> 
  mutate(exp = count > 0) |> 
  group_by(seurat_clusters, gene) |> 
  summarize(pct_exp = sum(exp)/n()*100)

# Average expression 
avg_exp <- AverageExpression(dat)$SCT |> 
  as.data.frame() |> 
  rownames_to_column("gene") |> 
  filter(gene %in% top_markers$gene) |> 
  pivot_longer(cols = -c(gene),
               names_to = "cluster",
               values_to = "expression") |> 
  mutate(cluster = str_remove_all(cluster, "g")) |> 
  group_by(gene) |> 
  mutate(norm_expression = scale(expression)) 

# Combine the average expression with percent expressed data
df <- full_join(
  avg_exp, 
  pct_exp,
  by = c("cluster" = "seurat_clusters", "gene")
)

# Get a pleasing gene order for plotting
gene_order <- df |> 
  mutate(cluster = as.numeric(cluster)) |> 
  group_by(cluster) |> 
  arrange(cluster, desc(norm_expression)) |> 
  slice_max(norm_expression, n = 3) |> 
  pull(gene)

# Marker gene dot plot
p_mark_dot <- df |> 
  mutate(cluster = factor(cluster, levels = as.character(0:11))) |> 
  mutate(gene = factor(gene, levels = top_markers$gene)) |> 
  ggplot(aes(x = gene, 
             y = cluster, 
             fill = norm_expression,
             size = pct_exp)) +
  geom_point(shape = 21, stroke = 0.25) +
  scale_fill_gradient(low = "grey95",
                      high = "firebrick",
                      name = "Normalized\nExpression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.98),
        axis.title = element_blank(),
        strip.background = element_rect(fill = "grey90", color = "grey90")) +
  ggh4x::force_panelsizes(cols = c(1, 0.4)) +
  scale_size_continuous(range = c(0.25, 6),
                        labels = scales::label_percent(scale = 1),
                        name = "Nuclei\nexpressing\nthe gene")

f2 <- cowplot::plot_grid(
  f2_row1,
  p_mark_dot,
  nrow = 2,
  labels = c("", "D"),
  label_colour = "grey20",
  label_y = 1.05,
  rel_heights = c(1, 1)
)

# Save
ggsave(here::here(fig_dir, "fig_2.png"),
       plot = f2,
       height = 20,
       width = 27,
       units = "cm")
ggsave2(here::here(fig_dir, "fig_2.pdf"),
       plot = f2,
       height = 20,
       width = 27,
       units = "cm")

# Figure 3 ---------------------------------------------------------------------

# Set cell types as identity to plot across
Idents(dat) <- "cell_type"
DefaultAssay(dat) <- "peaks"

# Btl coverage plot -----

# Set gene for coverage plot
gene <- "btl"

# Accessibility peaks
p_cov <- CoveragePlot(
  object = dat,
  region = gene,
  features = gene,
  links = FALSE,
  peaks = FALSE,
  annotation = FALSE,
  extend.upstream = 0,
  extend.downstream = 0) &
  scale_fill_manual(values = cell_type_colors) &
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 6),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) 

# Expression plot
p_exp <- ExpressionPlot(
  dat,
  assay = "RNA",
  features = gene) &
  scale_fill_manual(values = cell_type_colors) &
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6))

# Gene annotation plot with boxes for gene
p_ann <- AnnotationPlot_TSO(
  dat,
  region = gene,
  size_gene_label = 3,
  extend.upstream = 0,
  extend.downstream = 0)  &
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank())

# Combine together
p_btl <- CombineTracks(
  list(p_cov, p_ann),
  expression.plot = p_exp,
  widths = c(9, 1),
  heights = c(5, 1)) &
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 8))

# Mef2 coverage plot -----

# Set gene for coverage plot
gene <- "Mef2"

# Accessibility peaks
p_cov <- CoveragePlot(
  object = dat,
  region = gene,
  features = gene,
  links = FALSE,
  peaks = FALSE,
  annotation = FALSE,
  extend.upstream = 0,
  extend.downstream = 0) &
  scale_fill_manual(values = cell_type_colors) &
  theme(legend.position = "none",
        strip.text = element_blank(),
        strip.text.y.left = element_blank(),
        strip.background = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) 

# Expression plot
p_exp <- ExpressionPlot(
  dat, 
  assay = "RNA",
  features = gene) &
  scale_fill_manual(values = cell_type_colors) &
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6))

# Gene annotation plot with boxes for gene
p_ann <- AnnotationPlot_TSO(
  dat,
  region = gene,
  size_gene_label = 3,
  extend.upstream = 0,
  extend.downstream = 0)  &
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank())

# Combine together
p_Mef2 <- CombineTracks(
  list(p_cov, p_ann),
  expression.plot = p_exp,
  widths = c(9, 1),
  heights = c(5, 1)) &
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 8))

# Mef2 coverage plot -----

# Set gene for coverage plot
gene <- "otp"

# Accessibility peaks
p_cov <- CoveragePlot(
  object = dat,
  region = gene,
  features = gene,
  links = FALSE,
  peaks = FALSE,
  annotation = FALSE,
  extend.upstream = 0,
  extend.downstream = 0) &
  scale_fill_manual(values = cell_type_colors) &
  theme(legend.position = "none",
        strip.text = element_blank(),
        strip.text.y.left = element_blank(),
        strip.background = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) 

# Expression plot
p_exp <- ExpressionPlot(
  dat, 
  assay = "RNA",
  features = gene) &
  scale_fill_manual(values = cell_type_colors) &
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6))

# Gene annotation plot with boxes for gene
p_ann <- AnnotationPlot_TSO(
  dat,
  region = gene,
  size_gene_label = 3,
  extend.upstream = 0,
  extend.downstream = 0)  &
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank())

# Combine together
p_otp <- CombineTracks(
  list(p_cov, p_ann),
  expression.plot = p_exp,
  widths = c(9, 1),
  heights = c(5, 1)) &
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 8))

# Get two genes combined together
p_top <- cowplot::plot_grid(
  p_btl, 
  p_Mef2,
  p_otp,
  rel_widths = c(1, 0.7, 0.7),
  nrow = 1
)

ggsave(here::here(fig_dir, "f3_top_test.png"),
       plot = p_top,
       height = 8,
       width = 30,
       units = "cm")
ggsave(here::here(fig_dir, "f3_top_test.pdf"),
       plot = p_top,
       height = 8,
       width = 30,
       units = "cm")


#### WAIT CAREFULL!!! Make sure the annotation plot lines up right....
# GetLinkedPeaks(dat, features = "btl")
# GetLinkedPeaks(dat, features = "Mef2")
# GetLinkedPeaks(dat, features = "otp")

# I think the move is to include the feature plots
p_btl_umap <- FeaturePlot(
  dat,
  features = "3L-14076419-14077188", # btl peak
  cols = c("grey95", "darkgreen")) & 
  theme_void() &
  theme(title = element_blank(),
        legend.position = c(0.6, 0.2),
        legend.key.size = unit(0.6, "lines"))
  

# I think the move is to include the feature plots
p_Mef2_umap <- FeaturePlot(
  dat,
  features = "2R-9929424-9930308", # Mef2 peak
  cols = c("grey95", "darkgreen")) & 
  theme_void() &
  theme(title = element_blank(),
        legend.position = c(0.6, 0.2),
        legend.key.size = unit(0.6, "lines"))

# I think the move is to include the feature plots
p_otp_umap <- FeaturePlot(
  dat,
  features = "2R-20897593-20898203", # Mef2 peak
  cols = c("grey95", "darkgreen")) & 
  theme_void() &
  theme(title = element_blank(),
        legend.position = c(0.6, 0.2),
        legend.key.size = unit(0.6, "lines"))

# Bottom
p_bottom <- cowplot::plot_grid(
  NULL, p_btl_umap, p_Mef2_umap, p_otp_umap,
  nrow = 1,
  rel_widths = c(0.2, 0.5, 0.5, 0.5)
)


# Combine all together -----
f3 <- cowplot::plot_grid(
  p_top, p_bottom, 
  rel_heights = c(1, 0.6),
  labels = c("A", "B"),
  nrow = 2
)

# Save
ggsave(here::here(fig_dir, "fig_3.png"),
       plot = f3,
       bg = "white",
       height = 20,
       width = 25,
       units = "cm")
ggsave(here::here(fig_dir, "fig_3.pdf"),
       plot = f3,
       bg = "white",
       height = 20,
       width = 25,
       units = "cm")

# Figure 4 ---------------------------------------------------------------------

# Load data
dotplot_dat <- read_csv(
  here::here("data/processed/scenic/dotplot.csv")) |> 
  mutate(eRegulon_type = factor(
    str_extract(eRegulon_name, "(?<=_)\\+|-(?=_)"), levels = c("+", "-"))) |> 
  mutate(Region_signature_name = eRegulon_name) |> 
  mutate(eRegulon_name = str_replace_all(eRegulon_name, "_[\\+|-]_", " ")) |> 
  mutate(eRegulon_type = ifelse(eRegulon_type == "+", "activator", "repressor"))

# Create eRegulon order based on TF expression per cell-type
eRegulon_order <- dotplot_dat |> 
  group_by(eRegulon_name, eRegulon_type) |> 
  distinct() |>
  slice_max(color_val, n = 1) |> 
  arrange(eRegulon_type, index) |> 
  pull(eRegulon_name)

degs |> 
  filter(gene %in% dotplot_dat$TF)

# Make the plot for heatmap
p_ereg <- dotplot_dat |> 
  distinct() |> 
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
             scales = "free",
             labeller = label_wrap_gen(multi_line = TRUE)) +
  ggh4x::force_panelsizes(cols = c(1, 0.4)) +
  scale_size_continuous(range = c(0.25, 6),
                        name = "Target region\nenrichment") +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, 
                                         hjust = 1, 
                                         margin = margin(r = 5)),
        strip.background.y = element_rect(fill ="grey98", color = "grey98"))


# Create density plots with rug plots beneath for diff. expressed TFs 
DefaultAssay(dat) <- "RNA"

p_blimp_leg <- FetchData(object = dat, 
          vars = "Blimp-1",
          layer = "data") |>
  rownames_to_column("cell") |> 
  pivot_longer(-cell, names_to = "feature", values_to = "data") |> 
  full_join(meta) |> 
  mutate(rug_position = ifelse(acc_temp == "18°C", -0.07, -0.17)) |> 
  filter(cell_type == "foregut & hindgut prim.") |> 
  ggplot() +
  geom_label(aes(x = 2.2, y = 1.2),
             label = "Blimp-1",
             fill = "#F9DCEE",
             color = "grey20",
             size = 5) +
  geom_density(aes(x = data,
                   fill = acc_temp,
                   color = acc_temp),
               alpha = 0.5) +
  geom_point(aes(y = rug_position, 
                 x = data,
                 color = acc_temp),
             shape = "|", 
             size = 2) +
  scale_x_continuous(name = "",
                     expand = c(0, 0.05)) +
  scale_y_continuous(name = "Density",
                     expand = c(0, 0.03),
                     limits = c(-0.2, 1.4)) +
  scale_fill_manual(values = acc_colors,
                    name = "") +
  scale_color_manual(values = acc_colors,
                     name = "") +
  coord_cartesian(clip = "off") +
  theme_cowplot() +
  theme(axis.text = element_text(color = "grey40"),
        axis.ticks = element_line(color = "grey40"),
        axis.line = element_line(color = "grey40"),
        axis.title = element_text(color = "grey40"))

p_blimp <- p_blimp_leg + theme(legend.position = "none")

legend_density <- ggdraw(get_legend(p_blimp_leg))

p_nub <- FetchData(object = dat, 
                     vars = "nub",
                     layer = "data") |>
  rownames_to_column("cell") |> 
  pivot_longer(-cell, names_to = "feature", values_to = "data") |> 
  full_join(meta) |> 
  mutate(rug_position = ifelse(acc_temp == "18°C", -0.2, -0.55)) |> 
  filter(cell_type == "peripheral nervous system prim.") |> 
  ggplot() +
  geom_label(aes(x = 2.2, y = 3.8),
             label = "nub",
             fill = "#ADD9F4",
             color = "grey20",
             size = 5) +
  geom_density(aes(x = data,
                   fill = acc_temp,
                   color = acc_temp),
               alpha = 0.5) +
  geom_point(aes(y = rug_position, 
                 x = data,
                 color = acc_temp),
             shape = "|", 
             size = 2) +
  scale_x_continuous(name = "",
                     expand = c(0, 0.05)) +
  scale_y_continuous(name = "",
                     expand = c(0, 0.05),
                     limits = c(-0.7, 4.5)) +
  scale_fill_manual(values = acc_colors) +
  scale_color_manual(values = acc_colors) +
  coord_cartesian(clip = "off") +
  theme_cowplot() +
  theme(legend.position = "none",
        axis.text = element_text(color = "grey40"),
        axis.ticks = element_line(color = "grey40"),
        axis.line = element_line(color = "grey40"))

# p_hth <- FetchData(object = dat, 
#                      vars = "hth",
#                      layer = "data") |>
#   rownames_to_column("cell") |> 
#   pivot_longer(-cell, names_to = "feature", values_to = "data") |> 
#   full_join(meta) |> 
#   mutate(rug_position = ifelse(acc_temp == "18°C", -0.0325, -0.0825)) |> 
#   filter(cell_type == "mesoderm prim.") |> 
#   ggplot() +
#   geom_label(aes(x = 3, y = 0.6),
#              label = "hth",
#              fill = "#D39C9E",
#              color = "grey20",
#              size = 5) +
#   geom_density(aes(x = data,
#                    fill = acc_temp,
#                    color = acc_temp),
#                alpha = 0.5) +
#   geom_point(aes(y = rug_position, 
#                  x = data,
#                  color = acc_temp),
#              shape = "|", 
#              size = 2) +
#   scale_x_continuous(name = "Log-normalized expression",
#                      expand = c(0, 0.05)) +
#   scale_y_continuous(name = "",
#                      expand = c(0, 0.015),
#                      limits = c(-0.10, 0.7)) +
#   scale_fill_manual(values = acc_colors) +
#   scale_color_manual(values = acc_colors) +
#   coord_cartesian(clip = "off") +
#   theme_cowplot() +
#   theme(legend.position = "none",
#         axis.text = element_text(color = "grey40"),
#         axis.ticks = element_line(color = "grey40"),
#         axis.line = element_line(color = "grey40"),
#         axis.title = element_text(color = "grey40"))
# 
# p_lmd <- FetchData(object = dat, 
#                    vars = "lmd",
#                    layer = "data") |>
#   rownames_to_column("cell") |> 
#   pivot_longer(-cell, names_to = "feature", values_to = "data") |> 
#   full_join(meta) |> 
#   mutate(rug_position = ifelse(acc_temp == "18°C",  -0.07, -0.17)) |> 
#   filter(cell_type == "mesoderm prim.") |> 
#   ggplot() +
#   geom_label(aes(x = 2.6, y = 1.2),
#              label = "lmd",
#              fill = "#D39C9E",
#              color = "grey20",
#              size = 5) +
#   geom_density(aes(x = data,
#                    fill = acc_temp,
#                    color = acc_temp),
#                alpha = 0.5) +
#   geom_point(aes(y = rug_position, 
#                  x = data,
#                  color = acc_temp),
#              shape = "|", 
#              size = 2) +
#   scale_x_continuous(name = "",
#                      expand = c(0, 0.05)) +
#   scale_y_continuous(name = "",
#                      expand = c(0, 0.03),
#                      limits = c(-0.2, 1.4)) +
#   scale_fill_manual(values = acc_colors) +
#   scale_color_manual(values = acc_colors) +
#   coord_cartesian(clip = "off") +
#   theme_cowplot() +
#   theme(legend.position = "none",
#         axis.text = element_text(color = "grey40"),
#         axis.ticks = element_line(color = "grey40"),
#         axis.line = element_line(color = "grey40"))


p_ab <- FetchData(object = dat, 
                     vars = "ab",
                     layer = "data") |>
  rownames_to_column("cell") |> 
  pivot_longer(-cell, names_to = "feature", values_to = "data") |> 
  full_join(meta) |> 
  mutate(rug_position = ifelse(acc_temp == "18°C", -0.0325, -0.0825)) |> 
  filter(cell_type == "ventral nerve cord prim.") |> 
  ggplot() +
  geom_label(aes(x = 2.6, y = 0.6),
             label = "ab",
             fill = "#57A4B2",
             color = "grey20",
             size = 5) +
  geom_density(aes(x = data,
                   fill = acc_temp,
                   color = acc_temp),
               alpha = 0.5) +
  geom_point(aes(y = rug_position, 
                 x = data,
                 color = acc_temp),
             shape = "|", 
             size = 2) +
  scale_x_continuous(name = "",
                     expand = c(0, 0.05)) +
  scale_y_continuous(name = "",
                     expand = c(0, 0.015),
                     limits = c(-0.10, 0.7)) +
  scale_fill_manual(values = acc_colors) +
  scale_color_manual(values = acc_colors) +
  coord_cartesian(clip = "off") +
  theme_cowplot() +
  theme(legend.position = "none",
        axis.text = element_text(color = "grey40"),
        axis.ticks = element_line(color = "grey40"),
        axis.line = element_line(color = "grey40"))


# Combine rug plots
p_rugs <- cowplot::plot_grid(
  p_blimp,
  p_nub,
  p_ab,
  legend_density,
  rel_widths = c(1, 1, 1, 0.4),
  nrow = 1
)

# Target site accessibility 
dars_tf <- eReg |> 
  distinct(TF, Region_signature_name, Region, Gene) |> 
  mutate(Region = str_remove_all(Region, "chr")) |> 
  mutate(Region = str_replace(Region, ":", "-")) |> 
  right_join(dars, by = c("Region" = "region"),
             relationship = "many-to-many")

p_blimp_region <- dars_tf |>
  filter(Region_signature_name == "Blimp-1_+_(166r)" & 
           cell_type == "foregut & hindgut prim.") |>
  ggplot() + 
  geom_histogram(aes(x = avg_log2FC_18_25),
                 color = "grey20",
                 fill = "grey80",
                 bins = 30) +
  geom_vline(xintercept = 0,
             color = "grey40",
             linewidth = 1.2,
             linetype = 2) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(limits = c(-2.01, 2.01)) +
  labs(x = "",
       y = "Number of peaks") +
  theme_cowplot() +
  theme(axis.text = element_text(color = "grey40"),
        axis.ticks = element_line(color = "grey40"),
        axis.line = element_line(color = "grey40"),
        axis.title = element_text(color = "grey40"))

p_nub_region <- dars_tf |>
  filter(Region_signature_name == "nub_extended_+_(305r)" & 
           cell_type == "peripheral nervous system prim.") |>
  ggplot() + 
  geom_histogram(aes(x = avg_log2FC_18_25),
                 color = "grey20",
                 fill = "grey80",
                 bins = 30) +
  geom_vline(xintercept = 0,
             color = "grey40",
             linewidth = 1.2,
             linetype = 2) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(limits = c(-1.5, 1.5)) +
  labs(x = "",
       y = "") +
  theme_cowplot() +
  theme(axis.text = element_text(color = "grey40"),
        axis.ticks = element_line(color = "grey40"),
        axis.line = element_line(color = "grey40"))


# p_hth_region <- dars_tf |>
#   filter(TF == "hth" & 
#            cell_type == "mesoderm prim.") |>
#   ggplot() + 
#   geom_histogram(aes(x = avg_log2FC_18_25,fill=sig),
#                  color = "grey20",
#                  #fill = "grey80",
#                  bins = 30) +
#   geom_vline(xintercept = 0,
#              color = "grey40",
#              linewidth = 1.2,
#              linetype = 2) +
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_x_continuous(limits = c(-1, 1),
#                      breaks = c(-1, 0, 1)) +
#   scale_fill_manual(values = c("grey80", "firebrick")) +
#   labs(x = "log2(fold-change)",
#        y = "") +
#   theme_cowplot() +
#   theme(axis.text = element_text(color = "grey40"),
#         axis.ticks = element_line(color = "grey40"),
#         axis.line = element_line(color = "grey40"),
#         axis.title = element_text(color = "grey40"),
#         legend.position = "none")
# 
# p_lmd_region <- dars_tf |>
#   filter(TF == "lmd" & 
#            cell_type == "mesoderm prim.") |>
#   ggplot() + 
#   geom_histogram(aes(x = avg_log2FC_18_25, 
#                      fill = sig),
#                  color = "grey20",
#                  bins = 30) +
#   geom_vline(xintercept = 0,
#              color = "grey40",
#              linewidth = 1.2,
#              linetype = 2) +
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_x_continuous(limits = c(-1, 1),
#                      breaks = c(-1, 0, 1)) +
#   scale_fill_manual(values = c("grey80", "firebrick")) +
#   labs(x = "",
#        y = "") +
#   theme_cowplot() +
#   theme(axis.text = element_text(color = "grey40"),
#         axis.ticks = element_line(color = "grey40"),
#         axis.line = element_line(color = "grey40"),
#         legend.position = "none")

p_ab_region <- dars_tf |>
  filter(Region_signature_name == "ab_+_(694r)" & 
           cell_type == "ventral nerve cord prim.") |>
  ggplot() + 
  geom_histogram(aes(x = avg_log2FC_18_25),
                 color = "grey20",
                 fill = "grey80",
                 bins = 30) +
  geom_vline(xintercept = 0,
             color = "grey40",
             linewidth = 1.2,
             linetype = 2) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(limits = c(-2, 2)) +
  scale_fill_manual(values = c("grey80", "firebrick")) +
  labs(x = "",
       y = "") +
  theme_cowplot() +
  theme(axis.text = element_text(color = "grey40"),
        axis.ticks = element_line(color = "grey40"),
        axis.line = element_line(color = "grey40"))


p_reg <- cowplot::plot_grid(
  p_blimp_region,
  p_nub_region,
  p_ab_region,
  NULL,
  rel_widths = c(1, 1, 1, 0.4),
  nrow = 1
)

# Target gene expression panel
blimp_targets <- eReg |> 
  filter(Region_signature_name == "Blimp-1_+_(166r)") |> 
  distinct(Gene) |> 
  pull()

p_blimp_targets <- degs |> 
  filter(cell_type == "foregut & hindgut prim.") |> 
  mutate(target = ifelse(gene %in% blimp_targets, TRUE, FALSE)) |> 
  arrange(target, sig) |> 
  ggplot() +
  geom_abline(intercept = 0, slope = 0) +
  geom_point(aes(x = pct.25,
                 y = avg_log2FC_18_25,
                 color = target,
                 size = sig,
                 #alpha = sig,
                 shape = sig)) + 
  scale_color_manual(values = c("grey90", "firebrick"),
                     labels = c("other", "Blimp-1 target"),
                     name = element_blank()) +
  scale_size_manual(values = c(1.5, 3),
                    labels = c("non-DE gene","DE Gene"),
                    name = element_blank()) +
  scale_shape_manual(values = c(19, 17),
                     labels = c("non-DE gene","DE Gene"),
                     name = element_blank()) +
  scale_y_continuous(name = "Log2 Fold Change\n 18°C vs. 25°C",
                     limits = c(-0.9, 0.9)) +
  scale_x_continuous(name = "",
                     labels = scales::label_percent()) +
  theme_minimal() +
  theme(legend.position = "none")

nub_targets <- eReg |> 
  filter(Region_signature_name == "nub_extended_+_(305r)") |> 
  distinct(Gene) |> 
  pull()

p_nub_targets <- degs |> 
  filter(cell_type == "peripheral nervous system prim.") |> 
  mutate(target = ifelse(gene %in% nub_targets, TRUE, FALSE)) |> 
  arrange(target, sig) |> 
  ggplot() +
  geom_abline(intercept = 0, slope = 0) +
  geom_point(aes(x = pct.25,
                 y = avg_log2FC_18_25,
                 color = target,
                 size = sig,
                 #alpha = sig,
                 shape = sig)) + 
  scale_color_manual(values = c("grey90", "firebrick"),
                     labels = c("other", "nub target"),
                     name = element_blank()) +
  scale_size_manual(values = c(1.5, 3),
                    labels = c("non-DE gene","DE Gene"),
                    name = element_blank()) +
  scale_shape_manual(values = c(19, 17),
                     labels = c("non-DE gene","DE Gene"),
                     name = element_blank()) +
  scale_y_continuous(name = "",
                     limits = c(-0.9, 0.9)) +
  scale_x_continuous(name = "",
                     labels = scales::label_percent()) +
  theme_minimal() +
  theme(legend.position = "none")

# hth_targets <- eReg |> 
#   filter(TF == "hth") |> 
#   distinct(Gene) |> 
#   pull()
# 
# p_hth_targets <- degs |> 
#   filter(cell_type == "mesoderm prim.") |> 
#   mutate(target = ifelse(gene %in% hth_targets, TRUE, FALSE)) |> 
#   arrange(target, sig) |> 
#   ggplot() +
#   geom_abline(intercept = 0, slope = 0) +
#   geom_point(aes(x = pct.25,
#                  y = avg_log2FC_18_25,
#                  color = target,
#                  size = sig,
#                  #alpha = sig,
#                  shape = sig)) + 
#   scale_color_manual(values = c("grey90", "firebrick"),
#                      labels = c("other", "hth target"),
#                      name = element_blank()) +
#   scale_size_manual(values = c(1.5, 3),
#                     labels = c("non-DE gene","DE Gene"),
#                     name = element_blank()) +
#   scale_shape_manual(values = c(19, 17),
#                      labels = c("non-DE gene","DE Gene"),
#                      name = element_blank()) +
#   scale_y_continuous(name = "",
#                      limits = c(-0.9, 0.9)) +
#   scale_x_continuous(name = "Percent nuclei expressed (25°C)",
#                      labels = scales::label_percent()) +
#   theme_minimal() +
#   theme(legend.position = "none")

# lmd_targets <- eReg |> 
#   filter(TF == "lmd") |> 
#   distinct(Gene) |> 
#   pull()
# 
# p_lmd_targets <- degs |> 
#   filter(cell_type == "mesoderm prim.") |> 
#   mutate(target = ifelse(gene %in% ab_targets, TRUE, FALSE)) |> 
#   arrange(target, sig) |> 
#   ggplot() +
#   geom_abline(intercept = 0, slope = 0) +
#   geom_point(aes(x = pct.25,
#                  y = avg_log2FC_18_25,
#                  color = target,
#                  size = sig,
#                  #alpha = sig,
#                  shape = sig)) + 
#   scale_color_manual(values = c("grey90", "firebrick"),
#                      labels = c("other", "lmd target"),
#                      name = element_blank()) +
#   scale_size_manual(values = c(1.5, 3),
#                     labels = c("non-DE gene","DE Gene"),
#                     name = element_blank()) +
#   scale_shape_manual(values = c(19, 17),
#                      labels = c("non-DE gene","DE Gene"),
#                      name = element_blank()) +
#   scale_y_continuous(name = "",
#                      limits = c(-0.9, 0.9)) +
#   scale_x_continuous(name = "",
#                      labels = scales::label_percent()) +
#   theme_minimal() +
#   theme(legend.position = "none")

ab_targets <- eReg |> 
  filter(Region_signature_name == "ab_+_(694r)") |> 
  distinct(Gene) |> 
  pull()

p_ab_targets_leg <- degs |> 
  filter(cell_type == "ventral nerve cord prim.") |> 
  mutate(target = ifelse(gene %in% ab_targets, TRUE, FALSE)) |> 
  arrange(target, sig) |> 
  ggplot() +
  geom_abline(intercept = 0, slope = 0) +
  geom_point(aes(x = pct.25,
                 y = avg_log2FC_18_25,
                 color = target,
                 size = sig,
                 #alpha = sig,
                 shape = sig)) + 
  scale_color_manual(values = c("grey90", "firebrick"),
                     labels = c("other", "TF target"),
                     name = element_blank()) +
  scale_size_manual(values = c(1.5, 3),
                    labels = c("non-sig gene", "DE Gene"),
                    name = element_blank()) +
  scale_shape_manual(values = c(19, 17),
                     labels = c("non-sig gene", "DE Gene"),
                     name = element_blank()) +
  scale_y_continuous(name = "",
                     limits = c(-0.9, 0.9)) +
  scale_x_continuous(name = "",
                     labels = scales::label_percent()) +
  theme_minimal()

p_ab_targets <- p_ab_targets_leg + 
  theme(legend.position = "none")

target_exp_leg <- ggdraw() + 
  get_legend(p_ab_targets_leg)

p_ma <- cowplot::plot_grid(
  p_blimp_targets,
  p_nub_targets,
  p_ab_targets,
  target_exp_leg,
  rel_widths = c(1, 1, 1, 0.4),
  nrow = 1
)


a_title <- ggplot() + 
  annotate("label", 
           x = 0, 
           y = 0, 
           color = "grey20",
           fill = "grey90",
           label = "Cell-type-specific gene regulatory networks") +
  theme_void()
b_title <- ggplot() + 
  annotate("label", 
           x = 0, 
           y = 0, 
           color = "grey20",
           fill = "grey90",
           label = "Transcription factor expression") +
  theme_void()
c_title <- ggplot() + 
  annotate("label", 
           x = 0, 
           y = 0, 
           color = "grey20",
           fill = "grey90",
           label = "Target region accessibility") +
  theme_void()
d_title <- ggplot() + 
  annotate("label", 
           x = 0, 
           y = 0, 
           color = "grey20",
           fill = "grey90",
           label = "Target gene expression") +
  theme_void()

# Put together the whole figure
f4 <- cowplot::plot_grid(
  a_title,
  p_ereg,
  b_title,
  p_rugs,
  c_title,
  p_reg,
  d_title,
  p_ma,
  labels = c("A", "", "B", "", "C", "", "D", ""),
  rel_heights = c(0.1, 2.25, 0.15, 0.7, 0.15, 0.7, 0.15, 0.7),
  greedy = FALSE,
  nrow = 8
)

# Save
ggsave(
  here::here(fig_dir, "fig_4_test.png"),
  plot = f4,
  bg = "white",
  height = 40,
  width = 35,
  units = "cm"
)
ggsave(
  here::here(fig_dir, "fig_4_test.pdf"),
  plot = f4,
  bg = "white",
  height = 40,
  width = 35,
  units = "cm"
)




########################### new idea
# blimp_region_targets <- dars_tf |>
#   filter(Region_signature_name == "Blimp-1_+_(166r)" & 
#            cell_type == "foregut & hindgut prim.") |> 
#   dplyr::rename("gene" = "Gene")
# 
# blimp_gene_targets <- degs |> 
#   filter(cell_type == "foregut & hindgut prim.") |> 
#   mutate(target = ifelse(gene %in% blimp_targets, TRUE, FALSE)) |> 
#   filter(target)
# 
# blimp <- full_join(
#   blimp_region_targets, 
#   blimp_gene_targets,
#   by = "gene",
#   suffix = c(".region", ".gene")) |> 
#   add_column(avg_log2FC_18_25.TF = 1.41185430) |>
#   pivot_longer(cols = contains("avg_log2FC"),
#                names_to = "group", 
#                values_to = "L2FC") |> 
#   mutate(group = str_remove_all(group, "avg_log2FC_18_25.")) |> 
#   mutate(group = factor(group, levels = c("TF", "region", "gene")))

# ------------------------------------------------------------------------------
GRN <- full_join(
  dars_tf |> dplyr::rename("gene" = "Gene"), 
  degs,
  by = c("gene", "cell_type"),
  suffix = c(".region", ".gene")) |> 
  right_join(degs |> 
               filter(sig) |> 
               filter(gene %in% unique(eReg$TF)) |> 
               rename_with(.cols = p_val:sig, function(x){paste0(x, ".TF")}),
             by = c("TF" = "gene", "cell_type")) |> 
  pivot_longer(cols = contains("avg_log2FC"),
               names_to = "group", 
               values_to = "L2FC") |> 
  mutate(group = str_remove_all(group, "avg_log2FC_18_25.")) |> 
  mutate(group = factor(group, levels = c("TF", "region", "gene")))

GRN |>
  mutate(cell_type_eGRN = paste(Region_signature_name, cell_type, sep = "\n")) |> 
  filter(cell_type_eGRN %in% diffGRN) |> 
  # filter((Region_signature_name == "Blimp-1_+_(166r)" &
  #          cell_type == "foregut & hindgut prim.") |
  #          (Region_signature_name == "nub_extended_+_(305r)" &
  #             cell_type == "peripheral nervous system prim.") |
  #          (Region_signature_name == "ab_+_(694r)" &
  #             cell_type == "ventral nerve cord prim.")) |>
  filter(pct.18.gene >= 0.1 | pct.25.gene >= 0.1,
         pct.18.region >= 0.1 | pct.25.region >= 0.1) |> 
  # mutate(Region_signature_name = factor(Region_signature_name,
  #                                       levels = c(
  #                                         "Blimp-1_+_(166r)",
  #                                         "nub_extended_+_(305r)",
  #                                         "ab_+_(694r)"))) |>
  ggplot(
    aes(x = group,
        y = L2FC)) +
  geom_hline(yintercept = 0,
             color = "grey20") +
  #geom_violin(alpha = 0.2) +
  geom_line(aes(group = paste(TF, Region, gene)),
            size = 0.2,
            color = "grey50") +
  geom_point(aes(size = group),
             color = "grey50",
             fill = "grey80",
             shape = 21) +

  scale_size_manual(values = c(4, 2, 2)) +
  scale_y_continuous(limits = c(-3, 3),
                     name = "log2(fold-change)") +
  scale_x_discrete(name = element_blank(),
                   labels = c("Target gene exp", 
                              "Motif peak accessibility",
                              "TF exp"),
                   limits = rev) +
  cowplot::theme_minimal_grid() +
  theme(axis.line.x = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(size = 10)) +
  facet_wrap(~cell_type_eGRN,
             nrow = 2,
             labeller = label_wrap_gen(width = 20)) +
  coord_flip()



# Figure 5 ---------------------------------------------------------------------

## WGCNA AND SHIT

library(hdWGCNA)

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/13_dat_wgcna.rds"))
DMEs <- readRDS(here::here("output/wgcna/DMEs.rds"))


mod_tally <- GetModules(dat) |> 
  group_by(module) |> 
  tally()

GetModules(dat) |> 
  filter(module == "green")

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

p_wgcna <- left_join(x, mod_tally) |> 
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
  scale_size_continuous(range = c(5, 10),
                        name = "Number of\ngenes",
                        breaks = c(75, 100, 250, 500)) +
  labs(y = element_blank()) +
  cowplot::theme_minimal_vgrid() +
  theme(panel.grid.major.y = element_line(color = "grey95"),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) +
  guides(color = "none",
         fill = "none")


# Set
DefaultAssay(dat) <- "RNA"


p_rpl27 <- FetchData(object = dat, 
                   vars = "RpL27A",
                   layer = "data") |>
  rownames_to_column("cell") |> 
  pivot_longer(-cell, names_to = "feature", values_to = "data") |> 
  full_join(meta) |> 
  mutate(rug_position = ifelse(acc_temp == "18°C", -0.02, -0.055)) |> 
  ggplot() +
  geom_label(aes(x = 2.5, y = 0.7),
             label = "Rpl27A",
             fill = c("red") |> colorspace::desaturate(amount = 0.5),
             color = "white",
             size = 5) +
  geom_density(aes(x = data,
                   fill = acc_temp,
                   color = acc_temp),
               alpha = 0.5) +
  geom_point(aes(y = rug_position, 
                 x = data,
                 color = acc_temp),
             shape = "|", 
             size = 1.5) +
  scale_x_continuous(name = "",
                     expand = c(0, 0.05)) +
  scale_y_continuous(name = "Density",
                     expand = c(0, 0.005),
                     limits = c(-0.07, 0.8)) +
  scale_fill_manual(values = acc_colors) +
  scale_color_manual(values = acc_colors) +
  coord_cartesian(clip = "off") +
  theme_cowplot() +
  theme(legend.position = "none",
        axis.text = element_text(color = "grey40"),
        axis.ticks = element_line(color = "grey40"),
        axis.line = element_line(color = "grey40"))

p_rplp1 <- FetchData(object = dat, 
                     vars = "RpLP1",
                     layer = "data") |>
  rownames_to_column("cell") |> 
  pivot_longer(-cell, names_to = "feature", values_to = "data") |> 
  full_join(meta) |> 
  mutate(rug_position = ifelse(acc_temp == "18°C", -0.02, -0.055)) |> 
  ggplot() +
  geom_label(aes(x = 2.5, y = 0.7),
             label = "RpLP1",
             fill = c("red") |> colorspace::desaturate(amount = 0.5),
             color = "white",
             size = 5) +
  geom_density(aes(x = data,
                   fill = acc_temp,
                   color = acc_temp),
               alpha = 0.5) +
  geom_point(aes(y = rug_position, 
                 x = data,
                 color = acc_temp),
             shape = "|", 
             size = 1.5) +
  scale_x_continuous(name = "Log-normalized expression",
                     expand = c(0, 0.05)) +
  scale_y_continuous(name = "",
                     expand = c(0, 0.005),
                     limits = c(-0.07, 0.8)) +
  scale_fill_manual(values = acc_colors) +
  scale_color_manual(values = acc_colors) +
  coord_cartesian(clip = "off") +
  theme_cowplot() +
  theme(legend.position = "none",
        axis.text = element_text(color = "grey40"),
        axis.ticks = element_line(color = "grey40"),
        axis.line = element_line(color = "grey40"))

p_rpl15 <- FetchData(object = dat, 
                     vars = "RpL15",
                     layer = "data") |>
  rownames_to_column("cell") |> 
  pivot_longer(-cell, names_to = "feature", values_to = "data") |> 
  full_join(meta) |> 
  mutate(rug_position = ifelse(acc_temp == "18°C", -0.02, -0.055)) |> 
  ggplot() +
  geom_label(aes(x = 2.5, y = 0.7),
             label = "RpL15",
             fill = c("red") |> colorspace::desaturate(amount = 0.5),
             color = "white",
             size = 5) +
  geom_density(aes(x = data,
                   fill = acc_temp,
                   color = acc_temp),
               alpha = 0.5) +
  geom_point(aes(y = rug_position, 
                 x = data,
                 color = acc_temp),
             shape = "|", 
             size = 1.5) +
  scale_x_continuous(name = "",
                     expand = c(0, 0.05)) +
  scale_y_continuous(name = "",
                     expand = c(0, 0.005),
                     limits = c(-0.07, 0.8)) +
  scale_fill_manual(values = acc_colors) +
  scale_color_manual(values = acc_colors) +
  coord_cartesian(clip = "off") +
  theme_cowplot() +
  theme(legend.position = "none",
        axis.text = element_text(color = "grey40"),
        axis.ticks = element_line(color = "grey40"),
        axis.line = element_line(color = "grey40"))


p_red <- cowplot::plot_grid(
  p_rpl27,
  p_rplp1,
  p_rpl15,
  nrow = 1
)

p_hmgz <- FetchData(object = dat, 
                     vars = "HmgZ",
                     layer = "data") |>
  rownames_to_column("cell") |> 
  pivot_longer(-cell, names_to = "feature", values_to = "data") |> 
  full_join(meta) |> 
  mutate(rug_position = ifelse(acc_temp == "18°C", -0.02, -0.055)) |> 
  ggplot() +
  geom_label(aes(x = 2.5, y = 0.7),
             label = "HmgZ",
             fill = c("brown") |> colorspace::desaturate(amount = 0.5),
             color = "white",
             size = 5) +
  geom_density(aes(x = data,
                   fill = acc_temp,
                   color = acc_temp),
               alpha = 0.5) +
  geom_point(aes(y = rug_position, 
                 x = data,
                 color = acc_temp),
             shape = "|", 
             size = 1.5) +
  scale_x_continuous(name = "",
                     expand = c(0, 0.05)) +
  scale_y_continuous(name = "Density",
                     expand = c(0, 0.005),
                     limits = c(-0.07, 0.8)) +
  scale_fill_manual(values = acc_colors) +
  scale_color_manual(values = acc_colors) +
  coord_cartesian(clip = "off") +
  theme_cowplot() +
  theme(legend.position = "none",
        axis.text = element_text(color = "grey40"),
        axis.ticks = element_line(color = "grey40"),
        axis.line = element_line(color = "grey40"))

p_cox <- FetchData(object = dat, 
                    vars = "COX5A",
                    layer = "data") |>
  rownames_to_column("cell") |> 
  pivot_longer(-cell, names_to = "feature", values_to = "data") |> 
  full_join(meta) |> 
  mutate(rug_position = ifelse(acc_temp == "18°C", -0.02, -0.055)) |> 
  ggplot() +
  geom_label(aes(x = 2.5, y = 0.7),
             label = "COX5A",
             fill = c("brown") |> colorspace::desaturate(amount = 0.5),
             color = "white",
             size = 5) +
  geom_density(aes(x = data,
                   fill = acc_temp,
                   color = acc_temp),
               alpha = 0.5) +
  geom_point(aes(y = rug_position, 
                 x = data,
                 color = acc_temp),
             shape = "|", 
             size = 1.5) +
  scale_x_continuous(name = "Log-normalized expression",
                     expand = c(0, 0.05)) +
  scale_y_continuous(name = "",
                     expand = c(0, 0.005),
                     limits = c(-0.07, 0.8)) +
  scale_fill_manual(values = acc_colors) +
  scale_color_manual(values = acc_colors) +
  coord_cartesian(clip = "off") +
  theme_cowplot() +
  theme(legend.position = "none",
        axis.text = element_text(color = "grey40"),
        axis.ticks = element_line(color = "grey40"),
        axis.line = element_line(color = "grey40"))

p_atpsyn <- FetchData(object = dat, 
                   vars = "ATPsynE",
                   layer = "data") |>
  rownames_to_column("cell") |> 
  pivot_longer(-cell, names_to = "feature", values_to = "data") |> 
  full_join(meta) |> 
  mutate(rug_position = ifelse(acc_temp == "18°C", -0.02, -0.055)) |> 
  ggplot() +
  geom_label(aes(x = 2.5, y = 0.7),
             label = "ATPsynE",
             fill = c("brown") |> colorspace::desaturate(amount = 0.5),
             color = "white",
             size = 5) +
  geom_density(aes(x = data,
                   fill = acc_temp,
                   color = acc_temp),
               alpha = 0.5) +
  geom_point(aes(y = rug_position, 
                 x = data,
                 color = acc_temp),
             shape = "|", 
             size = 1.5) +
  scale_x_continuous(name = "",
                     expand = c(0, 0.05)) +
  scale_y_continuous(name = "",
                     expand = c(0, 0.005),
                     limits = c(-0.07, 0.8)) +
  scale_fill_manual(values = acc_colors) +
  scale_color_manual(values = acc_colors) +
  coord_cartesian(clip = "off") +
  theme_cowplot() +
  theme(legend.position = "none",
        axis.text = element_text(color = "grey40"),
        axis.ticks = element_line(color = "grey40"),
        axis.line = element_line(color = "grey40"))

p_brown <- cowplot::plot_grid(
  p_hmgz,
  p_cox,
  p_atpsyn,
  nrow = 1
)

p_how <- FetchData(object = dat, 
                    vars = "how",
                    layer = "data") |>
  rownames_to_column("cell") |> 
  pivot_longer(-cell, names_to = "feature", values_to = "data") |> 
  full_join(meta) |> 
  mutate(rug_position = ifelse(acc_temp == "18°C", -0.02, -0.055)) |> 
  ggplot() +
  geom_label(aes(x = 2.5, y = 0.7),
             label = "how",
             fill = c("green") |> colorspace::desaturate(amount = 0.5),
             color = "grey20",
             size = 5) +
  geom_density(aes(x = data,
                   fill = acc_temp,
                   color = acc_temp),
               alpha = 0.5) +
  geom_point(aes(y = rug_position, 
                 x = data,
                 color = acc_temp),
             shape = "|", 
             size = 1.5) +
  scale_x_continuous(name = "",
                     expand = c(0, 0.05)) +
  scale_y_continuous(name = "Density",
                     expand = c(0, 0.005),
                     limits = c(-0.07, 0.8)) +
  scale_fill_manual(values = acc_colors) +
  scale_color_manual(values = acc_colors) +
  coord_cartesian(clip = "off") +
  theme_cowplot() +
  theme(legend.position = "none",
        axis.text = element_text(color = "grey40"),
        axis.ticks = element_line(color = "grey40"),
        axis.line = element_line(color = "grey40"))

p_jeb <- FetchData(object = dat, 
                   vars = "jeb",
                   layer = "data") |>
  rownames_to_column("cell") |> 
  pivot_longer(-cell, names_to = "feature", values_to = "data") |> 
  full_join(meta) |> 
  mutate(rug_position = ifelse(acc_temp == "18°C", -0.02, -0.055)) |> 
  ggplot() +
  geom_label(aes(x = 2.5, y = 0.7),
             label = "jeb",
             fill = c("green") |> colorspace::desaturate(amount = 0.5),
             color = "grey20",
             size = 5) +
  geom_density(aes(x = data,
                   fill = acc_temp,
                   color = acc_temp),
               alpha = 0.5) +
  geom_point(aes(y = rug_position, 
                 x = data,
                 color = acc_temp),
             shape = "|", 
             size = 1.5) +
  scale_x_continuous(name = "Log-normalized expression",
                     expand = c(0, 0.05)) +
  scale_y_continuous(name = "",
                     expand = c(0, 0.005),
                     limits = c(-0.07, 0.8)) +
  scale_fill_manual(values = acc_colors) +
  scale_color_manual(values = acc_colors) +
  coord_cartesian(clip = "off") +
  theme_cowplot() +
  theme(legend.position = "none",
        axis.text = element_text(color = "grey40"),
        axis.ticks = element_line(color = "grey40"),
        axis.line = element_line(color = "grey40"))

p_osp <- FetchData(object = dat, 
                      vars = "osp",
                      layer = "data") |>
  rownames_to_column("cell") |> 
  pivot_longer(-cell, names_to = "feature", values_to = "data") |> 
  full_join(meta) |> 
  mutate(rug_position = ifelse(acc_temp == "18°C", -0.02, -0.055)) |> 
  ggplot() +
  geom_label(aes(x = 2.5, y = 0.7),
             label = "osp",
             fill = c("green") |> colorspace::desaturate(amount = 0.5),
             color = "grey20",
             size = 5) +
  geom_density(aes(x = data,
                   fill = acc_temp,
                   color = acc_temp),
               alpha = 0.5) +
  geom_point(aes(y = rug_position, 
                 x = data,
                 color = acc_temp),
             shape = "|", 
             size = 1.5) +
  scale_x_continuous(name = "",
                     expand = c(0, 0.05)) +
  scale_y_continuous(name = "",
                     expand = c(0, 0.005),
                     limits = c(-0.07, 0.8)) +
  scale_fill_manual(values = acc_colors) +
  scale_color_manual(values = acc_colors) +
  coord_cartesian(clip = "off") +
  theme_cowplot() +
  theme(legend.position = "none",
        axis.text = element_text(color = "grey40"),
        axis.ticks = element_line(color = "grey40"),
        axis.line = element_line(color = "grey40"))

p_green <- cowplot::plot_grid(
  p_how,
  p_jeb,
  p_osp,
  nrow = 1
)


p_exp_mods <- cowplot::plot_grid(
  p_red,
  p_brown,
  p_green,
  labels = c("B", "C", "D"),
  nrow = 3
)

# ggsave(here::here(fig_dir, "fig_5_test.png"),
#        plot = p_exp_mods,
#        height = 30,
#        width = 30,
#        bg = "white",
#        units = "cm")

# Module genes

color_cell_type <- tibble(
  cell_type = c(
    "germ cell", 
    "peripheral nervous system prim.",
    "ectoderm prim.",
    "mesoderm prim.",
    "endoderm prim.",
    "foregut & hindgut prim.",
    "ventral nerve cord prim.",
    "tracheal prim.",
    "amnioserosa"
  ), 
  colors = c(
    "grey90", 
    "#ADD9F4",
    "#57A4B2",
    "#D39C9E",
    "#FEF29A",
    "#F9DCEE",
    "#819FC5",
    "#A7BF9B",
    "#bfa3a4"
  )
)


degs_df <- degs |>  
  filter(sig == TRUE) |> 
  group_by(cell_type) |> 
  add_tally() |> 
  full_join(color_cell_type) 

dars_df <- dars |>  
  filter(sig == TRUE) |> 
  group_by(cell_type) |> 
  add_tally() |> 
  full_join(color_cell_type) 

cell_type_order <- degs_df |>
  arrange(n) |> 
  distinct(cell_type) |> 
  pull()
  
p_degs <- degs_df |> 
  ggplot(aes(x = avg_log2FC_18_25,
             y = factor(cell_type, levels = cell_type_order),
             fill = colors,
             color = colors)) +
  ggbeeswarm::geom_beeswarm(
    shape = 21,
    color = "grey60",
    size = 0.9,
    cex = 0.9,
    side = 1,
    alpha = 0.8) +
  geom_vline(xintercept = 0,
             color = "grey30",
             linetype = 5) +
  geom_label(data = degs_df |> distinct(cell_type, n, colors),
             aes(x = -1.2,
                 y = cell_type,
                 label = n),
             color = "grey20",
             vjust = -0.5,
             alpha = 0.4) +
  scale_y_discrete(name = element_blank(), 
                   expand = c(0, 0.05)) +
  scale_x_continuous(name = "Log2 fold-change",
                     limits = c(-1.25, 1.25),
                     breaks = c(-1, -0.5, 0, 0.5, 1)) +
  scale_fill_identity() +
  scale_color_identity() +
  ggridges::geom_density_ridges(alpha = 0.5) +
  theme_minimal_hgrid() +
  theme(axis.text.y = element_blank())

p_dars <- dars_df |> 
  ggplot(aes(x = avg_log2FC_18_25,
             y = factor(cell_type, levels = cell_type_order),
             fill = colors,
             color = colors)) +
  ggbeeswarm::geom_beeswarm(
    shape = 21,
    color = "grey60",
    size = 1.5,
    cex = 1.05,
    side = 1,
    alpha = 0.8) +
  geom_vline(xintercept = 0,
             color = "grey30",
             linetype = 5) +
  geom_label(data = dars_df |> distinct(cell_type, n, colors) |> 
               mutate(n = ifelse(is.na(n), 0, n)),
             aes(x = -3,
                 y = cell_type,
                 label = n),
             color = "grey20",
             vjust = -0.5,
             alpha = 0.4) +
  scale_y_discrete(name = element_blank(),
                   expand = c(0, 0.05)) +
  scale_x_continuous(name = "Log2 fold-change",
                     limits = c(-3.25, 3.25),
                     breaks = c(-3:3)) +
  scale_fill_identity() +
  scale_color_identity() +
  ggridges::geom_density_ridges(alpha = 0.5) +
  theme_minimal_hgrid()

p_degs_h <- cowplot::plot_grid(
  NULL,
  p_degs,
  rel_heights = c(0.03, 1),
  ncol = 1
)

p_ridge <- cowplot::plot_grid(
  p_dars, 
  p_degs_h, 
  nrow = 1,
  labels = c("E", "F"),
  hjust = c(-18, -0.5),
  rel_widths = c(1.25, 0.75)
)

p_wgcna_all <- cowplot::plot_grid(
  p_wgcna, 
  p_exp_mods, 
  labels = c("A", NULL),
  nrow = 1,
  rel_widths = c(1, 1)
)

f5 <- cowplot::plot_grid(
  p_wgcna_all, 
  p_ridge, 
  nrow = 2
)


dars_motif_df <- readRDS(here::here("output/tf_motif/dars_df_motif.rds"))

zelda_motif <- "MA1462.1"
hsf_motif <- "MA1458.1"

rp_motif <- "UN0756.1"

dat@assays$peaks@motifs@positions[rp_motif]

zelda_motif_ranges <- dat@assays$peaks@motifs@positions[zelda_motif]
hsf_motif_ranges <- dat@assays$peaks@motifs@positions[hsf_motif]

# Set of GRanges object of the pseudobulk dars
dars_granges <- dars |> 
  tibble() |> 
  distinct(region) |> 
  filter(str_detect(region, "mitochondrion", negate = TRUE)) |> 
  separate(region, into = c("seqnames", "start", "end"), sep = "-") |> 
  mutate(strand = "*") |> 
  GenomicRanges::makeGRangesFromDataFrame()

# Find the dars that include the motifs
zelda_target <- subsetByOverlaps(dars_granges, zelda_motif_ranges) |> 
  as_tibble() |> 
  mutate(peaks = paste(seqnames, start, end, sep = "-")) 
hsf_target <- subsetByOverlaps(dars_granges, hsf_motif_ranges) |> 
  as_tibble() |> 
  mutate(peaks = paste(seqnames, start, end, sep = "-")) 
  

readRDS(here::here("output/dars/dars_bulk.rds")) |> 
  as.data.frame() |> 
  rownames_to_column("region") |> 
  mutate(zelda = ifelse(region %in% zelda_target$peaks, TRUE, FALSE)) |> 
  filter(zelda) |> 
  filter(sig) |> 
  ggplot() +
  geom_histogram(aes(x = avg_log2FC)) 


readRDS(here::here("output/dars/dars_bulk.rds")) |> 
  as.data.frame() |> 
  rownames_to_column("region") |> 
  mutate(zelda = ifelse(region %in% hsf_target$peaks, TRUE, FALSE)) |> 
  filter(zelda) |> 
  #filter(sig) |> 
  ggplot() +
  geom_histogram(aes(x = avg_log2FC)) 

dars |> 
  mutate(zelda = ifelse(region %in% hsf_target$peaks, TRUE, FALSE)) |> 
  filter(zelda) |> 
  filter(sig) |> 
  ggplot() +
  geom_histogram(aes(x = avg_log2FC_18_25)) +
  facet_wrap(~cell_type)

# Save
ggsave(here::here(fig_dir, "fig_5_test.png"),
       plot = f5,
       height = 30,
       width = 30,
       bg = "white",
       units = "cm")
ggsave(here::here(fig_dir, "fig_5_test.pdf"),
       plot = f5,
       height = 30,
       width = 30,
       bg = "white",
       units = "cm")

# ------------------------------------------------------------------------------
# Supplementary figures --------------------------------------------------------
# ------------------------------------------------------------------------------

# Figure S1 --------------------------------------------------------------------

# Define count cutoffs
low_ATAC <- 800
low_RNA <- 200

# Knee-plot for ATAC 
p_knee_atac <- tibble(nCount_ATAC = sort(dat_raw$nCount_ATAC,
                                         decreasing = TRUE)) |>
  rownames_to_column("rank") |>
  mutate(rank = as.numeric(rank)) |>
  head(100000) |>
  ggplot() +
  geom_point(aes(x = rank,
                 y = nCount_ATAC),
             color = "grey50",
             shape = 21) +
  geom_hline(aes(yintercept = low_ATAC),
             linetype = 2,
             color = "grey20",
             alpha = 0.8,
             linewidth = 1.1) +
  scale_y_continuous(trans = "log10", 
                     breaks = c(100, 1000, 10000, 100000),
                     labels = c("100", "1k", "10k", "100k"),
                     name = "ATAC counts") +
  scale_x_continuous(trans = "log10",
                     breaks = c(1, 10, 100, 1000, 10000, 100000),
                     labels = c("1", "10", "100", "1k", "10k", "100k"),
                     name = "Barcodes in rank order") +
  theme_minimal_grid()

# Knee-plot for RNA 
p_knee_rna <- tibble(nCount_RNA = sort(dat_raw$nCount_RNA,
                                       decreasing = TRUE)) |>
  rownames_to_column("rank") |>
  mutate(rank = as.numeric(rank)) |>
  head(100000) |>
  ggplot() +
  geom_point(aes(x = rank,
                 y = nCount_RNA),
             color = "grey50",
             shape = 21) +
  geom_hline(aes(yintercept = low_RNA),
             linetype = 2,
             color = "grey20",
             alpha = 0.8,
             linewidth = 1.1) +
  scale_y_continuous(trans = "log10", 
                     breaks = c(100, 1000, 10000, 100000),
                     labels = c("100", "1k", "10k", "100k"),
                     name = "RNA counts") +
  scale_x_continuous(trans = "log10",
                     breaks = c(1, 10, 100, 1000, 10000, 100000),
                     labels = c("1", "10", "100", "1k", "10k", "100k"),
                     name = "Barcodes in rank order") +
  theme_minimal_grid()


# Knee plots together
p_knees <- cowplot::plot_grid(
  p_knee_atac, 
  p_knee_rna,
  labels = c("A", "B"),
  nrow = 1
)

# Scatter plot
p_scatter <- dat_raw@meta.data |>
  filter(nCount_RNA >= 1 &
           nCount_ATAC >= 1) |>
  rownames_to_column("cells") |>
  mutate(cells_10x = ifelse(cells %in% Cells(dat) &
                              nCount_RNA >= low_RNA &
                              nCount_ATAC >= low_ATAC, 
                            "Nuclei", 
                            "Non-nuclei barcodes")) |>
  mutate(cells_10x = factor(cells_10x,
                            levels = c("Non-nuclei barcodes",
                                       "Nuclei"))) |>
  arrange(cells_10x) |> 
  ggplot(aes(x = nCount_RNA,
             y = nCount_ATAC)) +
  geom_point(aes(color = cells_10x),
             shape = 21) +
  geom_segment(aes(x = low_RNA, xend = low_RNA,
                   y = low_ATAC, yend = 150000),
               linetype = 2,
               color = "grey50",
               alpha = 0.8,
               linewidth = 1.1) +
  geom_segment(aes(y = low_ATAC, yend = low_ATAC,
                   x = low_RNA, xend = 50000),
               linetype = 2,
               color = "grey50",
               alpha = 0.8,
               linewidth = 1.1) +
  scale_x_continuous(trans = "log10", 
                     breaks = c(1, 10, 100, 1000, 10000, 100000),
                     labels = c("1", "10", "100", "1k", "10k", "100k"),
                     name = "RNA counts") +
  scale_y_continuous(trans = "log10", 
                     breaks = c(1, 10, 100, 1000, 10000, 100000),
                     labels = c("1", "10", "100", "1k", "10k", "100k"),
                     name = "ATAC counts") +
  scale_color_manual(values = c("grey80", "#91AB73"),
                     name = element_blank()) +
  theme_minimal_grid() +
  theme(legend.position = "bottom")

# Number of cells per sample
p_cells_sample <- dat@meta.data |>
  group_by(sample_name) |>
  mutate(sample_name = factor(sample_name, levels = c(c("25C_Rep2",
                                                        "25C_Rep1",
                                                        "18C_Rep2",
                                                        "18C_Rep1")))) |>
  tally() |>
  ggplot(aes(x = n/10^3,
             y = sample_name,
             label = scales::comma_format()(n))) +
  geom_col(color = "grey20",
           fill = "grey80",
           width = 0.8)+
  scale_x_continuous(expand = c(0, 0),
                     name = "Number of nuclei",
                     labels = 
                       scales::label_number(scale = 1, suffix = "k")) +
  geom_label(nudge_x = -0.5,
             label.size = 0.1) +
  scale_y_discrete(name = element_blank(),
                   labels = c("25°C Rep 2",
                              "25°C Rep 1",
                              "18°C Rep 2",
                              "18°C Rep 1")) +
  cowplot::theme_minimal_vgrid() 

# Total number of ATAC reads
p_reads_atac <- dat@meta.data |>
  group_by(sample_name) |>
  mutate(sample_name = factor(sample_name, levels = c(c("25C_Rep2",
                                                        "25C_Rep1",
                                                        "18C_Rep2",
                                                        "18C_Rep1")))) |>
  summarize(n = sum(nCount_ATAC)) |>
  ggplot(aes(x = n/10^6,
             y = sample_name,
             label = scales::comma_format()(n))) +
  geom_col(color = "grey20",
           fill = "grey80",
           width = 0.8)+
  scale_x_continuous(expand = c(0, 0),
                     name = "Number of ATAC reads",
                     labels = 
                       scales::label_number(scale = 1, suffix = "M")) +
  geom_label(nudge_x = -5,
             label.size = 0.1) +
  scale_y_discrete(name = element_blank(),
                   labels = c("25°C Rep 2",
                              "25°C Rep 1",
                              "18°C Rep 2",
                              "18°C Rep 1")) +
  cowplot::theme_minimal_vgrid() 

# Total number of ATAC reads
p_reads_rna <- dat@meta.data |>
  group_by(sample_name) |>
  mutate(sample_name = factor(sample_name, levels = c(c("25C_Rep2",
                                                        "25C_Rep1",
                                                        "18C_Rep2",
                                                        "18C_Rep1")))) |>
  summarize(n = sum(nCount_RNA)) |>
  ggplot(aes(x = n/10^6,
             y = sample_name,
             label = scales::comma_format()(n))) +
  geom_col(color = "grey20",
           fill = "grey80",
           width = 0.8)+
  scale_x_continuous(expand = c(0, 0),
                     name = "Number of RNA reads",
                     labels = 
                       scales::label_number(scale = 1, suffix = "M")) +
  geom_label(nudge_x = -1,
             label.size = 0.1) +
  scale_y_discrete(name = element_blank(),
                   labels = c("25°C Rep 2",
                              "25°C Rep 1",
                              "18°C Rep 2",
                              "18°C Rep 1")) +
  cowplot::theme_minimal_vgrid() 

# Per sample metrics
p_samples <- cowplot::plot_grid(
  p_cells_sample,
  p_reads_atac,
  p_reads_rna,
  labels = c("D", "E", "F"),
  vjust = 1,
  nrow = 3
)

# Second row of plots
p_row_2 <- cowplot::plot_grid(
  p_scatter,
  p_samples,
  labels = c("C", ""),
  rel_widths = c(1.3, 1),
  nrow = 1
)

# All plots together
cowplot::plot_grid(
  p_knees, 
  p_row_2,
  rel_heights = c(1, 1.3),
  nrow = 2
)

# Save
ggsave(here::here(fig_dir, "fig_s1.png"),
       height = 25,
       width = 25,
       units = "cm")
ggsave(here::here(fig_dir, "fig_s1.pdf"),
       height = 25,
       width = 25,
       units = "cm")


# Figure S2 --------------------------------------------------------------------

# Count ATAC
p_c_atac <- dat@meta.data |>
  ggplot(aes(x = nCount_ATAC,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
              fill = "grey90",
              width = 0.1,
              outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     labels = scales::comma_format(),
                     name = "ATAC reads per cell",
                     limits = c(0, 22000),
                     position = "top") +
  cowplot::theme_minimal_grid() 

# Feature ATAC
p_f_atac <- dat@meta.data |>
  ggplot(aes(x = nFeature_ATAC,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
               fill = "grey90",
               width = 0.1,
               outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     labels = scales::comma_format(),
                     name = "ATAC features per cell",
                     limits = c(0, 6500),
                     position = "top") +
  cowplot::theme_minimal_grid() 

# Count RNA
p_c_rna <- dat@meta.data |>
  ggplot(aes(x = nCount_RNA,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
               fill = "grey90",
               width = 0.1,
               outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     labels = scales::comma_format(),
                     name = "RNA reads per cell",
                     limits = c(0, 7000),
                     position = "top") +
  cowplot::theme_minimal_grid() 

# Feature RNA
p_f_rna <- dat@meta.data |>
  ggplot(aes(x = nFeature_RNA,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
               fill = "grey90",
               width = 0.1,
               outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     labels = scales::comma_format(),
                     name = "Expressed genes detected per cell",
                     limits = c(0, 2250),
                     position = "top") +
  cowplot::theme_minimal_grid() 

# Count and feature plots together
p_c_f <- cowplot::plot_grid(
  p_c_rna,
  p_c_atac,
  p_f_rna,
  p_f_atac,
  nrow = 2,
  labels = c("A", "C", "B", "D")
)


# Plot transcription start site enrichment for each sample
p_tss <- TSSPlot(
  dat,
  assay = "ATAC",
  group.by = "sample_name_fig") +
  scale_color_manual(values = rep("grey50", 4)) +
  cowplot::theme_minimal_grid() +
  theme(legend.position = "none",
        plot.title = element_blank(),
        strip.background = element_rect(fill = "grey90"))

# Detailed RNA QC metrics -----

# Percent mito
p_mt <- dat@meta.data |>
  ggplot(aes(x = percent.mt/100,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
               fill = "grey90",
               width = 0.1,
               outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     name = "RNA reads mapped to mitochondrial genome",
                     limits = c(0, .0525),
                     labels = scales::label_percent(),
                     position = "top") +
  cowplot::theme_minimal_grid() 

# Percent ribo
p_ribo <- dat@meta.data |>
  ggplot(aes(x = percent.ribo/100,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
               fill = "grey90",
               width = 0.1,
               outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     name = "RNA reads from ribosomal genes",
                     limits = c(0, 0.27),
                     labels = scales::label_percent(),
                     position = "top") +
  cowplot::theme_minimal_grid() 

# Count and feature plots together
p_rna_qc <- cowplot::plot_grid(
  p_mt,
  p_ribo,
  nrow = 2,
  labels = c("F", "G"),
  vjust = 1
)

# Detailed ATAC QC metrics -----

# Nucleosome signal
p_nuc <- dat@meta.data |>
  ggplot(aes(x = nucleosome_signal,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
               fill = "grey90",
               width = 0.1,
               outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     name = "Nucleosome signal",
                     limits = c(0, 1.4),
                     breaks = seq(from = 0, to = 1.25, by = 0.25),
                     position = "top") +
  cowplot::theme_minimal_grid() 

# FRiP 
p_frip <- dat@meta.data |>
  ggplot(aes(x = FRiP,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
               fill = "grey90",
               width = 0.1,
               outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     name = "Fraction of reads in peaks (FRiP)",
                     limits = c(0, 1.05),
                     breaks = seq(from = 0, to = 1, by = 0.25),
                     position = "top") +
  cowplot::theme_minimal_grid() 


# Count and feature plots together
p_atac_qc <- cowplot::plot_grid(
  p_nuc,
  p_frip,
  nrow = 2,
  labels = c("H", "I"),
  vjust = 1
)

p_qc <- cowplot::plot_grid(
  p_rna_qc,
  p_atac_qc,
  nrow = 1
)

# Combine all plots
cowplot::plot_grid(
  p_c_f,
  p_tss,
  p_qc,
  nrow = 3,
  labels = c(NA, "E", NA)
)

# Save
ggsave(here::here(fig_dir, "fig_s2.png"),
       height = 30,
       width = 25,
       units = "cm")
ggsave(here::here(fig_dir, "fig_s2.pdf"),
       height = 30,
       width = 25,
       units = "cm")

