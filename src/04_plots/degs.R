# ------------------------------------------------------------------------------
# Find differentially expressed genes between conditions
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(tidyverse)
require(Seurat)
require(Signac)

# Load plot themes 
source(here::here("src/04_plots/00_plot_themes.R"))

# Output fig dir
fig_dir <- "output/figs/degs"

# Load data
res <- readRDS(here::here("output/degs/pseudobulk_DESeq_res.rds")) |> 
  as_tibble(rownames = "gene") 
dat <- readRDS(here::here("data/processed/seurat_object/10_dat_linked.rds"))
degs <- readRDS(here::here("output/degs/degs_cell-type.rds")) |> 
  filter(p_val_adj < 0.05)

# Subset of genes to label ----
# Top pvals
res_label_padj <- res |> 
  slice_min(n = 7, padj)

# Top lfc
res_label_lfc <- res |> 
  filter(padj < 0.05) |> 
  slice_max(n = 7, abs(log2FoldChange))

# Combine both
res_label <- bind_rows(res_label_padj, res_label_lfc)

# Create plot with labels
res |>
  ggplot(aes(label = gene,
             y = -log10(padj),
             x = log2FoldChange,
             fill = padj < 0.05,
             size = padj < 0.05)) +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_hline(yintercept = -log10(0.05),
             color = "grey50",
             linetype = 2) +
  geom_point(color = "grey90",
             stroke = 0.2,
             shape = 21,
             alpha = 0.9) +
  annotate(geom = "segment",
           x = 1, 
           xend = 3, 
           y = 15, 
           yend = 15,
           linewidth = 1,
           lineend = "round",
           linejoin = "round",
           color = "#43aa8b",
           alpha = 0.5,
           arrow = arrow(length = unit(0.2, "cm"))) +
  annotate(geom = "segment",
           x = -1, 
           xend = -3, 
           y = 15, 
           yend = 15,
           linewidth = 1,
           lineend = "round",
           linejoin = "round",
           color = "#f3722c",
           alpha = 0.5,
           arrow = arrow(length = unit(0.2, "cm"))) +
  annotate(geom = "text", 
           x = 2, 
           y = 15.5, 
           alpha = 0.5,
           label = "Higher in 18°C", 
           color =  "#43aa8b") +
  annotate(geom = "text", 
           x = -2, 
           y = 15.5, 
           alpha = 0.5,
           label = "Higher in 25°C", 
           color = "#f3722c") +
  ggrepel::geom_label_repel(data = res_label,
                            color = "grey20",
                            fill = "grey90",
                            alpha = 0.8,
                            fontface = "bold",
                            size = 3,
                            max.iter = 10000,
                            min.segment.length = 0.1,
                            max.overlaps = 100) +
  scale_size_manual(values = c(1, 2)) +
  scale_fill_manual(values = c("grey80", "firebrick")) +
  scale_y_continuous(expand = c(0, 0.5),
                     breaks = c(0, 5, 10, 15),
                     labels = c("0", "5", "10", "15")) +
  scale_x_continuous(name = "18°C vs. 25°C\nlog2(fold change)") +
  cowplot::theme_minimal_grid() +
  theme(legend.position = "none") 

# Save volcano
ggsave(here::here(fig_dir, "pseudobulk_volcano.pdf"),
       height = 20,
       width = 20,
       units = "cm")
ggsave(here::here(fig_dir, "pseudobulk_volcano.png"),
       height = 20,
       width = 20,
       units = "cm")



# Number of degs per cell type
degs |> 
  group_by(cell_type) |> 
  tally() |> 
  full_join(color_cell_type) |> 
  mutate(colors = colorspace::desaturate(colors, amount = 0.6)) |> 
  mutate(cell_type = fct_reorder(cell_type, n, .fun = "identity")) |>
  ggplot() +
  geom_col(aes(x = n,
               y = cell_type,
               fill = colors),
           na.rm = TRUE,
           color = "grey20",
           alpha = 0.5) +
  labs(y = "") +
  scale_x_continuous(position = "top",
                     name = "Number of differentially expressed genes",
                     expand = c(0, 0)) +
  scale_fill_identity() +
  cowplot::theme_minimal_vgrid()

ggsave(here::here(fig_dir, "degs_count.pdf"),
       height = 15,
       width = 30,
       units = "cm")
ggsave(here::here(fig_dir, "degs_count.png"),
       height = 15,
       width = 30,
       units = "cm")

n_genes <- readRDS(here::here(out_dir, "degs_cell-type_lfc_0_min.pct_0.02.rds")) |> 
  group_by(cell_type) |> 
  tally(name = "n_genes")

readRDS(here::here("output/degs/degs_cell-type.rds")) |> 
  filter(p_val_adj < 0.05) |> 
  group_by(cell_type) |> 
  tally(name = "n_degs") |> 
  full_join(n_genes) |> 
  mutate(proportion_degs = n_degs/n_genes) |> 
  select(cell_type, proportion_degs) |> 
  distinct(cell_type, .keep_all = TRUE) |> 
  ungroup() |> 
  full_join(color_cell_type) |> 
  mutate(colors = colorspace::desaturate(colors, amount = 0.6)) |> 
  mutate(cell_type = fct_reorder(cell_type, proportion_degs)) |>
  ggplot() +
  geom_col(aes(x = proportion_degs,
               y = cell_type,
               fill = colors),
           na.rm = TRUE,
           color = "grey20") +
  labs(y = "") +
  scale_x_continuous(position = "top",
                     name = "Percent differentially expressed genes",
                     labels = scales::percent_format(accuracy = 0.01),
                     expand = c(0, 0)) +
  scale_fill_identity() +
  cowplot::theme_minimal_vgrid()

ggsave(here::here(fig_dir, "degs_percent.pdf"),
       height = 15,
       width = 30,
       units = "cm")
ggsave(here::here(fig_dir, "degs_percent.png"),
       height = 15,
       width = 30,
       units = "cm")

# Set default assay
DefaultAssay(dat) <- "RNA"

# Reset primary idents to the acclimation temperatures
Idents(dat) <- "acc_temp"

# Average expression of RNA across all cells
avg_exp <- log10(AverageExpression(dat, verbose = FALSE)$RNA)
avg_exp <- avg_exp %>%
  as_tibble(rownames = "gene")

# Scatter plot of the average expression for all cells -------------------------
p <- avg_exp %>%
  ggplot() +
  geom_point(aes(x = `25°C`,
                 y = `18°C`,
                 label = gene)) +
  labs(x = "18°C",
       y = "25°C") +
  theme_minimal()

ggsave(here::here(fig_dir, "scatter.pdf"),
       plot = p,
       width = 20,
       height = 20,
       units = "cm")
ggsave(here::here(fig_dir, "scatter.png"),
       plot = p,
       width = 20,
       height = 20,
       units = "cm")

# Quick plotly with gene labels
plotly::ggplotly(p)

degs |>
  ggplot(aes(label = gene)) +
  geom_hline(yintercept = -log10(0.05),
             color = "grey50",
             linetype = 2) +
  geom_point(aes(y = -log10(padj),
                 x = avg_log2FC_18_25,
                 fill = padj < 0.05,
                 size = padj < 0.05),
             color = "grey90",
             stroke = 0.2,
             shape = 21) +
  scale_size_manual(values = c(1, 2)) +
  scale_fill_manual(values = c("grey80", "firebrick")) +
  cowplot::theme_cowplot() +
  theme(legend.position = "none") +
  facet_wrap(~cell_type)

ggsave(here::here(fig_dir, "celltype_volcano.pdf"),
       height = 30,
       width = 30,
       units = "cm")
ggsave(here::here(fig_dir, "celltype_volcano.png"),
       height = 30,
       width = 30,
       units = "cm")



# Coverage plot example -----
DefaultAssay(dat) <- "ATAC"
Idents(dat) <- "cell_type"
gene <- "CG13427"

dars <- readRDS(here::here("output/dars/dars.rds"))



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

ggsave(here::here(fig_dir, gene, "coverage.png"),
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
  scale_fill_manual(values = acc_colors) +
  coord_flip() +
  cowplot::theme_cowplot()


ggsave(here::here(fig_dir, gene, "vln.png"),
       height = 25,
       width = 20,
       units = "cm")

################################################################################
# CLEAN ########################################################################
################################################################################

degs <- readRDS(here::here("output/degs/degs_cell-type.rds"))

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

degs |>  
  group_by(cell_type) |> 
  add_tally() |> 
  summarise(n = mean(n),
            avg_log2FC = mean(avg_log2FC_18_25)) |> 
  full_join(color_cell_type) |> 
  mutate(cell_type = fct_reorder(cell_type, avg_log2FC)) |> 
  ggplot(aes(y = cell_type,
             x = avg_log2FC)) +
  geom_segment(aes(xend = 0, 
                   yend = cell_type,
                   color = colors),
               linewidth = 1) +
  geom_vline(xintercept = 0, color = "grey50") +
  geom_point(aes(fill = colors, 
                 size = n), 
             color = "grey80",
             shape = 21) +
  scale_x_continuous(name = "log2(fold-change)",
                     limits = c(-0.25, 0.25)) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_size_continuous(range = c(2, 10)) +
  labs(y = element_blank()) +
  cowplot::theme_minimal_vgrid() +
  theme(panel.grid.major.y = element_line(color = "grey95"),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

# Save
ggsave(here::here(fig_dir, "exp_lollipop.png"),
       height = 12,
       width = 20,
       units = "cm")

degs |> 
  filter(p_val_adj < 0.05) |> 
  group_by(cell_type) |> 
  add_tally() |> 
  summarise(n = mean(n),
            avg_log2FC = mean(avg_log2FC_18_25)) |> 
  full_join(color_cell_type) |> 
  mutate(cell_type = fct_reorder(cell_type, avg_log2FC)) |> 
  ggplot(aes(y = cell_type,
             x = avg_log2FC)) +
  geom_segment(aes(xend = 0, 
                   yend = cell_type,
                   color = colors),
               linewidth = 1) +
  geom_vline(xintercept = 0, color = "grey50") +
  geom_point(aes(fill = colors, 
               size = n), 
           color = "grey80",
           shape = 21) +
  scale_x_continuous(name = "log2(fold-change)",
                     limits = c(-0.6, 0.6)) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_size_continuous(range = c(2, 10),
                        breaks = c(1, 5, 10, 20)) +
  labs(y = element_blank()) +
  cowplot::theme_minimal_vgrid() +
  theme(panel.grid.major.y = element_line(color = "grey95"),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

# Save
ggsave(here::here(fig_dir, "degs_lollipop_p_val_adj.png"),
       height = 12,
       width = 20,
       units = "cm")

degs |> 
  filter(padj < 0.05) |> 
  group_by(cell_type) |> 
  add_tally() |> 
  summarise(n = mean(n),
            avg_log2FC = mean(avg_log2FC_18_25)) |> 
  full_join(color_cell_type) |> 
  mutate(cell_type = fct_reorder(cell_type, avg_log2FC)) |> 
  ggplot(aes(y = cell_type,
             x = avg_log2FC)) +
  geom_segment(aes(xend = 0, 
                   yend = cell_type,
                   color = colors),
               linewidth = 1) +
  geom_vline(xintercept = 0, color = "grey50") +
  geom_point(aes(fill = colors, 
                 size = n), 
             color = "grey80",
             shape = 21) +
  scale_x_continuous(name = "log2(fold-change)",
                     limits = c(-0.75, 0.75)) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_size_continuous(range = c(2, 10)) +
  labs(y = element_blank()) +
  cowplot::theme_minimal_vgrid() +
  theme(legend.position = "none",
        panel.grid.major.y = element_line(color = "grey95"),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

# Save
ggsave(here::here(fig_dir, "degs_lollipop_padj.png"),
       height = 12,
       width = 20,
       units = "cm")

dars <- readRDS(here::here("output/dars/dars_cell-type.rds"))


dars |> 
  group_by(cell_type) |> 
  add_tally() |> 
  summarise(n = mean(n),
            avg_log2FC = mean(avg_log2FC_18_25)) |> 
  full_join(color_cell_type) |> 
  filter(!is.na(n)) |> 
  mutate(cell_type = fct_reorder(cell_type, avg_log2FC)) |> 
  ggplot(aes(y = cell_type,
             x = avg_log2FC)) +
  geom_segment(aes(xend = 0, 
                   yend = cell_type,
                   color = colors),
               linewidth = 1) +
geom_point(aes(fill = colors, 
               size = n), 
           color = "grey80",
           shape = 21) +
  geom_vline(xintercept = 0, color = "grey50") +
  scale_x_continuous(name = "log2(fold-change)",
                     limits = c(-0.3, 0.3)) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_size_continuous(range = c(5, 10)) +
  labs(y = element_blank()) +
  cowplot::theme_minimal_vgrid() +
  theme(legend.position = "none",
        panel.grid.major.y = element_line(color = "grey95"),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

# Save
ggsave(here::here(fig_dir, "atac_lollipop.png"),
       height = 12,
       width = 20,
       units = "cm")

# atac
dars |> 
  filter(p_val_adj < 0.05) |> 
  group_by(cell_type) |> 
  add_tally() |> 
  summarise(n = mean(n),
            avg_log2FC = mean(avg_log2FC_18_25)) |> 
  full_join(color_cell_type) |> 
  filter(!is.na(n)) |> 
  mutate(cell_type = fct_reorder(cell_type, avg_log2FC)) |> 
  ggplot(aes(y = cell_type,
             x = avg_log2FC)) +
  geom_segment(aes(xend = 0, 
                   yend = cell_type,
                   color = colors),
               linewidth = 1) +
  geom_point(aes(fill = colors, 
                 size = n), 
             color = "grey80",
             shape = 21) +
  geom_vline(xintercept = 0, color = "grey50") +
  scale_x_continuous(name = "log2(fold-change)",
                     limits = c(-0.8, 0.8)) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_size_continuous(range = c(2, 10)) +
  labs(y = element_blank()) +
  cowplot::theme_minimal_vgrid() +
  theme(legend.position = "none",
        panel.grid.major.y = element_line(color = "grey95"),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

# Save
ggsave(here::here(fig_dir, "atac_lollipop_p_val_adj.png"),
       height = 12,
       width = 20,
       units = "cm")


# atac
dars |> 
  filter(padj < 0.05) |> 
  group_by(cell_type) |> 
  add_tally() |> 
  summarise(n = mean(n),
            avg_log2FC = mean(avg_log2FC_18_25)) |> 
  full_join(color_cell_type) |> 
  filter(!is.na(n)) |> 
  mutate(cell_type = fct_reorder(cell_type, avg_log2FC)) |> 
  ggplot(aes(y = cell_type,
             x = avg_log2FC)) +
  geom_segment(aes(xend = 0, 
                   yend = cell_type,
                   color = colors),
               linewidth = 1) +
  geom_point(aes(fill = colors, 
                 size = n), 
             color = "grey80",
             shape = 21) +
  geom_vline(xintercept = 0, color = "grey50") +
  scale_x_continuous(name = "log2(fold-change)",
                     limits = c(-1.25, 1.25)) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_size_continuous(range = c(2, 10)) +
  labs(y = element_blank()) +
  cowplot::theme_minimal_vgrid() +
  theme(legend.position = "none",
        panel.grid.major.y = element_line(color = "grey95"),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

# Save
ggsave(here::here(fig_dir, "atac_lollipop_padj.png"),
       height = 12,
       width = 20,
       units = "cm")

